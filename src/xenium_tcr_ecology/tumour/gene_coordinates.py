"""Gene genomic coordinate lookup (prerequisite for `07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R`'s CNV
appendix).

CNV inference requires ordering genes by chromosomal position, which
this project's panel metadata does not carry. Ensembl's REST API
(`rest.ensembl.org`, GRCh38) is queried, confirmed reachable
before use, via its batch `/lookup/symbol/homo_sapiens`
endpoint (accepts up to `MAX_BATCH_SIZE` symbols per POST request), rather
than adding a new local genome-annotation package dependency
(biomaRt/TxDb are heavy Bioconductor installs this environment does not
otherwise need -- consistent with `environment/conda/main.yml`'s own
stated policy of adding dependencies only when the phase that needs them
is actually implemented). Results are cached to
`references/gene_genomic_coordinates.tsv` so repeat runs do not re-hit the
API.

Only the 399 `biological_gene` features (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`) are looked up: CDR3/
HPV probes are excluded on principle, not merely convenience -- CDR3
probes are patient-specific TCR-sequence-derived oligos with no single
fixed genomic locus in the conventional sense, and HPV16 genes are viral
(not present on any human chromosome at all), so neither is a meaningful
CNV-inference input regardless of panel coverage.
"""

from __future__ import annotations

import time
from pathlib import Path

import pandas as pd
import requests

from xenium_tcr_ecology.infra.exceptions import PipelineError

ENSEMBL_LOOKUP_URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
MAX_BATCH_SIZE = 300
REQUEST_TIMEOUT_S = 30
_RETRIES = 5
_RETRY_DELAY_SECONDS = 5


def _post_batch_with_retries(batch: list[str]) -> dict:
    """A cache-miss here means a real, synchronous Ensembl REST call on the
    critical path of a clean-room run -- found to have no retry logic at
    all (a bare `requests.post`, one PipelineError on any non-200), unlike
    this project's other external-dependency helper
    (`xenium_tcr_ecology.infra.download.download_file`, which retries).
    Ensembl's REST API returning a transient 503/500 or timing out under
    load is expected occasional behaviour for a public service, not a
    sign of anything wrong locally -- confirmed: both failures
    observed during a clean-room run resolved on retry with no other
    change, and basic connectivity from the same host/container was
    unaffected throughout."""
    last_exc: Exception | None = None
    for attempt in range(1, _RETRIES + 1):
        try:
            response = requests.post(
                ENSEMBL_LOOKUP_URL,
                headers={"Content-Type": "application/json"},
                json={"symbols": batch},
                timeout=REQUEST_TIMEOUT_S,
            )
            if response.status_code == 200:
                return response.json()
            last_exc = PipelineError(
                f"Ensembl REST API returned status {response.status_code} for a batch of {len(batch)} gene(s): "
                f"{response.text[:200]}"
            )
        except requests.RequestException as exc:
            last_exc = exc
        if attempt < _RETRIES:
            time.sleep(_RETRY_DELAY_SECONDS)
    raise PipelineError(
        f"Ensembl REST API failed for a batch of {len(batch)} gene(s) after {_RETRIES} attempts: {last_exc}"
    ) from last_exc


def fetch_gene_coordinates(
    gene_symbols: list[str], batch_size: int = MAX_BATCH_SIZE
) -> pd.DataFrame:
    """Queries Ensembl's REST API for chromosome/start/end/strand for each
    gene symbol. Genes not found (withdrawn symbols, non-protein-coding
    entries Ensembl does not resolve by this exact symbol, etc.) are
    dropped, not fabricated -- the returned DataFrame may have fewer rows
    than `len(gene_symbols)`, and callers must not assume completeness."""
    rows = []
    for i in range(0, len(gene_symbols), batch_size):
        batch = gene_symbols[i : i + batch_size]
        payload = _post_batch_with_retries(batch)
        for symbol in batch:
            entry = payload.get(symbol)
            if entry is None:
                continue
            rows.append(
                {
                    "gene": symbol,
                    "chromosome": str(entry["seq_region_name"]),
                    "start": int(entry["start"]),
                    "end": int(entry["end"]),
                    "strand": int(entry["strand"]),
                }
            )
    return pd.DataFrame(rows)


def build_gene_coordinate_reference(project_root: Path, force_refetch: bool = False) -> dict:
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "references" / "gene_genomic_coordinates.tsv"

    if not feature_annotation_path.is_file():
        raise PipelineError(
            f"'{feature_annotation_path}' not found. Run `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py` first."
        )

    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    gene_symbols = sorted(
        feature_annotation.loc[
            feature_annotation["feature_class"] == "biological_gene", "feature_name"
        ]
    )

    if output_path.is_file() and not force_refetch:
        cached = pd.read_csv(output_path, sep="\t")
        if set(gene_symbols).issubset(set(cached["gene"])):
            return {
                "n_genes_requested": len(gene_symbols),
                "n_genes_resolved": int(cached["gene"].isin(gene_symbols).sum()),
                "output_path": str(output_path),
                "source": "cache",
            }

    coords = fetch_gene_coordinates(gene_symbols)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    coords.to_csv(output_path, sep="\t", index=False)

    return {
        "n_genes_requested": len(gene_symbols),
        "n_genes_resolved": len(coords),
        "n_genes_unresolved": len(gene_symbols) - len(coords),
        "unresolved_genes": sorted(set(gene_symbols) - set(coords["gene"])),
        "output_path": str(output_path),
        "source": "ensembl_rest_api",
    }
