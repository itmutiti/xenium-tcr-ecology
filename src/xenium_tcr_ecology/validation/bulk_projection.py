"""Projects Niche and Ecosystem Discovery ecosystem-derived gene signatures onto the
TCGA-HNSC bulk RNA-seq cohort and cautiously tests correlation with an
independent immune-infiltration proxy (Phase 16.04) -- the
`ecosystem_signature_bulk_validation` claim (`governance/
validation_plan.tsv`, `16_external_validation_and_generalisation/00_define_validation_claims.py`).

**Signature construction, not a generic deconvolution tool:** Niche and Ecosystem Discovery's
ecosystem labels are defined by spatial neighbourhood cell-type
composition (`10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`'s blinded rule-based annotation), not
directly a gene expression signature -- bulk RNA-seq has no spatial
information at all, so an ecosystem is projected onto bulk data via the
marker genes of the lineage(s) that define that ecosystem, reusing
`06_cell_type_annotation/00_compile_marker_and_reference_registry.py`'s already-curated, panel-checked
`CELL_TYPE_MARKER_REGISTRY` directly (not a fresh, independently-chosen
gene list). `Mixed/non-specific niche` is explicitly excluded -- it has
no coherent defining lineage by its name and definition, so no
signature can be built for it; this limitation is stated directly
rather than forcing a signature anyway.

**Scoring method, matching Phase 12.00's precedent for a
cross-platform/reference comparison:** `compute_module_score` is
imported directly from `external_checkpoint.bulk_reference` (`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`) -- the same simplified
per-sample mean z-scored expression proxy already justified there for
this kind of situation (a full scanpy-style control-matched background
score is not warranted for a cautious, exploratory bulk projection),
reused rather than reimplemented.

**Cautious success criterion, not a strong confirmatory claim:** tests
whether each ecosystem signature score correlates in the
biologically-expected direction with `PTPRC` (CD45, the canonical,
independent pan-leukocyte marker -- confirmed not a member of
any of this project's ecosystem signature gene sets, avoiding
circularity) -- immune-associated ecosystems (`NK-cell/T-cell/Myeloid
niche`, `Plasma-cell/Mast-cell niche`) are expected to correlate
positively; the `Tumour niche` signature is expected to correlate
negatively (the well-known immune-exclusion/tumour-purity trade-off in
bulk RNA-seq). Reported as qualitative support (per `governance/
validation_plan.tsv`'s stated cautious success criterion), not treated
as a confirmatory statistical test.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from scipy.stats import spearmanr

from xenium_tcr_ecology.annotation.marker_registry import CELL_TYPE_MARKER_REGISTRY
from xenium_tcr_ecology.external_checkpoint.bulk_reference import compute_module_score
from xenium_tcr_ecology.infra.exceptions import PipelineError

IMMUNE_PROXY_GENE = "PTPRC"

_LINEAGE_MARKERS = {
    entry["cell_identity"]: entry["markers"]
    for entry in CELL_TYPE_MARKER_REGISTRY
    if entry["hierarchy_level"] == "major_lineage"
}

ECOSYSTEM_SIGNATURE_GENES: dict[str, list[str]] = {
    "Fibroblast niche": _LINEAGE_MARKERS["Fibroblast"],
    "Tumour niche": _LINEAGE_MARKERS["Epithelial_Tumour"],
    "NK-cell/T-cell/Myeloid niche": sorted(
        set(_LINEAGE_MARKERS["T_cell"])
        | set(_LINEAGE_MARKERS["NK_cell"])
        | set(_LINEAGE_MARKERS["Myeloid"])
    ),
    "Perivascular/Endothelial niche": sorted(
        set(_LINEAGE_MARKERS["Perivascular_SmoothMuscle"]) | set(_LINEAGE_MARKERS["Endothelial"])
    ),
    "Plasma-cell/Mast-cell niche": sorted(
        set(_LINEAGE_MARKERS["Plasma_cell"]) | set(_LINEAGE_MARKERS["Mast_cell"])
    ),
}
EXPECTED_DIRECTION: dict[str, str] = {
    "Fibroblast niche": "unspecified",
    "Tumour niche": "negative",
    "NK-cell/T-cell/Myeloid niche": "positive",
    "Perivascular/Endothelial niche": "unspecified",
    "Plasma-cell/Mast-cell niche": "positive",
}


def compute_signature_immune_correlation(
    signature_score: pd.Series, immune_proxy: pd.Series
) -> dict:
    """Pure, testable: Spearman correlation between a per-sample
    signature score and a per-sample immune-proxy expression value."""
    rho, pvalue = spearmanr(signature_score, immune_proxy)
    return {"rho": float(rho), "pvalue": float(pvalue)}


def build_bulk_projection(project_root: Path) -> dict:
    tcga_dir = project_root / "data" / "external" / "bulk" / "TCGA-HNSC"
    counts_path = tcga_dir / "TCGA-HNSC.star_counts.tsv.gz"
    probemap_path = tcga_dir / "gencode.v36.annotation.gtf.gene.probemap"
    output_path = project_root / "data" / "derived" / "bulk_ecosystem_projection_results.parquet"

    for path, phase in [
        (
            counts_path,
            "`16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`",
        ),
        (
            probemap_path,
            "`16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`",
        ),
    ]:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    probemap = pd.read_csv(probemap_path, sep="\t")
    ensembl_to_symbol = dict(zip(probemap["id"], probemap["gene"]))

    all_genes_needed = sorted(
        set(sum(ECOSYSTEM_SIGNATURE_GENES.values(), [])) | {IMMUNE_PROXY_GENE}
    )
    symbol_to_ensembl = {v: k for k, v in ensembl_to_symbol.items()}
    missing = [g for g in all_genes_needed if g not in symbol_to_ensembl]
    if missing:
        raise PipelineError(f"Gene(s) {missing} not found in '{probemap_path}'.")
    ensembl_ids_needed = [symbol_to_ensembl[g] for g in all_genes_needed]

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    counts = counts.loc[counts.index.isin(ensembl_ids_needed)]
    counts.index = counts.index.map(ensembl_to_symbol)
    expr_df = counts.T  # samples x genes orientation, matching compute_module_score's expectation

    immune_proxy = expr_df[IMMUNE_PROXY_GENE]

    rows = []
    for ecosystem_label, genes in ECOSYSTEM_SIGNATURE_GENES.items():
        present_genes = [g for g in genes if g in expr_df.columns]
        signature_score = compute_module_score(expr_df, present_genes)
        correlation = compute_signature_immune_correlation(signature_score, immune_proxy)
        expected = EXPECTED_DIRECTION[ecosystem_label]
        if expected == "unspecified":
            matches_expected = None
        else:
            observed_sign = "positive" if correlation["rho"] > 0 else "negative"
            matches_expected = observed_sign == expected
        rows.append(
            {
                "ecosystem_label": ecosystem_label,
                "n_signature_genes": len(present_genes),
                "rho_vs_immune_proxy": correlation["rho"],
                "pvalue": correlation["pvalue"],
                "expected_direction": expected,
                "matches_expected_direction": matches_expected,
            }
        )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_tcga_samples": len(expr_df),
        "n_ecosystems_tested": len(result),
        "n_matching_expected_direction": int(result["matches_expected_direction"].sum(skipna=True)),
        "n_directional_predictions_made": int(result["matches_expected_direction"].notna().sum()),
        "output_path": str(output_path),
    }
