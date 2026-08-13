"""External bulk/scRNA reference projection (`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`).

Projects this project's T-cell-state marker methodology
(`src/xenium_tcr_ecology/preprocess/program_scores.py`'s
cytotoxicity/exhaustion/proliferation gene sets,
`src/xenium_tcr_ecology/annotation/t_cell_substates.py`'s Treg markers
and CD4/CD8A lineage rule, `scripts/06_cell_type_annotation/04_resolve_t_cell_
substates.R`'s classification priority order) onto GSE103322 (Puram et
al. 2017, `data/external/GSE103322/README.md`) -- an independent,
full-transcriptome HNSCC scRNA-seq cohort with 1,237 annotated T cells
and 100% marker-gene coverage (checked before selecting this
reference).

**Documented simplification from the primary pipeline's scoring
method:** `scanpy.tl.score_genes` (Seurat `AddModuleScore`-equivalent)
subtracts an expression-matched random control gene set, which needs a
broad background gene pool -- fetching the full 23,686-gene matrix for
this purpose was judged unnecessary overhead for a directional sanity
check, so only the marker genes plus CD4/CD8A were fetched. Module
scores here are instead a simple per-cell mean z-score across each
programme's marker genes -- a standard, simplified module-score proxy,
not presented as an exact reproduction of the primary pipeline's
control-matched score (which would not be a like-for-like comparison
across two different assay platforms -- Smart-seq2 full-transcriptome
vs. a 623-gene targeted Xenium panel -- in any case).

The discrete-state classification rule applied to these scores exactly
matches `04_resolve_t_cell_substates.R`'s priority order (Treg >
Cycling > Exhausted > Cytotoxic > CD4 > CD8 > Ambiguous), not an
approximation.

**Acquisition**: GSE103322 is downloaded automatically if not already
present -- public, unauthenticated NCBI GEO FTP-over-HTTPS, no license
restriction (`data/external/GSE103322/README.md`). Verified by SHA-256
checksum after every download, real or skipped-because-already-present.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.annotation.t_cell_substates import TREG_MARKERS
from xenium_tcr_ecology.infra.download import download_file, verify_checksums
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.preprocess.program_scores import PROGRAM_GENE_SETS

MODULE_SCORE_PROGRAMS = ["cytotoxicity", "exhaustion", "proliferation"]
LINEAGE_MARKERS = ["CD4", "CD8A"]
NON_CANCER_CELL_TYPE_ROW = "non-cancer cell type"
CANCER_CLASSIFICATION_ROW = "classified  as cancer cell"
METADATA_ROW_COUNT = 5

GSE103322_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103322/suppl/"
    "GSE103322_HNSCC_all_data.txt.gz"
)
GSE103322_EXPECTED_SIZE_BYTES = 90_226_510


def ensure_gse103322_acquired(project_root: Path) -> Path:
    """Downloads `GSE103322_HNSCC_all_data.txt.gz` if not already present
    at the expected size, then verifies its SHA-256 checksum. Safe to
    call every run -- a no-op download when the file already matches."""
    dataset_dir = project_root / "data" / "external" / "GSE103322"
    dest = dataset_dir / "GSE103322_HNSCC_all_data.txt.gz"
    download_file(GSE103322_URL, dest, expected_size_bytes=GSE103322_EXPECTED_SIZE_BYTES)

    checksum_results = verify_checksums(dataset_dir)
    failed = [f for f, ok in checksum_results.items() if not ok]
    if failed:
        raise PipelineError(
            f"Checksum verification failed for {failed} in '{dataset_dir}' -- "
            "re-download required, not silently trusted."
        )
    return dest


def required_gse103322_genes() -> set[str]:
    """Pure, testable: the exact gene set this module needs from
    GSE103322 -- every program marker gene plus Treg markers plus
    CD4/CD8A."""
    genes: set[str] = set()
    for program in MODULE_SCORE_PROGRAMS:
        genes.update(PROGRAM_GENE_SETS[program])
    genes.update(TREG_MARKERS)
    genes.update(LINEAGE_MARKERS)
    return genes


def list_gse103322_genes(path: Path) -> list[str]:
    """Fast single-pass read of every gene name in GSE103322 (first
    field of each gene row only, values not parsed) -- used to build a
    candidate pool for background/null-model gene sampling
    (`12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`), not for the required marker genes themselves.
    """
    genes = []
    with gzip.open(path, "rt") as f:
        for _ in range(1 + METADATA_ROW_COUNT):
            f.readline()
        for line in f:
            gene = line.split("\t", 1)[0]
            genes.append(gene.strip("'"))
    return genes


def parse_gse103322_t_cells(path: Path, extra_genes: set[str] | None = None) -> pd.DataFrame:
    """Memory-efficient two-pass parse of the Puram et al. 2017
    supplementary matrix: reads the header, the 5 metadata rows, and
    only the specific marker-gene rows this module needs (skipping
    most other gene rows entirely), restricted to annotated T cells
    (`non-cancer cell type == "T cell"`).

    `extra_genes` (e.g. a random background gene sample for a null
    model, `12_external_checkpoint_validation/01_test_transcriptional_program_transfer.py`) are fetched in the same file pass as the
    required marker genes -- not a second, wasteful re-read of the
    90MB compressed file."""
    all_wanted = required_gse103322_genes() | (extra_genes or set())
    wanted_quoted = {f"'{g}'" for g in all_wanted}

    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]

        metadata: dict[str, list[str]] = {}
        for _ in range(METADATA_ROW_COUNT):
            fields = f.readline().rstrip("\n").split("\t")
            metadata[fields[0]] = fields[1 : len(cell_ids) + 1]

        gene_rows: dict[str, np.ndarray] = {}
        for line in f:
            fields = line.rstrip("\n").split("\t", 1)
            gene = fields[0]
            if gene in wanted_quoted:
                values = fields[1].split("\t")[: len(cell_ids)]
                gene_rows[gene.strip("'")] = np.array(values, dtype=float)
                if len(gene_rows) == len(wanted_quoted):
                    break

    missing = required_gse103322_genes() - set(gene_rows)
    if missing:
        raise PipelineError(f"Required marker gene(s) not found in GSE103322: {sorted(missing)}")
    if NON_CANCER_CELL_TYPE_ROW not in metadata:
        raise PipelineError(
            f"Expected metadata row '{NON_CANCER_CELL_TYPE_ROW}' not found in GSE103322 -- file format may have changed."
        )

    metadata_df = pd.DataFrame(metadata, index=cell_ids)
    expr_df = pd.DataFrame(gene_rows, index=cell_ids)

    is_t_cell = metadata_df[NON_CANCER_CELL_TYPE_ROW] == "T cell"
    return expr_df.loc[is_t_cell].join(metadata_df.loc[is_t_cell])


def parse_gse103322_full_matrix_t_cells(path: Path) -> pd.DataFrame:
    """Full-transcriptome parse of GSE103322 (all 23,686 gene rows, not
    just the marker-gene subset `parse_gse103322_t_cells` reads),
    restricted to annotated T cells (`non-cancer cell type == "T
    cell"`).

    Added for `12_external_checkpoint_validation/05_rescore_cycling_
    state_with_primary_method.py`: `scanpy.tl.score_genes` needs a
    broad background gene pool for expression-matched control-gene
    sampling, which the original marker-only parse deliberately did not
    fetch (see this module's docstring) -- this function provides it.
    The memory cost is modest (1,237 T cells x 23,686 genes, float32,
    ~117MB), not the overhead the original simplification was avoiding
    for a lighter directional check; a proper re-score justifies the
    cost."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]

        metadata: dict[str, list[str]] = {}
        for _ in range(METADATA_ROW_COUNT):
            fields = f.readline().rstrip("\n").split("\t")
            metadata[fields[0]] = fields[1 : len(cell_ids) + 1]

        gene_names: list[str] = []
        rows: list[np.ndarray] = []
        for line in f:
            fields = line.rstrip("\n").split("\t", 1)
            gene_names.append(fields[0].strip("'"))
            values = fields[1].split("\t")[: len(cell_ids)]
            rows.append(np.array(values, dtype=np.float32))

    if NON_CANCER_CELL_TYPE_ROW not in metadata:
        raise PipelineError(
            f"Expected metadata row '{NON_CANCER_CELL_TYPE_ROW}' not found in GSE103322 -- file format may have changed."
        )

    metadata_df = pd.DataFrame(metadata, index=cell_ids)
    is_t_cell = metadata_df[NON_CANCER_CELL_TYPE_ROW] == "T cell"

    expr_matrix = np.stack(rows, axis=1)  # cells x genes
    expr_df = pd.DataFrame(expr_matrix, index=cell_ids, columns=gene_names)

    return expr_df.loc[is_t_cell].join(metadata_df.loc[is_t_cell, [NON_CANCER_CELL_TYPE_ROW]])


def compute_module_score(expr_df: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Pure, testable: mean z-scored expression across `genes` -- a
    simplified module-score proxy (module docstring explains why this
    differs from `scanpy.tl.score_genes`)."""
    sub = expr_df[genes]
    std = sub.std(ddof=0)
    z = (sub - sub.mean()) / std.replace(0, np.nan)
    return z.mean(axis=1)


def assign_t_cell_substate(
    treg_score: np.ndarray,
    cycling_score: np.ndarray,
    exhaustion_score: np.ndarray,
    cytotoxicity_score: np.ndarray,
    cd4_expr: np.ndarray,
    cd8a_expr: np.ndarray,
) -> np.ndarray:
    """Pure, testable: exact priority-order match to `scripts/06_
    annotation/04_resolve_t_cell_substates.R`'s `assign_t_cell_substate`
    (Treg > Cycling > Exhausted > Cytotoxic > CD4 > CD8 > Ambiguous)."""
    n = len(treg_score)
    state = np.full(n, "Ambiguous", dtype=object)

    is_cd8 = (cd8a_expr > cd4_expr) & (cd8a_expr > 0)
    is_cd4 = (cd4_expr > cd8a_expr) & (cd4_expr > 0)
    is_cytotoxic = cytotoxicity_score > 0
    is_exhausted = exhaustion_score > 0
    is_cycling = cycling_score > 0
    is_treg = (treg_score > 0) & (cd4_expr > cd8a_expr)

    state[is_cd8] = "CD8"
    state[is_cd4] = "CD4"
    state[is_cytotoxic] = "Cytotoxic"
    state[is_exhausted] = "Exhausted"
    state[is_cycling] = "Cycling"
    state[is_treg] = "Treg"
    return state


def compare_state_proportions(
    reference_states: pd.Series, project_states: pd.Series
) -> pd.DataFrame:
    """Pure, testable: side-by-side state-proportion comparison table
    between the external reference and this project's data."""
    ref_counts = reference_states.value_counts(normalize=True)
    proj_counts = project_states.value_counts(normalize=True)
    all_states = sorted(set(ref_counts.index) | set(proj_counts.index))
    return pd.DataFrame(
        {
            "state": all_states,
            "reference_fraction": [float(ref_counts.get(s, 0.0)) for s in all_states],
            "project_fraction": [float(proj_counts.get(s, 0.0)) for s in all_states],
        }
    )


def build_bulk_projection(project_root: Path) -> dict:
    gse103322_path = ensure_gse103322_acquired(project_root)
    t_cell_states_path = project_root / "data" / "derived" / "t_cell_states.parquet"
    output_path = project_root / "data" / "derived" / "bulk_reference_projection.parquet"

    if not t_cell_states_path.exists():
        raise PipelineError(
            f"'{t_cell_states_path}' not found. Run "
            "`06_cell_type_annotation/04_resolve_t_cell_substates.R` first."
        )

    reference = parse_gse103322_t_cells(gse103322_path)

    reference_scores = {}
    for program in MODULE_SCORE_PROGRAMS:
        reference_scores[program] = compute_module_score(reference, PROGRAM_GENE_SETS[program])
    reference_scores["treg"] = compute_module_score(reference, TREG_MARKERS)

    reference["state"] = assign_t_cell_substate(
        treg_score=reference_scores["treg"].to_numpy(),
        cycling_score=reference_scores["proliferation"].to_numpy(),
        exhaustion_score=reference_scores["exhaustion"].to_numpy(),
        cytotoxicity_score=reference_scores["cytotoxicity"].to_numpy(),
        cd4_expr=reference["CD4"].to_numpy(),
        cd8a_expr=reference["CD8A"].to_numpy(),
    )

    project_states = pd.read_parquet(t_cell_states_path)["t_cell_state"]

    comparison = compare_state_proportions(reference["state"], project_states)

    result = reference[["state"]].copy()
    result.index.name = "reference_cell_id"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.reset_index().to_parquet(output_path)

    comparison_path = project_root / "data" / "derived" / "bulk_reference_state_comparison.parquet"
    comparison.to_parquet(comparison_path)

    return {
        "n_reference_t_cells": len(reference),
        "n_project_t_cells": len(project_states),
        "reference_state_counts": reference["state"].value_counts().to_dict(),
        "comparison_path": str(comparison_path),
        "output_path": str(output_path),
    }
