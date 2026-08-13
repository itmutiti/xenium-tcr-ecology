"""Primary analysis matrix freeze (`05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py`).

Consolidates `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py`'s analysis_ready.h5ad (raw counts + lognorm +
pearson_residuals + detected layers, the `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`-selected primary
layer recorded in uns, and all existing per-cell technical covariates:
transcript_counts, control_probe_ratio, size_factor, etc.) with Phase
5.03's program scores (joined into obs) into a single frozen release
directory, alongside `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`, `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`, `05_preprocessing_and_normalisation/04_model_technical_covariates.R`'s supporting tables, with a
SHA256 manifest recording the exact content of every included file, fixing
byte-for-byte what a given analysis used from this release.

This step makes no new scientific decisions: `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` already selected
the primary normalisation layer per its documented evidence-based comparison, and `05_preprocessing_and_normalisation/03_calculate_program_scores.py`, `05_preprocessing_and_normalisation/04_model_technical_covariates.R`
already computed what they compute. This is a mechanical, verifiable
freeze, not a fresh analysis choice.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

RELEASE_NAME = "v1_primary_analysis"

# Supporting tables copied into the release directory if present --
# diagnostic outputs (5.02, 5.04), not required for the frozen matrix
# itself to be valid, so their absence does not block the freeze.
OPTIONAL_SUPPORTING_FILES = [
    "metadata/feature_annotation.tsv",
    "data/derived/normalisation_benchmark_summary.tsv",
    "data/derived/variance_partition_summary.tsv",
    "config/qc_thresholds.yaml",
]


def compute_file_hash(path: Path) -> str:
    sha256 = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def merge_program_scores_into_obs(adata: ad.AnnData, program_scores: pd.DataFrame) -> ad.AnnData:
    missing = adata.obs_names.difference(program_scores.index)
    if len(missing) > 0:
        raise PipelineError(
            f"{len(missing)} cell(s) in the analysis-ready object have no entry in program_scores.parquet."
        )
    score_cols = program_scores.reindex(adata.obs_names)
    overlap = set(score_cols.columns) & set(adata.obs.columns)
    if overlap:
        raise PipelineError(
            f"Program score column(s) {overlap} would collide with existing obs columns."
        )
    for col in score_cols.columns:
        adata.obs[col] = score_cols[col].to_numpy()
    return adata


def freeze_primary_analysis_matrix(project_root: Path, release_dir: Path) -> dict:
    analysis_ready_path = project_root / "data" / "objects" / "analysis_ready.h5ad"
    program_scores_path = project_root / "data" / "derived" / "program_scores.parquet"

    if not analysis_ready_path.is_file():
        raise PipelineError(
            f"'{analysis_ready_path}' not found. Run `05_preprocessing_and_normalisation/01_construct_analysis_count_layers.py` first."
        )
    if not program_scores_path.is_file():
        raise PipelineError(
            f"'{program_scores_path}' not found. Run `05_preprocessing_and_normalisation/03_calculate_program_scores.py` first."
        )

    adata = ad.read_h5ad(analysis_ready_path)
    if "primary_normalization_layer" not in adata.uns:
        raise PipelineError(
            "analysis_ready.h5ad has no uns['primary_normalization_layer'] -- run `05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R` first."
        )
    program_scores = pd.read_parquet(program_scores_path)
    adata = merge_program_scores_into_obs(adata, program_scores)

    release_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = release_dir / "primary_analysis_matrix.h5ad"
    adata.write_h5ad(matrix_path)

    copied_files = [matrix_path]
    for rel_path in OPTIONAL_SUPPORTING_FILES:
        source = project_root / rel_path
        if source.is_file():
            dest = release_dir / source.name
            shutil.copy2(source, dest)
            copied_files.append(dest)

    checksums = {f.name: compute_file_hash(f) for f in copied_files}
    manifest = {
        "release_name": RELEASE_NAME,
        "frozen_at": datetime.now(timezone.utc).isoformat(),
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "layers": sorted(adata.layers.keys()),
        "primary_normalization_layer": adata.uns["primary_normalization_layer"],
        "obs_columns": list(adata.obs.columns),
        "files": checksums,
    }
    manifest_path = release_dir / "MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=False))

    checksums_path = release_dir / "checksums.sha256"
    checksums_path.write_text("".join(f"{h}  {name}\n" for name, h in checksums.items()))

    return {
        "release_dir": str(release_dir),
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "primary_normalization_layer": adata.uns["primary_normalization_layer"],
        "n_files": len(copied_files),
        "matrix_hash": checksums[matrix_path.name],
    }
