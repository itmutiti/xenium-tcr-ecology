"""Transcript-level QC metric computation (`04_quality_control/01_compute_transcript_level_qc_metrics.py`).

Reads from each section's `transcripts` Points element (the SpatialData
object built in `03_spatialdata_import/01_import_each_section_to_spatialdata.py`), not from `table` -- the cell x gene
matrix excludes negative-control/unassigned-codeword features, but the
transcripts element is not gene-expression-filtered and carries every
decoded transcript regardless of type, which is what this phase actually
needs (Q-values, unassigned transcripts, and control-probe behaviour are
all transcript-level, not cell-level, questions).

Uses dask's lazy aggregation (groupby/mean over the Points dask DataFrame)
rather than materialising all ~6.8M transcripts per section into memory at
once as a pandas DataFrame.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Q20 is 10x's own documented minimum-confidence threshold for a decoded
# transcript (99% decode accuracy) -- used here as a QC profiling reference
# point, not (yet) as an exclusion threshold; thresholding is `04_quality_control/06_define_qc_thresholds_hierarchically.R`'s
# job, with its own documented rationale and sensitivity analysis.
QV_REFERENCE_THRESHOLD = 20.0

# Feature-type classification by name pattern -- matches the categories
# (verified against the h5's feature_types field for Gene Expression;
# negative-control/codeword naming conventions confirmed against the same
# per-section h5 files' var index for the other three categories).
_NEG_PROBE_PREFIX = "NegControlProbe_"
_NEG_CODEWORD_PREFIX = "NegControlCodeword_"
_UNASSIGNED_CODEWORD_PREFIX = "UnassignedCodeword_"

# A transcript not assigned to any segmented cell has cell_id == "UNASSIGNED"
# (a literal sentinel string, confirmed against the data), not an empty
# string. An earlier version of this module checked `cell_id != ""`, which
# does not match "UNASSIGNED" and therefore silently counted every
# unassigned transcript as assigned (n_unassigned was 0 for all 18 sections
# in the resulting report -- an implausible value that should have been
# investigated at the time and was not; caught and fixed while implementing
# `04_quality_control/07_apply_qc_filters_with_audit_trail.py`, section 6).
_UNASSIGNED_CELL_ID_SENTINELS = ("", "UNASSIGNED")


def _is_assigned_to_cell(cell_id: str) -> bool:
    """Pure, testable membership check -- factored out so the sentinel
    handling can be unit-tested with plain strings, not only via a full
    SpatialData store."""
    return cell_id not in _UNASSIGNED_CELL_ID_SENTINELS


def _classify_feature_name(name: str) -> str:
    if name.startswith(_NEG_PROBE_PREFIX):
        return "negative_control_probe"
    if name.startswith(_NEG_CODEWORD_PREFIX):
        return "negative_control_codeword"
    if name.startswith(_UNASSIGNED_CODEWORD_PREFIX):
        return "unassigned_codeword"
    return "gene_expression"


def compute_transcript_qc_metrics(zarr_path: Path) -> dict:
    if not zarr_path.exists():
        raise PipelineError(
            f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    sdata = sd.read_zarr(zarr_path)
    pts = sdata["transcripts"][["feature_name", "qv", "overlaps_nucleus", "cell_id"]]

    n_total = int(pts.shape[0].compute())
    n_assigned = (
        int(pts["cell_id"].map(_is_assigned_to_cell, meta=("cell_id", "bool")).sum().compute())
        if n_total
        else 0
    )
    mean_qv = float(pts["qv"].mean().compute()) if n_total else 0.0
    median_qv = float(pts["qv"].compute().median()) if n_total else 0.0
    frac_low_qv = float((pts["qv"] < QV_REFERENCE_THRESHOLD).mean().compute()) if n_total else 0.0
    frac_nucleus = float(pts["overlaps_nucleus"].mean().compute()) if n_total else 0.0

    feature_names = pts["feature_name"].compute()
    type_counts = feature_names.map(_classify_feature_name).value_counts().to_dict()

    return {
        "n_transcripts_total": n_total,
        "n_assigned_to_cell": n_assigned,
        "n_unassigned": n_total - n_assigned,
        "fraction_overlaps_nucleus": round(frac_nucleus, 4),
        "mean_qv": round(mean_qv, 3),
        "median_qv": round(median_qv, 3),
        "fraction_qv_below_20": round(frac_low_qv, 4),
        "n_gene_expression": type_counts.get("gene_expression", 0),
        "n_negative_control_probe": type_counts.get("negative_control_probe", 0),
        "n_negative_control_codeword": type_counts.get("negative_control_codeword", 0),
        "n_unassigned_codeword": type_counts.get("unassigned_codeword", 0),
    }


def build_transcript_qc_report(spatialdata_root: Path, output_path: Path) -> dict:
    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(
            f"No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    rows = []
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        metrics = compute_transcript_qc_metrics(zarr_path)
        rows.append({"section_id": section_id, **metrics})

    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output_path)

    return {
        "sections_processed": len(zarr_paths),
        "total_transcripts": int(df["n_transcripts_total"].sum()),
        "median_fraction_qv_below_20": float(df["fraction_qv_below_20"].median()),
        "median_fraction_overlaps_nucleus": float(df["fraction_overlaps_nucleus"].median()),
    }
