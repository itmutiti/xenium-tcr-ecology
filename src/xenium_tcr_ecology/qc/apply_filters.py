"""Apply QC filters with a full per-cell audit trail (`04_quality_control/07_apply_qc_filters_with_audit_trail.py`).

Applies `04_quality_control/06_define_qc_thresholds_hierarchically.R`'s active threshold profile (`config/qc_thresholds.yaml`)
plus `04_quality_control/02_detect_spatial_qc_artifacts.py`'s FOV-artifact flags to every cell, records a reason code
for every excluded cell (not just a pass/fail bit), and writes a physically
filtered analysis object -- while the untouched `03_spatialdata_import/05_build_combined_analysis_object.py` combined object
(`data/objects/hnscc_xenium_combined.h5ad`) remains on disk as the
immutable, reconstructable source of truth, so no raw data is destroyed by
this step even though the derived `qc_filtered.h5ad` is a smaller subset.

Transcript-level QV filtering was considered and explicitly not applied to
the primary count matrix -- this module only filters
cells, consistent with that decision.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd
import spatialdata as sd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.spatial_artifacts import _modified_z_scores

REQUIRED_CELL_QC_COLUMNS = [
    "section_id",
    "patient_id",
    "transcript_counts",
    "n_genes_detected",
    "control_probe_ratio",
    "control_codeword_ratio",
    "nucleus_area",
]

# Every criterion that can contribute to excluding a cell. Order here is
# also the order reason codes are reported in, for readability.
EXCLUSION_REASON_COLUMNS = [
    "low_transcript_count",
    "low_genes_detected",
    "high_control_probe_ratio",
    "high_control_codeword_ratio",
    "section_relative_low_count_outlier",
    "fov_artifact_candidate",
]


def assign_fov_to_cells(zarr_path: Path) -> pd.Series:
    """Majority-vote FOV membership per cell, derived from the transcripts
    table. `cells.parquet` carries no `fov_name` field directly (confirmed
    against the data), so FOV membership must be inferred from the
    per-molecule table. Confirmed against the data that ~8.9% of cells
    have transcripts split across more than one FOV (tissue-tile boundary
    effects) -- a majority vote (not simply the first transcript's FOV) is
    required to resolve this correctly, not just conveniently.
    """
    sdata = sd.read_zarr(zarr_path)
    pts = sdata["transcripts"][["cell_id", "fov_name"]]
    pair_counts = pts.groupby(["cell_id", "fov_name"], observed=True).size().compute()
    top_fov = pair_counts.groupby(level=0, observed=True).idxmax().map(lambda pair: pair[1])
    return top_fov.drop(index="UNASSIGNED", errors="ignore")


def _local_cell_id(index: pd.Index, section_ids: pd.Series) -> pd.Series:
    """Recovers the per-section raw cell_id from the combined object's
    globally-unique obs_names ("{section_id}_{local_cell_id}", confirmed
    against the data) -- needed to join against assign_fov_to_cells()'s
    per-section, locally-indexed output.
    """
    return pd.Series(
        [idx[len(sec) + 1 :] for idx, sec in zip(index, section_ids)],
        index=index,
    )


def evaluate_thresholds(cell_qc: pd.DataFrame, profile: dict) -> pd.DataFrame:
    """Pure function: given cell QC metrics (with a 'z_counts' column
    already attached) and one threshold profile, returns a boolean
    reason-code DataFrame -- factored out from the exclusion-log builder so
    the core exclusion logic is directly testable without a
    SpatialData store or config file on disk.
    """
    reasons = pd.DataFrame(index=cell_qc.index)
    reasons["low_transcript_count"] = (
        cell_qc["transcript_counts"] < profile["min_transcript_counts"]
    )
    reasons["low_genes_detected"] = cell_qc["n_genes_detected"] < profile["min_genes_detected"]
    reasons["high_control_probe_ratio"] = (
        cell_qc["control_probe_ratio"] > profile["max_control_probe_ratio"]
    )
    reasons["high_control_codeword_ratio"] = (
        cell_qc["control_codeword_ratio"] > profile["max_control_codeword_ratio"]
    )
    reasons["section_relative_low_count_outlier"] = (
        cell_qc["z_counts"] < profile["section_relative_min_counts_z"]
    )
    return reasons


def build_exclusion_log(
    cell_qc_metrics_path: Path,
    thresholds_yaml_path: Path,
    fov_artifact_path: Path,
    spatialdata_root: Path,
) -> pd.DataFrame:
    if not cell_qc_metrics_path.is_file():
        raise PipelineError(
            f"'{cell_qc_metrics_path}' not found. Run `04_quality_control/00_compute_cell_level_qc_metrics.py` first."
        )
    if not thresholds_yaml_path.is_file():
        raise PipelineError(
            f"'{thresholds_yaml_path}' not found. Run `04_quality_control/06_define_qc_thresholds_hierarchically.R` first."
        )
    if not fov_artifact_path.is_file():
        raise PipelineError(
            f"'{fov_artifact_path}' not found. Run `04_quality_control/02_detect_spatial_qc_artifacts.py` first."
        )

    cell_qc = pd.read_parquet(cell_qc_metrics_path)
    missing = [c for c in REQUIRED_CELL_QC_COLUMNS if c not in cell_qc.columns]
    if missing:
        raise PipelineError(f"'{cell_qc_metrics_path}' is missing required column(s) {missing}.")

    with open(thresholds_yaml_path) as f:
        config = yaml.safe_load(f)
    profile = config["profiles"][config["active_profile"]]

    cell_qc["z_counts"] = cell_qc.groupby("section_id", observed=True)[
        "transcript_counts"
    ].transform(_modified_z_scores)
    reasons = evaluate_thresholds(cell_qc, profile)

    fov_df = pd.read_csv(fov_artifact_path, sep="\t")
    flagged_fovs = set(
        zip(
            fov_df.loc[fov_df["flagged_artifact_candidate"], "section_id"],
            fov_df.loc[fov_df["flagged_artifact_candidate"], "fov_name"],
        )
    )

    local_cell_id = _local_cell_id(cell_qc.index, cell_qc["section_id"])
    fov_name_per_cell = pd.Series(index=cell_qc.index, dtype=object)
    for section_id in cell_qc["section_id"].unique():
        zarr_path = spatialdata_root / f"{section_id}.zarr"
        if not zarr_path.exists():
            raise PipelineError(
                f"'{zarr_path}' not found. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
            )
        cell_to_fov = assign_fov_to_cells(zarr_path)
        mask = cell_qc["section_id"] == section_id
        fov_name_per_cell.loc[mask] = local_cell_id.loc[mask].map(cell_to_fov)

    reasons["fov_artifact_candidate"] = [
        (sec, fov) in flagged_fovs for sec, fov in zip(cell_qc["section_id"], fov_name_per_cell)
    ]

    excluded = reasons[EXCLUSION_REASON_COLUMNS].any(axis=1)
    exclusion_reason_codes = reasons[EXCLUSION_REASON_COLUMNS].apply(
        lambda row: ";".join(c for c in EXCLUSION_REASON_COLUMNS if row[c]), axis=1
    )

    log = pd.DataFrame(
        {
            "section_id": cell_qc["section_id"],
            "patient_id": cell_qc["patient_id"],
            "qc_pass": ~excluded,
            "excluded": excluded,
            "exclusion_reason_codes": exclusion_reason_codes,
            # Informational only, per `04_quality_control/06_define_qc_thresholds_hierarchically.R`'s config -- a missing
            # nucleus segmentation does not by itself indicate the
            # transcript signal for that cell is unreliable, so this does
            # not feed into `excluded` above.
            "no_nucleus_detected_flag": cell_qc["nucleus_area"].isna(),
            "fov_name": fov_name_per_cell,
        },
        index=cell_qc.index,
    )
    return log


def apply_qc_filters(
    combined_h5ad_path: Path, exclusion_log: pd.DataFrame, output_h5ad_path: Path
) -> ad.AnnData:
    if not combined_h5ad_path.is_file():
        raise PipelineError(
            f"'{combined_h5ad_path}' not found. Run `03_spatialdata_import/05_build_combined_analysis_object.py` first."
        )

    adata = ad.read_h5ad(combined_h5ad_path)
    missing = adata.obs_names.difference(exclusion_log.index)
    if len(missing) > 0:
        raise PipelineError(
            f"{len(missing)} cell(s) in '{combined_h5ad_path}' have no entry in the "
            "exclusion log -- the log must cover every cell in the combined object."
        )

    keep = exclusion_log.loc[adata.obs_names, "qc_pass"]

    # Quality Control completion gate: "no entire patient removed without a
    # documented scientific decision." A patient losing *all* its cells to
    # the numeric threshold profile would be exactly that, happening
    # silently -- surfaced as a hard error requiring human review, not
    # auto-applied.
    n_before = adata.obs.groupby("patient_id", observed=True).size()
    n_after = adata.obs.loc[keep.to_numpy()].groupby("patient_id", observed=True).size()
    n_after = n_after.reindex(n_before.index, fill_value=0)
    fully_excluded_patients = n_before.index[n_after == 0].tolist()
    if fully_excluded_patients:
        raise PipelineError(
            f"Patient(s) {fully_excluded_patients} would be entirely excluded by the "
            "active QC threshold profile. Per the Quality Control completion gate, no entire "
            "patient may be removed without an explicit, documented scientific "
            "decision -- this requires human review, not automatic exclusion."
        )

    filtered = adata[keep.to_numpy()].copy()
    output_h5ad_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(output_h5ad_path)
    return filtered


def build_qc_filter_report(project_root: Path) -> dict:
    cell_qc_metrics_path = project_root / "data" / "derived" / "cell_qc_metrics.parquet"
    thresholds_yaml_path = project_root / "config" / "qc_thresholds.yaml"
    fov_artifact_path = (
        project_root / "reports" / "qc" / "spatial_artifact_masks" / "fov_artifact_candidates.tsv"
    )
    spatialdata_root = project_root / "data" / "objects" / "spatialdata"
    combined_h5ad_path = project_root / "data" / "objects" / "hnscc_xenium_combined.h5ad"
    exclusion_log_path = project_root / "data" / "derived" / "exclusion_log.tsv"
    output_h5ad_path = project_root / "data" / "objects" / "qc_filtered.h5ad"

    exclusion_log = build_exclusion_log(
        cell_qc_metrics_path, thresholds_yaml_path, fov_artifact_path, spatialdata_root
    )
    exclusion_log_path.parent.mkdir(parents=True, exist_ok=True)
    exclusion_log.to_csv(exclusion_log_path, sep="\t", index_label="cell_id")

    apply_qc_filters(combined_h5ad_path, exclusion_log, output_h5ad_path)

    n_total = len(exclusion_log)
    n_excluded = int(exclusion_log["excluded"].sum())
    reason_counts = {
        col: int(exclusion_log["exclusion_reason_codes"].str.contains(col, regex=False).sum())
        for col in EXCLUSION_REASON_COLUMNS
    }

    with open(thresholds_yaml_path) as f:
        active_profile = yaml.safe_load(f)["active_profile"]

    return {
        "n_cells_total": n_total,
        "n_cells_excluded": n_excluded,
        "fraction_excluded": round(n_excluded / n_total, 4),
        "n_cells_retained": n_total - n_excluded,
        "active_profile": active_profile,
        "reason_counts": reason_counts,
        "transcript_level_qv_filtering": (
            "Not applied to the primary count matrix -- documented, human-approved "
            "scope decision."
        ),
    }
