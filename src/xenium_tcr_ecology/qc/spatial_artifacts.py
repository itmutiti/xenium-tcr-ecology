"""Spatial QC artefact detection (`04_quality_control/02_detect_spatial_qc_artifacts.py`).

Xenium images are acquired as a mosaic of fields of view (FOVs, confirmed
against the data: 27 FOVs per section, named by a grid coordinate like
'M6', 'L7' -- a genuine tiling structure, not assumed). Two classic Xenium
artefact types are FOV-level: local decoding failures (one FOV's chemistry/
imaging degrades) and striping (systematic differences at FOV boundaries).
Both manifest as one or a few FOVs looking statistically different from
their neighbours on the same section, which is what this module tests for,
per FOV per section, using a median-absolute-deviation (MAD) based outlier
score -- robust to the fact that a few bad FOVs would otherwise
skew a standard mean/SD-based z-score.

This flags candidate artefact FOVs for review; it does not exclude them --
exclusion with a documented rationale is `04_quality_control/06_define_qc_thresholds_hierarchically.R`, `04_quality_control/07_apply_qc_filters_with_audit_trail.py`'s job.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter

REPORT_FIELDS = [
    "section_id",
    "fov_name",
    "n_transcripts",
    "mean_qv",
    "fraction_qv_below_20",
    "negative_control_fraction",
    "mad_score_mean_qv",
    "mad_score_neg_control_fraction",
    "flagged_artifact_candidate",
]

MAD_FLAG_THRESHOLD = 3.5  # a commonly used robust-outlier cutoff (Iglewicz & Hoaglin, 1993)


def _modified_z_scores(values: pd.Series) -> pd.Series:
    median = values.median()
    mad = (values - median).abs().median()
    if mad == 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return 0.6745 * (values - median) / mad


def compute_fov_qc(zarr_path: Path) -> pd.DataFrame:
    sdata = sd.read_zarr(zarr_path)
    pts = sdata["transcripts"][["fov_name", "qv", "feature_name"]].compute()

    is_neg_control = pts["feature_name"].str.startswith(("NegControlProbe_", "NegControlCodeword_"))

    grouped = pts.groupby("fov_name", observed=True).agg(
        n_transcripts=("qv", "size"),
        mean_qv=("qv", "mean"),
        fraction_qv_below_20=("qv", lambda s: (s < 20).mean()),
    )
    grouped["negative_control_fraction"] = is_neg_control.groupby(
        pts["fov_name"], observed=True
    ).mean()

    grouped["mad_score_mean_qv"] = _modified_z_scores(grouped["mean_qv"])
    grouped["mad_score_neg_control_fraction"] = _modified_z_scores(
        grouped["negative_control_fraction"]
    )
    grouped["flagged_artifact_candidate"] = (
        grouped["mad_score_mean_qv"].abs() > MAD_FLAG_THRESHOLD
    ) | (grouped["mad_score_neg_control_fraction"].abs() > MAD_FLAG_THRESHOLD)

    return grouped.reset_index()


def build_spatial_artifact_report(
    spatialdata_root: Path, output_dir: Path, project_root: Path
) -> dict:
    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(
            f"No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "fov_artifact_candidates.tsv"
    if output_path.exists():
        output_path.unlink()
    writer = InventoryWriter(
        output_path, project_root=project_root, fields=REPORT_FIELDS, delimiter="\t"
    )

    flagged_by_section: dict[str, int] = {}
    total_fovs = 0
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        fov_df = compute_fov_qc(zarr_path)
        total_fovs += len(fov_df)
        flagged_by_section[section_id] = int(fov_df["flagged_artifact_candidate"].sum())

        for _, row in fov_df.iterrows():
            writer.write_row(
                section_id=section_id,
                fov_name=row["fov_name"],
                n_transcripts=int(row["n_transcripts"]),
                mean_qv=round(float(row["mean_qv"]), 3),
                fraction_qv_below_20=round(float(row["fraction_qv_below_20"]), 4),
                negative_control_fraction=round(float(row["negative_control_fraction"]), 4),
                mad_score_mean_qv=round(float(row["mad_score_mean_qv"]), 2),
                mad_score_neg_control_fraction=round(
                    float(row["mad_score_neg_control_fraction"]), 2
                ),
                flagged_artifact_candidate=bool(row["flagged_artifact_candidate"]),
            )

    return {
        "sections_processed": len(zarr_paths),
        "total_fovs": total_fovs,
        "total_flagged_fovs": sum(flagged_by_section.values()),
        "sections_with_flagged_fovs": sum(1 for v in flagged_by_section.values() if v > 0),
    }
