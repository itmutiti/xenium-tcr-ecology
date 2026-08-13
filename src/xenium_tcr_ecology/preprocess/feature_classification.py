"""Feature classification (`05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py`).

Builds a complete feature dictionary covering every distinct feature name
observed across all 18 sections' raw transcripts tables -- not just the 623
genes in the analysis matrix (`adata.X` is gene-expression-only by Phase
3.01's `gex_only` design decision, so negative
control probes/codewords and unassigned codewords are invisible there and
must be recovered from the transcripts element, matching `04_quality_control/01_compute_transcript_level_qc_metrics.py`'s own
approach).

Classifies each feature into an explicit class using rules confirmed
against the data, not assumed:
  - negative_control_probe / negative_control_codeword / unassigned_codeword:
    reuses `04_quality_control/01_compute_transcript_level_qc_metrics.py`'s `_classify_feature_name` prefix rules -- these are a
    small, fixed technical set (20 / 41 / 3 features respectively),
    confirmed identical across all 18 sections (i.e. panel-design elements,
    not patient-specific).
  - tcr_cdr3_probe: matches the CDR3 probe naming convention (a date
    prefix, 6 digits with an optional single batch letter -- e.g.
    "231004B_CASRDSPSTDTQYF_TRB" -- confirmed against all 11 distinct
    date/batch prefixes in the panel; an earlier, narrower version of
    this pattern used in `04_quality_control/08_assess_replicate_concordance.R` missed the "231004B" batch entirely).
  - hpv_probe: the 8 "HPV16_*" features, confirmed patient-specific (only
    present in 6 of 18 sections, i.e. HPV+ patients' panels) -- the only
    non-CDR3 genes in the panel that are not present in all 18 sections
    (399 "core" genes are; 216 CDR3 + 8 HPV = 224 "variable" genes are not).
  - biological_gene: everything else.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import spatialdata as sd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.qc.transcript_metrics import _classify_feature_name

HPV_PROBE_PREFIX = "HPV16_"
# Confirmed against all 11 date/batch prefixes in the combined panel
# (see module docstring); the optional single uppercase letter after the
# 6-digit date is required to match the full probe set (216/216), not
# the narrower 199/216 an earlier version of this pattern matched.
CDR3_PROBE_PATTERN = re.compile(r"^[0-9]{6}[A-Z]?_[A-Z]+_TR[AB]$")

FEATURE_CLASSES_IN_ANALYSIS_MATRIX = {"biological_gene", "hpv_probe", "tcr_cdr3_probe"}


def classify_feature(name: str) -> str:
    """Pure, testable classification for one feature name."""
    coarse = _classify_feature_name(name)
    if coarse != "gene_expression":
        return coarse
    if CDR3_PROBE_PATTERN.match(name):
        return "tcr_cdr3_probe"
    if name.startswith(HPV_PROBE_PREFIX):
        return "hpv_probe"
    return "biological_gene"


def collect_feature_presence(spatialdata_root: Path) -> pd.DataFrame:
    """Every distinct feature name observed across all sections, with the
    set of sections it appears in -- the complete feature universe."""
    zarr_paths = sorted(spatialdata_root.glob("*.zarr"))
    if not zarr_paths:
        raise PipelineError(
            f"No .zarr stores found under '{spatialdata_root}'. Run `03_spatialdata_import/01_import_each_section_to_spatialdata.py` first."
        )

    presence: dict[str, set[str]] = {}
    for zarr_path in zarr_paths:
        section_id = zarr_path.stem
        sdata = sd.read_zarr(zarr_path)
        names = sdata["transcripts"]["feature_name"].unique().compute()
        for name in names:
            presence.setdefault(name, set()).add(section_id)

    rows = [
        {
            "feature_name": name,
            "n_sections_present": len(sections),
            "sections_present": ";".join(sorted(sections)),
        }
        for name, sections in presence.items()
    ]
    return pd.DataFrame(rows).sort_values("feature_name").reset_index(drop=True)


def build_feature_annotation_report(spatialdata_root: Path, output_path: Path) -> dict:
    df = collect_feature_presence(spatialdata_root)
    n_sections_total = len(sorted(spatialdata_root.glob("*.zarr")))

    df["feature_class"] = df["feature_name"].map(classify_feature)
    df["in_analysis_matrix"] = df["feature_class"].isin(FEATURE_CLASSES_IN_ANALYSIS_MATRIX)
    df["patient_specific"] = df["n_sections_present"] < n_sections_total

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)

    class_counts = df["feature_class"].value_counts().to_dict()
    return {
        "n_features_total": len(df),
        "n_sections": n_sections_total,
        "class_counts": class_counts,
        "n_patient_specific_features": int(df["patient_specific"].sum()),
    }
