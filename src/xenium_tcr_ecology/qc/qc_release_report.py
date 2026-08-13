"""QC release report and go/no-go decision (`04_quality_control/09_generate_qc_release_report.py`).

Aggregates `04_quality_control/00_compute_cell_level_qc_metrics.py`-4.03 and 4.06-4.08's outputs into one section-by-
section report and a single structured go/no-go decision -- the Quality Control
completion gate the blueprint calls for before biological analysis begins.

That decision is deliberately a CONDITIONAL GO, not an unconditional one:
`04_quality_control/04_estimate_transcript_spillover.py` (transcript spillover) and 4.05 (resegmentation robustness) were
explicitly deferred until Cell Type Annotation cell-type annotation exists, and the blueprint's own Phase
4 completion gate requires both before any headline tumour-engagement result
(Clone Spatial Descriptors+/13) can be considered protected against segmentation-spillover
artefacts. Phases 5 and 6 do not depend on 4.04/4.05 and are correctly
unblocked by this report; Clone Spatial Descriptors+/13 remain explicitly blocked until a
follow-up QC release amends this decision once 4.04/4.05 complete.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError

REQUIRED_INPUTS = {
    "cell_qc_metrics": "data/derived/cell_qc_metrics.parquet",
    "transcript_qc_metrics": "data/derived/transcript_qc_metrics.parquet",
    "fov_artifact_candidates": "reports/qc/spatial_artifact_masks/fov_artifact_candidates.tsv",
    "segmentation_review": "reports/qc/segmentation_review.parquet",
    "qc_thresholds": "config/qc_thresholds.yaml",
    "exclusion_log": "data/derived/exclusion_log.tsv",
    "replicate_concordance": "data/derived/replicate_concordance.tsv",
}

# Sections flagged by 2+ independent QC checks across SpatialData Import/4 as warranting
# extra scrutiny, carried forward here rather than re-derived.
SECTIONS_WITH_CONVERGENT_QC_FLAGS = {
    "P13": (
        "Flagged by 3 independent checks: unusually dense/closely-packed nuclei "
        "(`03_spatialdata_import/02_validate_coordinate_systems.py` coordinate-validation profiling), elevated invalid-nucleus-"
        "polygon rate (`04_quality_control/03_assess_segmentation_quality.py`, 6.8%, among the highest of all 18 sections), "
        "and the weakest replicate-concordance pair on 3 metrics though below "
        "the formal MAD-outlier threshold (`04_quality_control/08_assess_replicate_concordance.R`)."
    ),
}


def _check_inputs_exist(project_root: Path) -> dict[str, Path]:
    paths = {}
    for name, rel_path in REQUIRED_INPUTS.items():
        path = project_root / rel_path
        if not path.is_file():
            raise PipelineError(
                f"'{path}' not found. Run the corresponding earlier Quality Control script first."
            )
        paths[name] = path
    return paths


def build_section_summary(paths: dict[str, Path]) -> pd.DataFrame:
    cell_qc = pd.read_parquet(paths["cell_qc_metrics"])
    per_section_cell = cell_qc.groupby("section_id", observed=True).agg(
        n_cells_raw=("transcript_counts", "size"),
        median_transcript_counts=("transcript_counts", "median"),
        median_genes_detected=("n_genes_detected", "median"),
    )

    transcript_qc = pd.read_parquet(paths["transcript_qc_metrics"]).set_index("section_id")

    fov_df = pd.read_csv(paths["fov_artifact_candidates"], sep="\t")
    per_section_fov = fov_df.groupby("section_id", observed=True).agg(
        n_fovs=("fov_name", "size"),
        n_fovs_flagged=("flagged_artifact_candidate", "sum"),
    )

    seg_df = pd.read_parquet(paths["segmentation_review"]).set_index("section_id")

    excl_log = pd.read_csv(paths["exclusion_log"], sep="\t")
    per_section_excl = excl_log.groupby("section_id", observed=True).agg(
        n_cells_after_qc=("qc_pass", "sum"),
        n_cells_excluded=("excluded", "sum"),
    )
    per_section_excl["fraction_excluded"] = (
        per_section_excl["n_cells_excluded"]
        / (per_section_excl["n_cells_after_qc"] + per_section_excl["n_cells_excluded"])
    ).round(4)

    summary = (
        per_section_cell.join(transcript_qc[["fraction_qv_below_20", "fraction_overlaps_nucleus"]])
        .join(per_section_fov)
        .join(
            seg_df[
                [
                    "fraction_invalid_cell_polygon",
                    "fraction_invalid_nucleus_polygon",
                    "fraction_nucleus_not_contained",
                ]
            ]
        )
        .join(per_section_excl)
    )
    summary["convergent_qc_flag"] = [
        SECTIONS_WITH_CONVERGENT_QC_FLAGS.get(pid, "")
        for pid in summary.index.str.split("_").str[0]
    ]
    return summary.reset_index().rename(columns={"index": "section_id"})


def build_go_no_go_decision(replicate_concordance_path: Path) -> dict:
    concordance = pd.read_csv(replicate_concordance_path, sep="\t")
    n_flagged = int(concordance["flagged_discordant"].sum())
    flagged_pairs = concordance.loc[concordance["flagged_discordant"], "patient_id"].tolist()

    return {
        "status": "CONDITIONAL GO",
        "completed_workstreams": [
            "4.00 Cell-level QC metrics",
            "4.01 Transcript-level QC metrics",
            "4.02 Spatial QC artefact detection",
            "4.03 Segmentation quality assessment",
            "4.06 QC threshold definition (documented hierarchical methodology)",
            "4.07 QC filter application with audit trail",
            "4.08 Replicate concordance assessment",
        ],
        "outstanding_workstreams": [
            "4.04 Transcript spillover estimation -- deferred until Cell Type Annotation cell-type "
            "labels exist",
            "4.05 Resegmentation robustness subset -- deferred alongside 4.04",
        ],
        "unblocked_for": [
            "Preprocessing and Normalisation (expression normalisation/preprocessing)",
            "Cell Type Annotation (cell-type annotation)",
        ],
        "blocked_for": [
            "Clone Spatial Descriptors+ / Clone Ecology Confirmatory Models tumour-engagement claims -- the Quality Control completion "
            "gate requires every headline tumour-engagement result to be cross-checked "
            "against 4.04's spillover-risk output and re-run on 4.05's resegmented "
            "subset before it can be reported. This QC release does not satisfy that "
            "gate; a follow-up release must amend this decision once 4.04/4.05 complete."
        ],
        "replicate_concordance_pairs_flagged": n_flagged,
        "replicate_concordance_flagged_patients": flagged_pairs,
        "sections_with_convergent_qc_flags": list(SECTIONS_WITH_CONVERGENT_QC_FLAGS.keys()),
    }


def build_qc_release_report(project_root: Path, output_path: Path) -> dict:
    paths = _check_inputs_exist(project_root)

    section_summary = build_section_summary(paths)
    decision = build_go_no_go_decision(paths["replicate_concordance"])

    with open(paths["qc_thresholds"]) as f:
        thresholds_config = yaml.safe_load(f)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    rows_html = "".join(
        f"<tr><td>{r.section_id}</td><td>{int(r.n_cells_raw)}</td>"
        f"<td>{int(r.n_cells_after_qc)}</td><td>{r.fraction_excluded:.4f}</td>"
        f"<td>{r.median_transcript_counts:.0f}</td>"
        f"<td>{int(r.n_fovs_flagged)}/{int(r.n_fovs)}</td>"
        f"<td>{r.fraction_invalid_nucleus_polygon:.4f}</td>"
        f"<td>{r.convergent_qc_flag}</td></tr>"
        for r in section_summary.itertuples()
    )

    outstanding_html = "".join(f"<li>{item}</li>" for item in decision["outstanding_workstreams"])
    blocked_html = "".join(f"<li>{item}</li>" for item in decision["blocked_for"])
    unblocked_html = "".join(f"<li>{item}</li>" for item in decision["unblocked_for"])

    html = f"""<html><head><title>Quality Control release report</title></head><body>
<h1>Quality Control release report</h1>
<h2>Decision: {decision['status']}</h2>
<p>Active QC threshold profile: <b>{thresholds_config['active_profile']}</b>
(see <code>config/qc_thresholds.yaml</code> for the full rationale).</p>
<p><b>Unblocked for:</b></p><ul>{unblocked_html}</ul>
<p><b>Blocked for:</b></p><ul>{blocked_html}</ul>
<p><b>Outstanding workstreams:</b></p><ul>{outstanding_html}</ul>
<p>Replicate concordance: {decision['replicate_concordance_pairs_flagged']} / 7 pairs
formally flagged as discordant (MAD threshold); see
<code>reports/qc/replicate_concordance.pdf</code> for the full per-pair
breakdown, including the P13 soft-outlier finding that did not cross the
formal threshold but is carried forward as a caveat below.</p>
<h2>Section-by-section summary</h2>
<table border='1' cellpadding='6' cellspacing='0'>
<tr><th>Section</th><th>N cells (raw)</th><th>N cells (post-QC)</th>
<th>Frac. excluded</th><th>Median transcript counts</th>
<th>FOVs flagged</th><th>Frac. invalid nucleus polygon</th>
<th>Convergent QC flag</th></tr>
{rows_html}</table>
<p>Full detail: <code>data/derived/cell_qc_metrics.parquet</code>,
<code>data/derived/transcript_qc_metrics.parquet</code>,
<code>reports/qc/spatial_artifact_masks/fov_artifact_candidates.tsv</code>,
<code>reports/qc/segmentation_review.html</code>,
<code>data/derived/exclusion_log.tsv</code>.</p>
</body></html>"""
    output_path.write_text(html)

    section_summary.to_parquet(output_path.with_name("qc_release_section_summary.parquet"))

    return {
        "status": decision["status"],
        "sections_processed": len(section_summary),
        "n_cells_raw_total": int(section_summary["n_cells_raw"].sum()),
        "n_cells_post_qc_total": int(section_summary["n_cells_after_qc"].sum()),
        "outstanding_workstreams": decision["outstanding_workstreams"],
        "sections_with_convergent_qc_flags": decision["sections_with_convergent_qc_flags"],
    }
