"""TCR probe-to-patient mapping and leakage audit (`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`).

`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py` found that every CDR3 probe is physically present on 3-4
patients' panels (a manufacturing-batch artefact, not a design choice ), so "patient-specific"
cannot be read off panel presence. This module determines each probe's
most likely intended patient from detection evidence within T cells
(the biologically appropriate population for a CDR3 signal -- checked
that restricting to T cells rather than all cell
types actually matters: background detection across the whole dataset is
far noisier than within T cells alone) and flags probes whose specificity
is statistically weak or absent, rather than silently assuming every
probe cleanly identifies one patient.

**Finding, checked before choosing a method:** signal is not always a
dramatic single-patient spike. A first look at several probes' per-patient
T-cell detection rates showed candidate patients within a batch group
often have broadly overlapping, low detection rates (e.g. one probe: P01
5.6%, P10 1.5%, P12 6.7%, P28 2.6% -- a gradient, not a clean peak), while
others show near-uniformly low noise-level detection across all
candidates. A simple "highest raw count wins" heuristic would be
unreliable at this scale, especially for low-count probes where sampling
noise dominates. **Method:** for each probe, a two-proportion Fisher's
exact test compares the top-detection-rate candidate patient against the
pooled detection rate of the other candidate patients; the resulting
p-values are FDR-corrected (Benjamini-Hochberg, `statsmodels`) across all
216 probes, consistent with this project's standing multiplicity
discipline (governance/analysis_registry.tsv). A probe is
`intended_patient_identified` only if this corrected test is significant
and in the expected direction (higher in the candidate patient); otherwise
it is `no_significant_specificity` -- a reported outcome, not forced to
resolve to a guess.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from xenium_tcr_ecology.infra.exceptions import PipelineError

FDR_ALPHA = 0.05


def compute_probe_patient_detection(
    probe_counts: np.ndarray, patient_ids: pd.Series, candidate_patients: list[str]
) -> pd.DataFrame:
    """Pure, testable per-(probe already fixed)-patient detection summary,
    restricted to `candidate_patients` (the patients whose panel includes
    this probe). `probe_counts` and `patient_ids` must be the same length
    (one row per T cell)."""
    df = pd.DataFrame({"patient_id": patient_ids, "count": probe_counts})
    df = df[df["patient_id"].isin(candidate_patients)]
    summary = df.groupby("patient_id", observed=True)["count"].agg(
        n_tcells="size", n_detected=lambda s: int((s > 0).sum())
    )
    summary["detection_rate"] = summary["n_detected"] / summary["n_tcells"]
    return summary.reindex(candidate_patients)


def evaluate_patient_specificity(detection_summary: pd.DataFrame) -> dict:
    """Fisher's exact test: the candidate patient with the highest
    detection rate vs. the pooled detection of all other candidates.
    Returns the top patient, its p-value, and the raw contingency counts."""
    valid = detection_summary.dropna(subset=["n_tcells"])
    valid = valid[valid["n_tcells"] > 0]
    if len(valid) < 2:
        return {"top_patient": None, "pvalue": None, "direction_consistent": False}

    top_patient = valid["detection_rate"].idxmax()
    top_row = valid.loc[top_patient]
    others = valid.drop(index=top_patient)
    other_detected = int(others["n_detected"].sum())
    other_n = int(others["n_tcells"].sum())
    if other_n == 0:
        return {"top_patient": top_patient, "pvalue": None, "direction_consistent": False}

    table = [
        [int(top_row["n_detected"]), int(top_row["n_tcells"] - top_row["n_detected"])],
        [other_detected, other_n - other_detected],
    ]
    _, pvalue = fisher_exact(table, alternative="greater")
    other_rate = other_detected / other_n
    direction_consistent = top_row["detection_rate"] > other_rate
    return {
        "top_patient": top_patient,
        "pvalue": float(pvalue),
        "direction_consistent": direction_consistent,
    }


def build_patient_probe_audit(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    registry_path = project_root / "metadata" / "tcr_probe_registry.tsv"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_path = project_root / "reports" / "tcr" / "patient_probe_audit.tsv"

    for p in (matrix_path, registry_path, final_annotations_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    registry = pd.read_csv(registry_path, sep="\t")
    final_annotations = pd.read_parquet(final_annotations_path)

    tcell_ids = final_annotations.index[final_annotations["final_lineage"] == "T_cell"]
    tcell_ids = adata.obs_names.intersection(tcell_ids)
    if len(tcell_ids) == 0:
        raise PipelineError("No T cells found in final_cell_annotations.parquet.")

    probe_names = registry["probe_name"].tolist()
    probe_names_present = [g for g in probe_names if g in adata.var_names]
    tcell_adata = adata[tcell_ids, probe_names_present]
    X = tcell_adata.layers["counts"]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    patient_ids = tcell_adata.obs["patient_id"]

    rows = []
    for i, probe_name in enumerate(probe_names_present):
        candidate_patients = (
            registry.loc[registry["probe_name"] == probe_name, "patients_with_probe"]
            .iloc[0]
            .split(";")
        )
        detection_summary = compute_probe_patient_detection(
            X[:, i], patient_ids, candidate_patients
        )
        test_result = evaluate_patient_specificity(detection_summary)
        rows.append(
            {
                "probe_name": probe_name,
                "n_candidate_patients": len(candidate_patients),
                "top_patient": test_result["top_patient"],
                "top_patient_detection_rate": (
                    round(
                        float(detection_summary.loc[test_result["top_patient"], "detection_rate"]),
                        6,
                    )
                    if test_result["top_patient"] is not None
                    else None
                ),
                "pvalue_raw": test_result["pvalue"],
                "direction_consistent": test_result["direction_consistent"],
            }
        )

    audit = pd.DataFrame(rows)
    testable = audit["pvalue_raw"].notna()
    audit["pvalue_fdr"] = np.nan
    if testable.sum() > 0:
        _, fdr_pvals, _, _ = multipletests(
            audit.loc[testable, "pvalue_raw"], alpha=FDR_ALPHA, method="fdr_bh"
        )
        audit.loc[testable, "pvalue_fdr"] = fdr_pvals

    audit["intended_patient_identified"] = (
        testable & (audit["pvalue_fdr"] < FDR_ALPHA) & audit["direction_consistent"]
    )
    audit["intended_patient"] = np.where(
        audit["intended_patient_identified"], audit["top_patient"], None
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    audit.to_csv(output_path, sep="\t", index=False)

    return {
        "n_probes_audited": len(audit),
        "n_probes_with_identified_patient": int(audit["intended_patient_identified"].sum()),
        "n_probes_no_significant_specificity": int((~audit["intended_patient_identified"]).sum()),
        "fdr_alpha": FDR_ALPHA,
        "n_tcells_used": len(tcell_ids),
        "output_path": str(output_path),
    }
