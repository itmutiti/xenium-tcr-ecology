"""False-positive TCR call estimation -- Python input preparation (Phase
8.04 helper).

`04_estimate_false_positive_tcr_calls.R` cannot read `.h5ad` directly (no
R HDF5/AnnData reader available),
so this module -- invoked from R via `system2()`, matching the
`_02_compute_normalization_benchmark_metrics.py` (`05_preprocessing_and_normalisation/02_evaluate_normalisation_strategies.R`) and
`_03_prepare_cnv_inputs.py` (`07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R`) precedent -- computes and exports
three independent empirical negative controls for each of Phase 8.01's
105 identified probes:

1. **Off-patient control:** the same other-candidate-patient detection
   rate already computed once in `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py` (`compute_probe_patient_detection`,
   reused directly here, not recomputed with different logic) -- a probe
   with signal should show near-zero detection in the same batch's
   non-intended patients' own T cells.
2. **Non-T-cell control:** detection rate for each probe among its own
   intended patient's non-T cells (all other lineages, `06_cell_type_annotation/06_integrate_annotation_evidence.py`). A
   CDR3 probe with signal should not be expressed outside the
   T-cell lineage; elevated non-T-cell detection indicates background/
   technical noise, not patient misassignment.
3. **Spatial control:** Moran's I spatial autocorrelation (squidpy, the
   same method and rationale already used and validated in Phase 7.02 --
   noise has no reason to be spatially clustered; a genuine clonal
   T-cell population, having expanded and/or migrated together, plausibly
   does) of each probe's binary detection status among its own intended
   patient's T cells, computed per section (not pooled across a multi-
   section patient's sections, which are physically different tissue
   pieces). Only computed for probes with at least
   `MIN_CELLS_FOR_SPATIAL_TEST` detected cells in a given section --
   checked against the data first: 90/105 probes have >=10
   detected cells overall (median 34, max 902), so this floor excludes a
   minority, not the majority, of probes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tcr.patient_mapping import compute_probe_patient_detection

MIN_CELLS_FOR_SPATIAL_TEST = 10
SPATIAL_KNN_NEIGHBORS = 6


def compute_non_tcell_detection_rate(
    probe_counts: np.ndarray, patient_ids: pd.Series, intended_patient: str
) -> dict:
    """Pure, testable non-T-cell background rate for one probe. Expects
    `probe_counts`/`patient_ids` restricted to non-T cells already."""
    mask = (patient_ids == intended_patient).to_numpy()
    n_cells = int(mask.sum())
    if n_cells == 0:
        return {"n_non_tcells": 0, "n_detected": 0, "detection_rate": None}
    detected = int((probe_counts[mask] > 0).sum())
    return {"n_non_tcells": n_cells, "n_detected": detected, "detection_rate": detected / n_cells}


def compute_detection_spatial_autocorrelation(
    x: np.ndarray, y: np.ndarray, detected: np.ndarray, n_neighs: int = SPATIAL_KNN_NEIGHBORS
) -> dict:
    """Pure(ish) Moran's I of a binary detection variable, matching Phase
    7.02's `compute_malignancy_spatial_autocorrelation` pattern exactly
    but on detection status rather than a continuous malignancy score."""
    n = len(x)
    if n < n_neighs + 1 or detected.sum() == 0 or detected.sum() == n:
        return {"n_cells": n, "morans_i": None, "pval_norm": None}

    import anndata as ad
    import squidpy as sq

    sub = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32))
    sub.obsm["spatial"] = np.column_stack([x, y])
    sub.obs["detected"] = detected.astype(float)

    sq.gr.spatial_neighbors_knn(sub, n_neighs=n_neighs)
    result = sq.gr.spatial_autocorr(
        sub, mode="moran", attr="obs", genes=["detected"], copy=True, show_progress_bar=False
    )
    row = result.loc["detected"]
    return {"n_cells": n, "morans_i": float(row["I"]), "pval_norm": float(row["pval_norm"])}


def prepare_false_positive_inputs(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    ascertainment_path = project_root / "metadata" / "clone_ascertainment.tsv"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_path = project_root / "data" / "derived" / "tcr_false_positive_controls.parquet"

    for p in (matrix_path, ascertainment_path, final_annotations_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    ascertainment = pd.read_csv(ascertainment_path, sep="\t")
    final_annotations = pd.read_parquet(final_annotations_path)

    identified = ascertainment[ascertainment["intended_patient_identified"]]
    tcell_ids = final_annotations.index[final_annotations["final_lineage"] == "T_cell"]
    tcell_ids = adata.obs_names.intersection(tcell_ids)
    non_tcell_ids = adata.obs_names.difference(tcell_ids)

    rows = []
    for _, probe_row in identified.iterrows():
        gene = probe_row["probe_name"]
        intended_patient = probe_row["intended_patient"]
        candidate_patients = probe_row["patients_with_probe"].split(";")
        if gene not in adata.var_names:
            continue

        X_all = adata[:, gene].layers["counts"]
        all_counts = (
            X_all.toarray().ravel() if hasattr(X_all, "toarray") else np.asarray(X_all).ravel()
        )

        tcell_mask = adata.obs_names.isin(tcell_ids)
        off_patient = compute_probe_patient_detection(
            all_counts[tcell_mask], adata.obs.loc[tcell_mask, "patient_id"], candidate_patients
        )
        other_rows = off_patient.drop(index=intended_patient, errors="ignore")
        off_patient_detected = int(other_rows["n_detected"].sum())
        off_patient_n = int(other_rows["n_tcells"].sum())

        non_tcell_mask = adata.obs_names.isin(non_tcell_ids)
        non_tcell_result = compute_non_tcell_detection_rate(
            all_counts[non_tcell_mask],
            adata.obs.loc[non_tcell_mask, "patient_id"],
            intended_patient,
        )

        own_mask = tcell_mask & (adata.obs["patient_id"] == intended_patient).to_numpy()
        own_detected = (all_counts[own_mask] > 0).astype(int)
        own_x = adata.obs.loc[own_mask, "x_centroid"].to_numpy()
        own_y = adata.obs.loc[own_mask, "y_centroid"].to_numpy()
        own_section = adata.obs.loc[own_mask, "section_id"]
        best_spatial = {"n_cells": 0, "morans_i": None, "pval_norm": None}
        for section_id, section_mask in own_section.groupby(
            own_section, observed=True
        ).groups.items():
            idx = own_section.index.get_indexer(section_mask)
            if own_detected[idx].sum() < MIN_CELLS_FOR_SPATIAL_TEST:
                continue
            spatial_result = compute_detection_spatial_autocorrelation(
                own_x[idx], own_y[idx], own_detected[idx]
            )
            if spatial_result["morans_i"] is not None:
                best_spatial = spatial_result
                break

        rows.append(
            {
                "probe_name": gene,
                "intended_patient": intended_patient,
                "off_patient_n_tcells": off_patient_n,
                "off_patient_n_detected": off_patient_detected,
                "off_patient_detection_rate": (
                    (off_patient_detected / off_patient_n) if off_patient_n else None
                ),
                "non_tcell_n_cells": non_tcell_result["n_non_tcells"],
                "non_tcell_n_detected": non_tcell_result["n_detected"],
                "non_tcell_detection_rate": non_tcell_result["detection_rate"],
                "own_n_tcells": int(own_mask.sum()),
                "own_n_detected": int(own_detected.sum()),
                "own_detection_rate": float(own_detected.mean()) if own_mask.sum() else None,
                "spatial_morans_i": best_spatial["morans_i"],
                "spatial_pval_norm": best_spatial["pval_norm"],
            }
        )

    result = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path)

    return {
        "n_probes": len(result),
        "n_probes_with_spatial_test": int(result["spatial_morans_i"].notna().sum()),
        "output_path": str(output_path),
    }
