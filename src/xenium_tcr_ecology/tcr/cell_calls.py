"""Cell-level TCR detection calling (`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`).

For each of `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s 105 empirically identified probes (probes with a
statistically significant intended patient), calls a per-cell clone
detection within that probe's intended patient's OWN T cells (Phase
6.06's `final_lineage == "T_cell"`) -- restricting evaluation to intended
specimens, per `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s own title, applied here for the first time as
an actual downstream constraint, not just an audit finding.

**Threshold:** `MIN_COUNT_THRESHOLD = 1` (any nonzero raw transcript
count), consistent with `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s methodology, which already
established statistically significant patient-specificity using this same
binary detected/not-detected definition. Checked before
finalising: 77.2% of all nonzero detections among the
105 identified probes' intended-patient T cells are exactly 1 transcript
(7,847 of 10,171) -- a classic sparse targeted-probe distribution, not
evidence that single-transcript detections are noise (the Fisher's-exact
patient-specificity test in `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py` was already computed at this same
threshold and found signal). The raw `detection_count` is retained
alongside the binary call so a stricter downstream threshold remains
possible without re-deriving from raw data.

**Controls:** the 111/216 probes without an identified intended patient
(`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`) are not called here at all -- there is no statistically
defensible patient assignment to restrict evaluation to, and calling them
anyway would silently reintroduce the exact panel-presence-is-not-intent
problem `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py` exists to correct. Off-patient background detection
(the same candidate-but-not-intended patients already characterised in
`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`) is not recomputed here; it is `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s
`top_patient_detection_rate`-vs-others comparison, already on record.

**Multi-probe ambiguity:** flags (does not resolve -- `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py`'s
explicit job) T cells with a detection for more than one distinct
in-scope probe simultaneously.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

MIN_COUNT_THRESHOLD = 1


def call_cell_detections(
    counts: np.ndarray,
    probe_names: list[str],
    probe_intended_patient: dict[str, str],
    cell_patient_ids: pd.Series,
    probe_chain: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Pure, testable core calling logic. `counts` is a (n_cells x
    n_probes) array aligned with `probe_names` (columns) and
    `cell_patient_ids` (rows, one per cell). For each cell, only the
    columns whose probe is intended for that cell's own patient are
    considered "in scope" -- a probe never contributes a detection call
    for a cell from a different patient, even if physically present on
    that patient's panel (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`, `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s batch-sharing finding).
    """
    n_cells, n_probes = counts.shape
    detected = counts >= MIN_COUNT_THRESHOLD

    probe_patient_arr = np.array([probe_intended_patient[p] for p in probe_names])
    cell_patient_arr = cell_patient_ids.to_numpy()
    in_scope = cell_patient_arr[:, None] == probe_patient_arr[None, :]

    in_scope_detected = detected & in_scope
    n_probes_evaluated = in_scope.sum(axis=1)
    n_probes_detected = in_scope_detected.sum(axis=1)

    detected_probe_names = []
    likely_single_clone_pair = np.zeros(n_cells, dtype=bool)
    for i in range(n_cells):
        names = [probe_names[j] for j in np.where(in_scope_detected[i])[0]]
        detected_probe_names.append(";".join(names))
        if probe_chain is not None and len(names) == 2:
            chains = {probe_chain[n] for n in names}
            likely_single_clone_pair[i] = chains == {"TRA", "TRB"}

    result = pd.DataFrame(
        {
            "n_probes_evaluated": n_probes_evaluated,
            "n_probes_detected": n_probes_detected,
            "detected_probes": detected_probe_names,
            "is_multi_probe_ambiguous": n_probes_detected > 1,
            "any_detection": n_probes_detected > 0,
        }
    )
    if probe_chain is not None:
        # A cell with exactly one detected TRA probe and one detected TRB
        # probe is consistent with a single clone's normally-paired
        # alpha-beta TCR, not genuine multi-clonal ambiguity -- checked
        # on the data before adding this field: 63.1% of
        # "is_multi_probe_ambiguous" cells fall into exactly this pattern.
        # This is a useful
        # distinction to carry forward for `08_tcr_clonal_analysis/06_resolve_multiclonal_and_ambiguous_cells.py` to actually resolve,
        # not a claim that pairing is confirmed -- no ground-truth
        # TRA/TRB linkage exists in this project's data (`08_tcr_clonal_analysis/02_document_clone_ascertainment.py`).
        result["is_multi_probe_ambiguous_excluding_likely_pairs"] = (
            result["is_multi_probe_ambiguous"] & ~likely_single_clone_pair
        )
        result["likely_single_clone_tra_trb_pair"] = likely_single_clone_pair
    return result


def build_tcr_cell_calls(project_root: Path) -> dict:
    matrix_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    ascertainment_path = project_root / "metadata" / "clone_ascertainment.tsv"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    output_path = project_root / "data" / "derived" / "tcr_cell_calls.parquet"

    for p in (matrix_path, ascertainment_path, final_annotations_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    import anndata as ad

    adata = ad.read_h5ad(matrix_path)
    ascertainment = pd.read_csv(ascertainment_path, sep="\t")
    final_annotations = pd.read_parquet(final_annotations_path)

    identified = ascertainment[ascertainment["intended_patient_identified"]]
    if len(identified) == 0:
        raise PipelineError(
            "No probes with an identified intended patient (`08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`) -- cannot call detections."
        )

    tcell_ids = final_annotations.index[final_annotations["final_lineage"] == "T_cell"]
    tcell_ids = adata.obs_names.intersection(tcell_ids)
    if len(tcell_ids) == 0:
        raise PipelineError("No T cells found in final_cell_annotations.parquet.")

    probe_names = [g for g in identified["probe_name"] if g in adata.var_names]
    probe_intended_patient = identified.set_index("probe_name")["intended_patient"].to_dict()
    probe_chain = identified.set_index("probe_name")["tcr_chain"].to_dict()

    tcell_adata = adata[tcell_ids, probe_names]
    X = tcell_adata.layers["counts"]
    counts = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    cell_patient_ids = tcell_adata.obs["patient_id"]

    calls = call_cell_detections(
        counts, probe_names, probe_intended_patient, cell_patient_ids, probe_chain
    )
    calls.index = tcell_adata.obs_names
    calls.index.name = "cell_id"
    calls["patient_id"] = cell_patient_ids.to_numpy()
    calls["section_id"] = tcell_adata.obs["section_id"].to_numpy()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    calls.to_parquet(output_path)

    return {
        "n_tcells": len(calls),
        "n_probes_used": len(probe_names),
        "n_tcells_with_any_detection": int(calls["any_detection"].sum()),
        "n_tcells_multi_probe_ambiguous": int(calls["is_multi_probe_ambiguous"].sum()),
        "n_tcells_likely_single_clone_tra_trb_pair": int(
            calls["likely_single_clone_tra_trb_pair"].sum()
        ),
        "n_tcells_multi_probe_ambiguous_excluding_likely_pairs": int(
            calls["is_multi_probe_ambiguous_excluding_likely_pairs"].sum()
        ),
        "fraction_tcells_with_detection": round(float(calls["any_detection"].mean()), 4),
        "detection_counts_by_patient": {
            str(k): int(v)
            for k, v in calls.groupby("patient_id", observed=True)["any_detection"].sum().items()
        },
        "output_path": str(output_path),
    }
