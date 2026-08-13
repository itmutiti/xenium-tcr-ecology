"""Validates this project's Xenium CDR3-probe-based clone-calling
(`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`-`01_map_
tcr_probes_to_patients.py`) against independent, paired scTCR-seq VDJ
data (`08_tcr_clonal_analysis/09_validate_probe_clones_against_
paired_vdj_ground_truth.py`).

**Supersedes a previously-documented data-availability limitation.**
An earlier check of this project's only companion dataset (GSE287301)
found it to contain a gene-expression matrix only, with no VDJ/TCR-contig
files in that GEO deposit -- correct at the time, since it
covered only the small subset of GSE287301's supplementary files
acquired at the time (`filtered_feature_bc_matrix`, `aggregation.csv`,
`patient_matrix.txt`, the Loupe file) -- GSE287301's full `RAW.tar`
(822MB, not previously pulled) contains per-sample Cell Ranger VDJ
output (16 pooled reactions, GSM8743474-GSM8743489:
`filtered_contig_annotations.csv`, `clonotypes.csv`) for the same
28-patient cohort. This module builds the per-patient VDJ clonotype
table this earlier check found absent.

**Demultiplexing method.** Each of the 16 pooled reactions multiplexes
~4 patients via cell hashing (4 antibody-capture features literally
named "1"-"4" in the aggregated GEX matrix's `features.tsv`,
corresponding to `hash1`-`hash4` in `GSE287301_patient_matrix.txt`).
Per-cell hash assignment: normalized-count argmax with a minimum-margin
confidence filter (ratio of top to second-highest normalized hash
count) -- a simple, standard-practice demultiplexing rule, not a full
HTODemux-equivalent mixture model. Patient-number-to-`P0N`-ID mapping
(the primary Xenium cohort's own `patient_id` convention,
`metadata/sample_manifest.tsv`) is by direct numeric correspondence: all
11 primary-cohort patient numbers extracted from `metadata/tcr_probe_
registry.tsv`'s `patients_with_probe` column (1, 9, 10, 12, 13, 15, 17,
19, 20, 23, 28) match the 11 distinct patient numbers appearing across
`GSE287301_patient_matrix.txt` restricted to hash/pool combinations
present in this project's probe registry. The join is verified again
downstream: every matched CDR3 sequence is checked for chain/amino-acid
identity, not just a numeric patient-ID match.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse

from xenium_tcr_ecology.infra.exceptions import PipelineError

HASH_FEATURE_NAMES = ["1", "2", "3", "4"]
MIN_HASH_CONFIDENCE_RATIO = 2.0
MIN_HASH_UMI_COUNT = 10
POOL_NAMES = [
    "chip1pool1",
    "chip1pool2",
    "chip1pool3",
    "chip1pool4",
    "chip1pool5",
    "chip1pool6",
    "chip1pool7",
    "chip1pool8",
    "chip2pool1",
    "chip2pool2",
    "chip2pool3",
    "chip2pool4",
    "chip2pool5",
    "chip2pool6",
    "chip2pool7",
    "chip2pool16",
]


def load_aggregation_order(aggregation_path: Path) -> dict[int, str]:
    """Barcode-suffix (`-N`) to physical sample-name mapping, from
    `cellranger aggr`'s convention (row order in the aggregation CSV =
    suffix order, 1-indexed)."""
    agg = pd.read_csv(aggregation_path)
    return {i + 1: sample_id for i, sample_id in enumerate(agg["sample_id"])}


def load_pool_to_patients(patient_matrix_path: Path) -> dict[str, dict[str, int]]:
    """`pool_key -> {hash1..hash4: patient_number}` mapping from
    `GSE287301_patient_matrix.txt`."""
    matrix = pd.read_csv(patient_matrix_path, sep="\t", index_col=0)
    return {
        pool: {hash_name: int(matrix.loc[hash_name, pool]) for hash_name in matrix.index}
        for pool in matrix.columns
    }


def pool_key_from_sample_name(sample_name: str) -> str:
    """Physical sample name (e.g. `chip1pool3`) to its patient-
    composition pool key (`pool3`) -- pools 1-7 are technical replicates
    across both chips (same 4 patients, re-run), so the `chipN` prefix
    is not part of the patient-composition identity."""
    return (
        sample_name.split("chip1")[-1].split("chip2")[-1]
        if "chip1" in sample_name or "chip2" in sample_name
        else sample_name
    )


def assign_cell_hashes(hash_counts: np.ndarray) -> np.ndarray:
    """Per-cell hash assignment (0-indexed hash number, or -1 for
    unassigned/ambiguous cells) from an `n_cells x 4` raw hash-UMI count
    matrix. A cell is assigned to its top hash only if the top count
    clears `MIN_HASH_UMI_COUNT` and exceeds the second-highest count by
    at least `MIN_HASH_CONFIDENCE_RATIO`-fold -- a simple,
    standard-practice demultiplexing rule (argmax with a confidence
    margin), not a full mixture-model demultiplexing algorithm."""
    sorted_counts = np.sort(hash_counts, axis=1)
    top = sorted_counts[:, -1]
    second = sorted_counts[:, -2]
    top_idx = np.argmax(hash_counts, axis=1)

    confident = (top >= MIN_HASH_UMI_COUNT) & (
        top >= MIN_HASH_CONFIDENCE_RATIO * np.maximum(second, 1)
    )
    assignment = np.where(confident, top_idx, -1)
    return assignment


def load_hash_assignments(matrix_dir: Path, aggregation_path: Path) -> pd.DataFrame:
    """Per-cell hash assignment for every cell in the aggregated
    GSE287301 GEX matrix, joined to its physical pool sample name via
    barcode suffix. Returns columns `barcode_prefix` (16bp, no `-N`
    suffix), `pool_sample_name`, `hash_index` (0-3, or -1 if
    unassigned)."""
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    features_path = matrix_dir / "features.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    for p in (barcodes_path, features_path, matrix_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found -- GSE287301 aggregated matrix is incomplete.")

    with gzip.open(barcodes_path, "rt") as f:
        barcodes = [line.strip() for line in f]
    with gzip.open(features_path, "rt") as f:
        features = [line.strip().split("\t") for line in f]

    hash_row_indices = [
        i
        for i, feat in enumerate(features)
        if len(feat) >= 3 and feat[2] == "Antibody Capture" and feat[0] in HASH_FEATURE_NAMES
    ]
    if len(hash_row_indices) != 4:
        raise PipelineError(
            f"Expected 4 'Antibody Capture' hash features, found {len(hash_row_indices)} -- matrix format may have changed."
        )

    matrix = scipy.io.mmread(str(matrix_path)).tocsr()

    hash_counts = np.asarray(matrix[hash_row_indices, :].todense()).T  # cells x 4 hashes

    suffix_to_sample = load_aggregation_order(aggregation_path)
    suffixes = np.array([int(b.split("-")[1]) for b in barcodes])
    barcode_prefixes = np.array([b.split("-")[0] for b in barcodes])
    pool_sample_names = np.array([suffix_to_sample[s] for s in suffixes])

    hash_index = assign_cell_hashes(hash_counts)

    return pd.DataFrame(
        {
            "barcode_prefix": barcode_prefixes,
            "pool_sample_name": pool_sample_names,
            "hash_index": hash_index,
        }
    )


def build_vdj_ground_truth_clonotypes(
    vdj_dir: Path, hash_assignments: pd.DataFrame, pool_to_patients: dict[str, dict[str, int]]
) -> pd.DataFrame:
    """Per-patient VDJ clonotype table: for every pool's
    `filtered_contig_annotations.csv`, joins each productive contig's
    cell barcode to its hash-derived patient assignment, then aggregates
    to per-(patient, chain, cdr3_amino_acid_sequence) cell counts -- the
    independent analogue of this project's Xenium probe-detection
    data."""
    hash_names = ["hash1", "hash2", "hash3", "hash4"]
    rows = []

    for pool_sample_name in POOL_NAMES:
        contig_path = vdj_dir / pool_sample_name / "filtered_contig_annotations.csv"
        if not contig_path.is_file():
            raise PipelineError(
                f"'{contig_path}' not found -- VDJ acquisition for '{pool_sample_name}' is incomplete."
            )

        pool_key = pool_key_from_sample_name(pool_sample_name)
        if pool_key not in pool_to_patients:
            raise PipelineError(
                f"Pool key '{pool_key}' (from '{pool_sample_name}') not found in the patient composition matrix."
            )

        cell_to_patient = hash_assignments[
            (hash_assignments["pool_sample_name"] == pool_sample_name)
            & (hash_assignments["hash_index"] >= 0)
        ].copy()
        cell_to_patient["patient_number"] = cell_to_patient["hash_index"].map(
            lambda i: pool_to_patients[pool_key][hash_names[i]]
        )
        barcode_to_patient = dict(
            zip(cell_to_patient["barcode_prefix"], cell_to_patient["patient_number"])
        )

        contigs = pd.read_csv(contig_path)
        contigs = contigs[
            (contigs["productive"] == True) & (contigs["is_cell"] == True)  # noqa: E712
        ]
        contigs["barcode_prefix"] = contigs["barcode"].str.split("-").str[0]
        contigs["patient_number"] = contigs["barcode_prefix"].map(barcode_to_patient)
        contigs = contigs.dropna(subset=["patient_number"])
        contigs["patient_number"] = contigs["patient_number"].astype(int)

        rows.append(contigs[["patient_number", "chain", "cdr3", "barcode"]])

    all_contigs = pd.concat(rows, ignore_index=True)
    all_contigs = all_contigs.rename(columns={"cdr3": "cdr3_amino_acid_sequence"})

    clonotype_counts = (
        all_contigs.groupby(["patient_number", "chain", "cdr3_amino_acid_sequence"])["barcode"]
        .nunique()
        .reset_index(name="n_cells")
    )
    clonotype_counts["patient_id"] = clonotype_counts["patient_number"].map(lambda n: f"P{n:02d}")
    clonotype_counts["rank_within_patient_chain"] = clonotype_counts.groupby(
        ["patient_id", "chain"]
    )["n_cells"].rank(ascending=False, method="min")
    return clonotype_counts


def compare_probe_detections_to_vdj_ground_truth(project_root: Path) -> dict:
    probe_registry_path = project_root / "metadata" / "tcr_probe_registry.tsv"
    patient_probe_audit_path = project_root / "reports" / "tcr" / "patient_probe_audit.tsv"
    matrix_dir = project_root / "data" / "external" / "GSE287301" / "filtered_feature_bc_matrix"
    aggregation_path = (
        project_root / "data" / "external" / "GSE287301" / "GSE287301_aggregation.csv.gz"
    )
    patient_matrix_path = (
        project_root / "data" / "external" / "GSE287301" / "GSE287301_patient_matrix.txt.gz"
    )
    vdj_dir = project_root / "data" / "external" / "GSE287301" / "vdj"
    output_path = project_root / "data" / "derived" / "probe_vdj_ground_truth_comparison.parquet"

    for p in (
        probe_registry_path,
        patient_probe_audit_path,
        matrix_dir,
        aggregation_path,
        patient_matrix_path,
        vdj_dir,
    ):
        if not p.exists():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    probe_registry = pd.read_csv(probe_registry_path, sep="\t")
    audit = pd.read_csv(patient_probe_audit_path, sep="\t")
    merged = probe_registry.merge(audit, on="probe_name")

    hash_assignments = load_hash_assignments(matrix_dir, aggregation_path)
    pool_to_patients = load_pool_to_patients(patient_matrix_path)
    ground_truth = build_vdj_ground_truth_clonotypes(vdj_dir, hash_assignments, pool_to_patients)

    identified = merged[merged["intended_patient_identified"] == True].copy()  # noqa: E712
    gt_lookup = ground_truth.set_index(["patient_id", "chain", "cdr3_amino_acid_sequence"])

    results = []
    for _, row in identified.iterrows():
        key = (row["intended_patient"], row["tcr_chain"], row["cdr3_amino_acid_sequence"])
        found = key in gt_lookup.index
        n_cells = int(gt_lookup.loc[key, "n_cells"]) if found else 0
        rank = float(gt_lookup.loc[key, "rank_within_patient_chain"]) if found else np.nan
        results.append(
            {
                "probe_name": row["probe_name"],
                "intended_patient": row["intended_patient"],
                "tcr_chain": row["tcr_chain"],
                "cdr3_amino_acid_sequence": row["cdr3_amino_acid_sequence"],
                "xenium_detection_rate": row["top_patient_detection_rate"],
                "found_in_real_vdj_ground_truth": found,
                "vdj_ground_truth_n_cells": n_cells,
                "vdj_ground_truth_rank_within_patient_chain": rank,
            }
        )

    comparison = pd.DataFrame(results)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_parquet(output_path)

    ground_truth_path = (
        project_root / "data" / "derived" / "gse287301_vdj_ground_truth_clonotypes.parquet"
    )
    ground_truth.to_parquet(ground_truth_path)

    n_found = int(comparison["found_in_real_vdj_ground_truth"].sum())
    found_subset = comparison[comparison["found_in_real_vdj_ground_truth"]]
    rank_corr = None
    if len(found_subset) >= 3:
        rank_corr = float(
            found_subset["xenium_detection_rate"].corr(
                found_subset["vdj_ground_truth_rank_within_patient_chain"], method="spearman"
            )
        )

    report_path = project_root / "reports" / "tcr" / "probe_vdj_ground_truth_validation.pdf"
    render_ground_truth_validation_report(comparison, report_path)

    return {
        "n_probes_with_identified_patient": len(identified),
        "n_found_in_real_vdj_ground_truth": n_found,
        "fraction_found": round(n_found / len(identified), 4) if len(identified) else None,
        "n_real_vdj_ground_truth_cells_total": int(ground_truth["n_cells"].sum()),
        "n_real_vdj_ground_truth_patients": int(ground_truth["patient_id"].nunique()),
        "xenium_detection_vs_vdj_rank_spearman": rank_corr,
        "output_path": str(output_path),
        "ground_truth_path": str(ground_truth_path),
        "report_path": str(report_path),
    }


def render_ground_truth_validation_report(comparison: pd.DataFrame, output_path: Path) -> None:
    """Two-panel manuscript figure for the probe-vs-paired-scTCR-seq
    comparison. Panel A: per-patient count of probes confirmed vs. not
    confirmed against paired scTCR-seq VDJ data from the same patients
    (a cross-modality concordance check, not a blinded independent
    validation -- see Methods for the circularity this comparison cannot
    rule out). Panel B: Xenium spatial detection rate vs. VDJ-derived
    clonal rank for every confirmed probe, with the Spearman correlation
    and its p-value, computed identically to Results."""
    import matplotlib
    from scipy import stats

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from xenium_tcr_ecology.viz.style import (
        COLORS,
        FS_ANNOTATION,
        apply_publication_style,
        panel_label,
        panel_title,
    )

    apply_publication_style()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    # Side-by-side, not stacked: for 2 panels a single moderate-width
    # row keeps the aspect ratio closer to square than either a wide
    # 2-in-a-row (too wide once shrunk to fit a word-processor page) or
    # a tall 2x1 stack (too long/elongated) would.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.6, 6.6))

    per_patient = (
        comparison.groupby(["intended_patient", "found_in_real_vdj_ground_truth"])
        .size()
        .unstack(fill_value=0)
    )
    for col in (True, False):
        if col not in per_patient.columns:
            per_patient[col] = 0
    per_patient["total"] = per_patient[True] + per_patient[False]
    per_patient = per_patient.sort_values("total", ascending=False)

    x = range(len(per_patient))
    ax1.bar(x, per_patient[True], label="Confirmed", color=COLORS["confirmed"], width=0.7)
    ax1.bar(
        x,
        per_patient[False],
        bottom=per_patient[True],
        label="Not confirmed",
        color=COLORS["not_confirmed"],
        width=0.7,
    )
    ax1.set_xticks(list(x))
    ax1.set_xticklabels(per_patient.index, rotation=45, ha="right")
    ax1.set_ylabel("Patient-identified CDR3 probes")
    ax1.tick_params(axis="both", which="major", pad=5)
    panel_title(ax1, "Confirmation by patient")
    ax1.legend(loc="upper right", handlelength=1.6, handletextpad=0.6)
    panel_label(ax1, "A")

    found = comparison[comparison["found_in_real_vdj_ground_truth"]]
    rho, pvalue = stats.spearmanr(
        found["xenium_detection_rate"], found["vdj_ground_truth_rank_within_patient_chain"]
    )
    ax2.scatter(
        found["vdj_ground_truth_rank_within_patient_chain"],
        found["xenium_detection_rate"],
        color=COLORS["confirmed"],
        s=42,
        alpha=0.75,
        linewidths=0,
    )
    ax2.set_xlabel("VDJ-derived rank within patient/chain\n(1 = most abundant)")
    ax2.set_ylabel("Xenium detection rate")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.tick_params(axis="both", which="major", pad=5)
    y_min, y_max = ax2.get_ylim()
    ax2.set_ylim(
        y_min, y_max * 2.4
    )  # headroom so the top-right point and annotation don't crowd the frame
    panel_title(ax2, f"Detection rate vs. VDJ rank (n = {len(found)})")
    ax2.annotate(
        f"Spearman ρ = {rho:.2f}, p = {pvalue:.2f}",
        xy=(0.97, 0.94),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=FS_ANNOTATION,
    )
    panel_label(ax2, "B")

    fig.tight_layout(w_pad=3.6)
    fig.savefig(output_path)
    fig.savefig(output_path.with_suffix(".png"), dpi=600)
    plt.close(fig)
