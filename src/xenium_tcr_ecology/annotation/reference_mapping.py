"""External scRNA-seq reference mapping with label transfer (`06_cell_type_annotation/03_map_external_scrna_reference.py`).

Reference: GSE287301 (McCord et al. 2026's own companion scRNA-seq dataset --
T cells from the same 28-patient HNSCC cohort as this study's Xenium data,
366,632 cells x 38,606 genes, Cell Ranger 8.0.1 filtered output). This is
the ideal available reference: matched cohort, matched disease context --
not an unrelated public atlas.

The paper's own cell-type/cluster labels (its reported "14 transcriptionally
distinct T-cell clusters") are embedded in a proprietary 10x Loupe Browser
file (.cloupe) in this GEO deposit, which is not a scriptable, parseable
format from outside Loupe Browser itself -- confirmed by inspecting the
GEO supplementary file listing before starting this module (no separate
tabular cell-type annotation file is provided). Reference labels are
therefore derived independently here, via marker-based scoring on the
reference's own whole transcriptome (using literature T-cell-state markers
far richer than the 623-gene Xenium panel could support alone -- the whole
point of using an external reference). This is an approximation of, not a
reproduction of, the paper's own 14-cluster structure -- stated explicitly,
not implied to be the same thing.

Label transfer uses a nearest-centroid Pearson-correlation classifier
(the same family of method as SingleR), not a joint-embedding/kNN approach:
robust to the substantial platform and normalisation differences between
whole-transcriptome scRNA-seq and a targeted spatial panel, which a shared
PCA space would be more sensitive to. Transfer-confidence degradation from
restricting to Xenium-panel-overlapping genes is quantified:
held-out reference cells (with known labels) are classified using
(a) the full reference gene set and (b) only the panel-overlap genes, and
the accuracy difference is reported -- not asserted, measured.
"""

from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import anndata as ad
import fast_matrix_market as fmm
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.model_selection import train_test_split

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.seeding import get_annotation_seed

MIN_GENES_PER_CELL = 200
MIN_COUNTS_PER_CELL = 500
MAX_MT_FRACTION = 0.20
RNG_SEED = get_annotation_seed()

# Literature T-cell-state markers, using whole-transcriptome access (richer
# than the Xenium panel alone could support) -- an independently-curated
# 7-state approximation of the source paper's own reported finer structure
# (14 clusters, inaccessible here -- see module docstring), not a
# reproduction of it.
T_CELL_STATE_MARKERS: dict[str, list[str]] = {
    "Naive_CM": ["CCR7", "TCF7", "SELL", "LEF1"],
    "Memory_TRM": ["IL7R", "ITGAE", "CXCR6", "ZNF683"],
    "Cytotoxic_effector": ["GZMB", "PRF1", "GNLY", "KLRG1", "FGFBP2"],
    "Exhausted": ["PDCD1", "HAVCR2", "LAG3", "TOX", "CXCL13", "ENTPD1", "TIGIT"],
    "Treg": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "Cycling": ["MKI67", "TOP2A", "TUBA1B"],
    "MAIT_unconventional": ["SLC4A10", "KLRB1", "ZBTB16"],
}
MIN_MARKERS_PER_STATE = 2


def extract_reference_matrix(archive_path: Path, extract_dir: Path) -> Path:
    if not archive_path.is_file():
        raise PipelineError(f"'{archive_path}' not found.")
    extract_dir.mkdir(parents=True, exist_ok=True)
    if not (extract_dir / "matrix.mtx.gz").is_file():
        with tarfile.open(archive_path) as tar:
            tar.extractall(extract_dir)
    return extract_dir


def load_reference_matrix(matrix_dir: Path) -> ad.AnnData:
    features_path = matrix_dir / "features.tsv.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    for p in (features_path, barcodes_path, matrix_path):
        if not p.is_file():
            raise PipelineError(f"'{p}' not found. Run extract_reference_matrix() first.")

    features = pd.read_csv(
        features_path, sep="\t", header=None, names=["gene_id", "gene_symbol", "feature_type"]
    )
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, names=["barcode"])

    # scipy.io.mmread's pure-Python parser is impractically slow on this
    # scale (706M nonzero entries for the GSE287301 matrix); fast_matrix_market
    # is a drop-in, multi-threaded C++ replacement -- see environment/conda/main.yml.
    with gzip.open(matrix_path, "rb") as f:
        X = fmm.mmread(f).T.tocsr()

    gex_mask = (features["feature_type"] == "Gene Expression").to_numpy()
    X = X[:, gex_mask]
    gex_features = features.loc[gex_mask].reset_index(drop=True)

    adata = ad.AnnData(X=X.astype(np.float32))
    adata.obs_names = barcodes["barcode"].to_numpy()
    adata.var_names = gex_features["gene_symbol"].to_numpy()
    adata.var["gene_id"] = gex_features["gene_id"].to_numpy()
    adata.var_names_make_unique()
    return adata


def qc_filter_reference(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    n_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    n_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    mt_fraction = adata.obs["pct_counts_mt"].to_numpy() / 100.0

    keep = (
        (n_genes >= MIN_GENES_PER_CELL)
        & (n_counts >= MIN_COUNTS_PER_CELL)
        & (mt_fraction <= MAX_MT_FRACTION)
    )
    return adata[keep].copy()


def normalize_reference(adata: ad.AnnData) -> ad.AnnData:
    adata = adata.copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def validate_t_cell_state_markers(
    available_genes: set[str], markers: dict[str, list[str]] = T_CELL_STATE_MARKERS
) -> dict[str, list[str]]:
    validated = {}
    for state, genes in markers.items():
        present = [g for g in genes if g in available_genes]
        if len(present) < MIN_MARKERS_PER_STATE:
            raise PipelineError(
                f"T-cell state '{state}' has only {len(present)} marker(s) present."
            )
        validated[state] = present
    return validated


def annotate_t_cell_states(adata: ad.AnnData, rng_seed: int = RNG_SEED) -> pd.Series:
    """Marker-score-based direct annotation (score_genes + argmax per cell),
    not a full clustering-then-manual-labelling workflow: a deliberate scope
    simplification given this reference's sole purpose here is to seed label
    transfer, not to be a definitive T-cell atlas in its own right."""
    validated = validate_t_cell_state_markers(set(adata.var_names))
    scores = pd.DataFrame(index=adata.obs_names)
    for state, genes in validated.items():
        score_name = f"{state}_score"
        sc.tl.score_genes(adata, gene_list=genes, score_name=score_name, random_state=rng_seed)
        scores[state] = adata.obs[score_name]
    return scores.idxmax(axis=1)


def standardize_for_transfer(expr: pd.DataFrame) -> pd.DataFrame:
    """Z-scores each gene using this dataframe's own population mean/std.

    Required before cross-platform correlation-based comparison: Pearson
    correlation row-centers each cell (scale-invariant to a *linear*
    rescaling of a cell's own profile), but does nothing to correct for
    *per-gene* scale differences between two platforms with very different
    depth/normalisation conventions -- confirmed as an actual, not
    hypothetical, problem on this project's own data: comparing raw
    log-normalised values directly between GSE287301's scRNA reference
    (target_sum=1e4) and the Xenium data's lognorm layer (target_sum=197)
    produced a degenerate result (65.3% of all Xenium cells predicted as
    "Cycling", vs a plausible 9.3% in the reference itself), because log1p
    is not linear, so a fixed additive/multiplicative scale mismatch in
    normalised-count space becomes a *shape*-distorting, non-removable
    difference after log-transformation. Standardising each gene within its
    own platform's population removes this before comparison.
    """
    mean = expr.mean(axis=0)
    std = expr.std(axis=0)
    std_safe = std.where(std > 0, 1.0)
    return (expr - mean) / std_safe


def compute_state_centroids(
    adata: ad.AnnData, labels: pd.Series, gene_subset: list[str] | None = None
) -> pd.DataFrame:
    genes = gene_subset if gene_subset is not None else list(adata.var_names)
    genes = [g for g in genes if g in adata.var_names]
    sub = adata[:, genes]
    X = sub.X.toarray() if sparse.issparse(sub.X) else np.asarray(sub.X)
    df = pd.DataFrame(X, index=adata.obs_names, columns=genes)
    df = standardize_for_transfer(df)
    df["label"] = labels.reindex(adata.obs_names).to_numpy()
    return df.groupby("label").mean()


def classify_by_nearest_centroid(query_expr: pd.DataFrame, centroids: pd.DataFrame) -> pd.DataFrame:
    query_expr = standardize_for_transfer(query_expr)
    shared_genes = [g for g in centroids.columns if g in query_expr.columns]
    if len(shared_genes) == 0:
        raise PipelineError("No shared genes between query expression and centroids.")
    q = query_expr[shared_genes].to_numpy()
    c = centroids[shared_genes].to_numpy()

    q_centered = q - q.mean(axis=1, keepdims=True)
    c_centered = c - c.mean(axis=1, keepdims=True)
    q_norm = np.linalg.norm(q_centered, axis=1, keepdims=True)
    c_norm = np.linalg.norm(c_centered, axis=1, keepdims=True)
    q_norm[q_norm == 0] = 1.0
    c_norm[c_norm == 0] = 1.0
    correlations = (q_centered / q_norm) @ (c_centered / c_norm).T

    best_idx = np.argmax(correlations, axis=1)
    predicted = centroids.index.to_numpy()[best_idx]
    confidence = correlations[np.arange(len(q)), best_idx]
    return pd.DataFrame(
        {"predicted_state": predicted, "confidence": confidence}, index=query_expr.index
    )


N_HVG_FOR_BROAD_GENE_SET = 2000
MAX_BENCHMARK_TEST_CELLS = 20000


def benchmark_transfer_confidence_degradation(
    adata: ad.AnnData,
    labels: pd.Series,
    panel_genes: list[str],
    rng_seed: int = RNG_SEED,
    test_size: float = 0.3,
) -> dict:
    """Compares Xenium-panel-restricted transfer against a "broad" gene set
    -- the top N_HVG_FOR_BROAD_GENE_SET highly-variable genes, not literally
    every gene in the reference's ~38,600-gene annotation. This is both a
    methodological correction and a memory-management necessity: most of
    those ~38,600 genes are near-zero/uninformative in this reference
    (confirmed on the data -- their inclusion measurably *hurt* nearest-
    centroid correlation accuracy relative to a curated/HVG-restricted gene
    set, the opposite of what "full gene set = better ceiling" implies),
    and densifying a 100,000+ cell x 38,600-gene matrix for classification
    caused an out-of-memory kill during development. Both the gene set
    and the test-cell count are bounded for this reason.
    """
    train_idx, test_idx = train_test_split(
        adata.obs_names,
        test_size=test_size,
        random_state=rng_seed,
        stratify=labels.reindex(adata.obs_names),
    )
    if len(test_idx) > MAX_BENCHMARK_TEST_CELLS:
        test_idx, _ = train_test_split(
            test_idx,
            train_size=MAX_BENCHMARK_TEST_CELLS,
            random_state=rng_seed,
            stratify=labels.reindex(test_idx),
        )
    train_adata = adata[train_idx]
    test_adata = adata[test_idx]
    train_labels = labels.reindex(train_idx)
    test_labels = labels.reindex(test_idx)

    hvg_adata = train_adata.copy()
    sc.pp.highly_variable_genes(
        hvg_adata, n_top_genes=min(N_HVG_FOR_BROAD_GENE_SET, hvg_adata.n_vars)
    )
    broad_genes = hvg_adata.var_names[hvg_adata.var["highly_variable"]].tolist()

    test_expr_broad = pd.DataFrame(
        (
            test_adata[:, broad_genes].X.toarray()
            if sparse.issparse(test_adata.X)
            else np.asarray(test_adata[:, broad_genes].X)
        ),
        index=test_adata.obs_names,
        columns=broad_genes,
    )

    broad_centroids = compute_state_centroids(train_adata, train_labels, gene_subset=broad_genes)
    broad_pred = classify_by_nearest_centroid(test_expr_broad, broad_centroids)
    broad_accuracy = float(
        (broad_pred["predicted_state"].to_numpy() == test_labels.to_numpy()).mean()
    )
    broad_mean_confidence = float(broad_pred["confidence"].mean())

    panel_overlap = [g for g in panel_genes if g in adata.var_names]
    test_expr_panel = pd.DataFrame(
        (
            test_adata[:, panel_overlap].X.toarray()
            if sparse.issparse(test_adata.X)
            else np.asarray(test_adata[:, panel_overlap].X)
        ),
        index=test_adata.obs_names,
        columns=panel_overlap,
    )
    panel_centroids = compute_state_centroids(train_adata, train_labels, gene_subset=panel_overlap)
    panel_pred = classify_by_nearest_centroid(test_expr_panel, panel_centroids)
    panel_accuracy = float(
        (panel_pred["predicted_state"].to_numpy() == test_labels.to_numpy()).mean()
    )
    panel_mean_confidence = float(panel_pred["confidence"].mean())

    return {
        "n_train": len(train_idx),
        "n_test": len(test_idx),
        "n_broad_genes": len(broad_genes),
        "n_panel_overlap_genes": len(panel_overlap),
        "broad_gene_set_accuracy": broad_accuracy,
        "broad_gene_set_mean_confidence": broad_mean_confidence,
        "panel_restricted_accuracy": panel_accuracy,
        "panel_restricted_mean_confidence": panel_mean_confidence,
        "accuracy_degradation": broad_accuracy - panel_accuracy,
        "confidence_degradation": broad_mean_confidence - panel_mean_confidence,
    }


def build_reference_mapping_report(project_root: Path) -> dict:
    reference_root = project_root / "data" / "external" / "GSE287301"
    archive_path = reference_root / "GSE287301_filtered_feature_bc_matrix.tar.gz"
    matrix_dir = reference_root / "filtered_feature_bc_matrix"
    matrix_release_path = (
        project_root / "data" / "releases" / "v1_primary_analysis" / "primary_analysis_matrix.h5ad"
    )
    output_path = project_root / "data" / "derived" / "reference_labels.parquet"

    if not matrix_release_path.is_file():
        raise PipelineError(
            f"'{matrix_release_path}' not found. Run `05_preprocessing_and_normalisation/05_create_primary_analysis_matrix.py` first."
        )

    extract_reference_matrix(archive_path, matrix_dir)
    reference = load_reference_matrix(matrix_dir)
    reference = qc_filter_reference(reference)
    reference = normalize_reference(reference)

    t_cell_labels = annotate_t_cell_states(reference)

    xenium = ad.read_h5ad(matrix_release_path)
    panel_genes = list(xenium.var_names)

    degradation = benchmark_transfer_confidence_degradation(reference, t_cell_labels, panel_genes)

    panel_overlap = [g for g in panel_genes if g in reference.var_names]
    centroids = compute_state_centroids(reference, t_cell_labels, gene_subset=panel_overlap)

    xenium_layer = xenium.uns["primary_normalization_layer"]
    xenium_expr = xenium.layers[xenium_layer]
    xenium_expr = xenium_expr.toarray() if sparse.issparse(xenium_expr) else np.asarray(xenium_expr)
    xenium_expr_df = pd.DataFrame(xenium_expr, index=xenium.obs_names, columns=xenium.var_names)[
        panel_overlap
    ]

    predicted = classify_by_nearest_centroid(xenium_expr_df, centroids)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    predicted.to_parquet(output_path)

    return {
        "reference_n_cells_raw": None,
        "reference_n_cells_after_qc": reference.n_obs,
        "reference_n_genes": reference.n_vars,
        "n_t_cell_states": len(T_CELL_STATE_MARKERS),
        "t_cell_state_counts_in_reference": t_cell_labels.value_counts().to_dict(),
        "transfer_confidence_degradation": degradation,
        "n_xenium_cells_labeled": len(predicted),
        "xenium_predicted_state_counts": predicted["predicted_state"].value_counts().to_dict(),
        "xenium_mean_confidence": float(predicted["confidence"].mean()),
    }
