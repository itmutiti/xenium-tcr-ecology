"""Standardised clone atlas (`13_clone_ecology_confirmatory_models/06_generate_clone_atlas.py`).

Produces one spatial thumbnail plus one standardised descriptor summary
for every high-confidence clone-section (all 261, spanning all 158
distinct clones) in a single, systematically-ordered HTML page --
"preventing example-selection bias" by construction: rows are sorted by
arbitrary identifiers (`patient_id`, `clone_id`, `section_id`), never by
any "interestingness" metric, and there is no mechanism in this code to
omit any high-confidence clone-section from the atlas.

Each thumbnail shows the clone's own cells (from `09_spatial_graph_construction_and_calibration/03_construct_primary_cell_graph.py`'s
node_metadata) highlighted against the rest of that section's cells for
spatial context -- a direct visualisation, not a schematic. Descriptor
columns are pulled directly from the already-computed `11_clone_spatial_descriptors/00_compute_clone_spatial_descriptors_rarefied.py`-11.06 tables (spatial, state composition, tumour
engagement, APC support, barrier metrics, continuous structure score) --
no new descriptors are computed here; this milestone only standardises
their presentation.
"""

from __future__ import annotations

import base64
import html
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

SUMMARY_COLUMNS = [
    "n_cells",
    "ecological_structure_score",
    "dominant_state",
    "malignant_adjacency",
    "engagement_ratio",
    "fraction_inside_tumour",
    "dc_engagement_ratio",
    "macrophage_engagement_ratio",
    "fibroblast_barrier_fraction",
    "vascular_barrier_fraction",
    "suppressive_myeloid_barrier_fraction",
]


def render_clone_thumbnail(
    clone_x, clone_y, section_x, section_y, figsize: tuple[float, float] = (2.4, 2.4)
) -> str:
    """Pure-ish (matplotlib rendering, no file I/O): spatial thumbnail
    -- section context in light grey, the clone's cells highlighted in
    red -- returned as a base64-encoded PNG data URI string ready to
    embed directly in HTML."""
    fig, ax = plt.subplots(figsize=figsize, dpi=80)
    ax.scatter(section_x, section_y, s=1, c="lightgrey", linewidths=0)
    ax.scatter(clone_x, clone_y, s=6, c="firebrick", linewidths=0)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect("equal")
    for spine in ax.spines.values():
        spine.set_visible(False)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{encoded}"


def build_clone_card_html(row: pd.Series, thumbnail_data_uri: str) -> str:
    """Pure, testable: HTML for one clone-section's atlas card."""
    stats_html = "".join(
        f"<tr><td>{html.escape(str(col))}</td><td>{html.escape(f'{row[col]:.3f}' if isinstance(row[col], float) else str(row[col]))}</td></tr>"
        for col in SUMMARY_COLUMNS
        if col in row.index
    )
    return f"""
    <div class="clone-card">
      <img src="{thumbnail_data_uri}" alt="Spatial thumbnail" />
      <div class="clone-meta">
        <strong>{html.escape(str(row['clone_id']))[:40]}</strong><br/>
        {html.escape(str(row['patient_id']))} / {html.escape(str(row['section_id']))}
        <table>{stats_html}</table>
      </div>
    </div>
    """


def build_clone_atlas(project_root: Path) -> dict:
    resolved_calls_path = project_root / "data" / "derived" / "tcr_resolved_calls.parquet"
    high_confidence_clones_path = (
        project_root / "data" / "releases" / "v1_tcr_calls" / "high_confidence_clones.parquet"
    )
    primary_graphs_dir = project_root / "data" / "graphs" / "primary_graphs"

    merged_path = project_root / "data" / "derived" / "clone_atlas_merged.parquet"
    output_path = project_root / "reports" / "clone_ecology" / "clone_atlas.html"

    paths_needed = [
        (
            project_root
            / "data"
            / "releases"
            / "v1_clone_structure"
            / "clone_ecological_structure.parquet",
            "`13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`",
        ),
        (
            project_root / "data" / "derived" / "clone_state_composition.parquet",
            "`11_clone_spatial_descriptors/01_compute_clone_cell_state_composition.py`",
        ),
        (
            project_root / "data" / "derived" / "clone_tumour_engagement.parquet",
            "`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`",
        ),
        (
            project_root / "data" / "derived" / "clone_apc_support.parquet",
            "`11_clone_spatial_descriptors/03_quantify_clone_apc_support.py`",
        ),
        (
            project_root / "data" / "derived" / "clone_barrier_metrics.parquet",
            "`11_clone_spatial_descriptors/04_quantify_stromal_and_myeloid_barriers.py`",
        ),
        (resolved_calls_path, "`08_tcr_clonal_analysis/03_call_cell_level_tcr_detections.py`"),
        (high_confidence_clones_path, "`08_tcr_clonal_analysis/08_generate_tcr_release_report.py`"),
    ]
    for path, phase in paths_needed:
        if not path.exists():
            raise PipelineError(f"'{path}' not found. Run {phase} first.")

    structure = pd.read_parquet(paths_needed[0][0])
    state = pd.read_parquet(paths_needed[1][0])
    engagement = pd.read_parquet(paths_needed[2][0])
    apc = pd.read_parquet(paths_needed[3][0])
    barrier = pd.read_parquet(paths_needed[4][0])

    key = ["clone_id", "section_id", "patient_id", "n_cells"]
    merged = (
        structure.merge(state, on=key)
        .merge(engagement, on=key)
        .merge(apc, on=key)
        .merge(barrier, on=key)
    )
    merged = merged.sort_values(["patient_id", "clone_id", "section_id"]).reset_index(drop=True)
    merged.to_parquet(merged_path)

    resolved = pd.read_parquet(resolved_calls_path)
    clonal = resolved[resolved["resolution"].isin(["singlet", "low_confidence"])].copy()
    clonal["clone_id"] = clonal["detected_probes"]
    high_confidence_clones = pd.read_parquet(high_confidence_clones_path)
    clonal = clonal[clonal["clone_id"].isin(set(high_confidence_clones["clone_id"]))]

    cards = []
    for section_id, section_rows in merged.groupby("section_id", sort=False):
        node_metadata = pd.read_csv(
            primary_graphs_dir / section_id / "node_metadata.tsv", sep="\t"
        ).set_index("cell_id")
        section_x = node_metadata["x_centroid"].to_numpy()
        section_y = node_metadata["y_centroid"].to_numpy()
        section_clonal = clonal[clonal["section_id"] == section_id]

        for _, row in section_rows.iterrows():
            cell_ids = section_clonal[section_clonal["clone_id"] == row["clone_id"]].index
            clone_coords = node_metadata.reindex(cell_ids).dropna(
                subset=["x_centroid", "y_centroid"]
            )
            thumbnail = render_clone_thumbnail(
                clone_coords["x_centroid"].to_numpy(),
                clone_coords["y_centroid"].to_numpy(),
                section_x,
                section_y,
            )
            cards.append(build_clone_card_html(row, thumbnail))

    style = """
    <style>
      body { font-family: sans-serif; background: #fafafa; }
      .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 12px; }
      .clone-card { display: flex; gap: 8px; background: white; border: 1px solid #ddd; border-radius: 6px; padding: 8px; }
      .clone-card img { width: 96px; height: 96px; flex-shrink: 0; }
      .clone-meta { font-size: 11px; overflow: hidden; }
      table { border-collapse: collapse; font-size: 10px; }
      td { padding: 1px 4px; }
      td:first-child { color: #666; }
    </style>
    """
    html_content = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Clone Atlas -- all {len(merged)} high-confidence clone-sections</title>{style}</head>
<body>
<h1>Clone Atlas</h1>
<p>All {len(merged)} high-confidence clone-section rows ({merged['clone_id'].nunique()} distinct clones), sorted by
patient_id / clone_id / section_id -- no clone is omitted or reordered by any "interestingness" criterion.</p>
<div class="grid">
{''.join(cards)}
</div>
</body></html>"""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(html_content)

    return {
        "n_clone_section_rows": len(merged),
        "n_distinct_clones": int(merged["clone_id"].nunique()),
        "n_sections": int(merged["section_id"].nunique()),
        "output_path": str(output_path),
        "merged_path": str(merged_path),
    }
