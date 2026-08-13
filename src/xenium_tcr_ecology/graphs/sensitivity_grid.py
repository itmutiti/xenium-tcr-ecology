"""Graph sensitivity grid (`09_spatial_graph_construction_and_calibration/06_run_graph_sensitivity_grid.py`).

Repeats core, already-established graph metrics over all plausible
radii/k values (`09_spatial_graph_construction_and_calibration/00_generate_candidate_spatial_graphs.py`'s candidates, gap-pruned in `09_spatial_graph_construction_and_calibration/01_prune_graphs_for_tissue_gaps.py`) to
check whether conclusions drawn from the calibrated primary graph (Phase
9.02/9.03, radius=30um) would change under a different reasonable
parameter choice -- the concrete purpose a sensitivity grid exists to
serve, not just tabulating raw metrics per candidate in isolation.

**Core metrics, reusing already-established, already-tested functions
directly rather than re-deriving them:**
1. `largest_component_fraction` (`09_spatial_graph_construction_and_calibration/02_calibrate_graph_parameters.py`'s own
   `compute_largest_component_fraction`) -- general graph topology.
2. `fraction_tcells_in_contact` -- the tumour-T-cell contact rate (Phase
   9.04's own `extract_bipartite_subgraph`, generalised here from the two
   chosen tiers to all six radius/k-NN candidates) -- the concrete,
   scientifically load-bearing metric this project's later clone-tumour
   engagement work (`11_clone_spatial_descriptors/02_quantify_clone_tumour_engagement.py`, Clone Ecology Confirmatory Models) will actually depend on, so it is the
   metric most worth stress-testing for parameter sensitivity now, not a
   generic placeholder.

**Robustness check, not just a raw metric table:** for
`fraction_tcells_in_contact`, the per-SECTION values at every candidate
are Spearman-rank-correlated against the calibrated `radius_30.0um`
choice -- a high rank correlation means "which sections have relatively
more/less tumour-T-cell contact" is a conclusion robust to the exact
radius chosen, not an artefact of the specific calibrated value; a low
one would mean this metric's conclusions are sensitive to graph
construction choice and should be reported with that caveat attached.
Only the three radius candidates (15/30/50um) are correlated against the
30um calibrated choice for this purpose -- k-NN graphs are a
qualitatively different construction (density-adaptive neighbour count,
not physical distance), so a k-NN-vs-radius rank correlation would
conflate "sensitive to radius choice" with "sensitive to graph TYPE",
two different questions.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.graphs.graph_calibration import compute_largest_component_fraction
from xenium_tcr_ecology.graphs.tumour_tcell_bipartite import extract_bipartite_subgraph

CANDIDATE_GRAPH_NAMES = [
    "radius_15.0um",
    "radius_30.0um",
    "radius_50.0um",
    "knn_6",
    "knn_10",
    "knn_15",
]
CALIBRATED_RADIUS_GRAPH = "radius_30.0um"
RADIUS_CANDIDATES_FOR_ROBUSTNESS_CHECK = ["radius_15.0um", "radius_30.0um", "radius_50.0um"]


def build_sensitivity_grid(project_root: Path) -> dict:
    pruned_dir = project_root / "data" / "graphs" / "pruned"
    final_annotations_path = project_root / "data" / "derived" / "final_cell_annotations.parquet"
    tumour_masks_dir = project_root / "data" / "derived" / "tumour_masks"
    output_path = project_root / "reports" / "graphs" / "sensitivity_grid.parquet"

    for p in (pruned_dir, final_annotations_path, tumour_masks_dir):
        if not p.exists():
            raise PipelineError(f"'{p}' not found. Run the corresponding earlier phase first.")

    final_annotations = pd.read_parquet(final_annotations_path)

    rows = []
    for section_dir in sorted(d for d in pruned_dir.iterdir() if d.is_dir()):
        section_id = section_dir.name
        mask_path = tumour_masks_dir / f"{section_id}.parquet"
        if not mask_path.is_file():
            raise PipelineError(
                f"'{mask_path}' not found. Run `07_tumour_epithelium_characterisation/04_construct_tumour_region_masks.py` first."
            )

        cell_order = pd.read_csv(section_dir / "cell_order.tsv", sep="\t")["cell_id"]
        tumour_mask_df = pd.read_parquet(mask_path)
        is_malignant = tumour_mask_df.reindex(cell_order)["in_tumour_region"].to_numpy(
            dtype=bool, na_value=False
        )
        is_tcell = (final_annotations.reindex(cell_order)["final_lineage"] == "T_cell").to_numpy()
        n_tcells = int(is_tcell.sum())

        for graph_name in CANDIDATE_GRAPH_NAMES:
            graph_path = section_dir / f"{graph_name}.npz"
            if not graph_path.is_file():
                continue
            graph = sparse.load_npz(graph_path)

            largest_frac = compute_largest_component_fraction(graph)
            bipartite = extract_bipartite_subgraph(graph, is_malignant, is_tcell)
            n_tcells_in_contact = int(
                ((np.asarray(bipartite.sum(axis=1)).ravel() > 0) & is_tcell).sum()
            )

            rows.append(
                {
                    "section_id": section_id,
                    "graph_name": graph_name,
                    "largest_component_fraction": largest_frac,
                    "n_tcells": n_tcells,
                    "n_tcells_in_contact": n_tcells_in_contact,
                    "fraction_tcells_in_contact": (
                        (n_tcells_in_contact / n_tcells) if n_tcells else None
                    ),
                }
            )

    grid = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    grid.to_parquet(output_path)

    calibrated = grid[grid["graph_name"] == CALIBRATED_RADIUS_GRAPH].set_index("section_id")[
        "fraction_tcells_in_contact"
    ]
    robustness: dict[str, float | None] = {}
    for graph_name in RADIUS_CANDIDATES_FOR_ROBUSTNESS_CHECK:
        candidate = grid[grid["graph_name"] == graph_name].set_index("section_id")[
            "fraction_tcells_in_contact"
        ]
        common = calibrated.index.intersection(candidate.index)
        if len(common) < 3:
            robustness[graph_name] = None
            continue
        rho, _ = spearmanr(calibrated.loc[common], candidate.loc[common])
        robustness[graph_name] = round(float(rho), 4)

    return {
        "n_sections": grid["section_id"].nunique(),
        "n_candidate_graph_types": grid["graph_name"].nunique(),
        "n_rows": len(grid),
        "median_fraction_tcells_in_contact_by_graph": grid.groupby("graph_name")[
            "fraction_tcells_in_contact"
        ]
        .median()
        .round(4)
        .to_dict(),
        "spearman_rho_vs_calibrated_radius": robustness,
        "output_path": str(output_path),
    }
