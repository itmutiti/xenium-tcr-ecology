#!/usr/bin/env python3
"""
`14_spatial_interactions_and_barriers/02_compute_spatially_constrained_scores.py`

Scores panel-complete, program-overlap-restricted ligand-receptor pairs
(`14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py` sender-receiver pairs x `14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py` LR pairs) across
graph-connected sender->receiver cell pairs in each primary section,
against a degree-preserving null (reusing `09_spatial_graph_construction_and_calibration/08_run_calibration_suite_on_synthetic_data.py`'s
calibrated degree-stratification methodology). Exploratory -- not a
prespecified confirmatory analysis (see governance/analysis_registry.tsv;
only `14_spatial_interactions_and_barriers/03_model_barrier_topology_by_structure.R`'s barrier-topology model is prespecified
as q3_barrier_topology_confirmatory).

See src/xenium_tcr_ecology/interactions/spatial_scores.py's module
docstring for the
method and findings.

Primary output: data/derived/spatial_interaction_scores.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.interactions.spatial_scores import build_spatial_interaction_scores


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "14_spatial_interactions_and_barriers",
        script_name="02_compute_spatially_constrained_scores",
        project_root=project_root,
        phase="14_spatial_interactions_and_barriers",
    )

    try:
        summary = build_spatial_interaction_scores(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_combinations_tested']} program-overlap combination(s) tested across "
        f"{summary['n_sections']} section(s): {summary['n_rows']} row(s), {summary['n_significant']} significant "
        f"(p<0.05). Wrote {summary['output_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
