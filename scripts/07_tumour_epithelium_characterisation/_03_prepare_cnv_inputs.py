#!/usr/bin/env python3
"""
`07_tumour_epithelium_characterisation/03_infer_cnv_appendix_only.R` helper -- _03_prepare_cnv_inputs.py

NOT its own numbered blueprint phase step: exports the epithelial cell
expression matrix and per-patient non-epithelial reference baseline that
03_infer_cnv_appendix_only.R needs, since that R script cannot read
.h5ad directly (no R HDF5/AnnData reader available). Also resolves gene
genomic coordinates via Ensembl's REST API (cached to
references/gene_genomic_coordinates.tsv). Invoked by the R script via a
subprocess call, not run standalone as a pipeline phase.

Primary output: data/derived/cnv_epithelial_expression.parquet,
data/derived/cnv_reference_baseline.parquet
"""

from __future__ import annotations

import sys

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.tumour.cnv_inference import prepare_cnv_inputs


def main() -> int:
    parser = base_parser(__doc__)
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "07_tumour_epithelium_characterisation",
        script_name="_03_prepare_cnv_inputs",
        project_root=project_root,
        phase="07_tumour_epithelium_characterisation",
    )

    try:
        summary = prepare_cnv_inputs(project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        logger.log_error(str(exc))
        logger.write(status="failed")
        return 1

    logger.log_event(**summary)
    logger.write(status="ok")
    print(
        f"[OK]   {summary['n_genes_used']} gene(s), {summary['n_epithelial_cells']} epithelial cell(s), "
        f"{summary['n_reference_cells']} reference cell(s) across "
        f"{summary['n_patients_with_reference_baseline']} patient(s). "
        f"Wrote {summary['epithelial_expr_path']}, {summary['reference_baseline_path']}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
