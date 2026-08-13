"""Assembles publication-ready sample, QC, model and validation tables
from already-computed source files (`17_statistical_closure_and_release/04_generate_results_tables.py`) -- every source
is copied/converted as-is (parquet sources are converted to TSV for
direct human readability; TSV sources are copied unchanged), no values
are recomputed or edited.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

TABLE_MANIFEST: list[dict] = [
    {
        "table_number": 1,
        "name": "sample_manifest",
        "source_path": "metadata/sample_manifest.tsv",
        "category": "sample",
    },
    {
        "table_number": 2,
        "name": "replicate_concordance",
        "source_path": "data/derived/replicate_concordance.tsv",
        "category": "QC",
    },
    {
        "table_number": 3,
        "name": "statistical_summary",
        "source_path": "results/statistical_summary.tsv",
        "category": "model",
    },
    {
        "table_number": 4,
        "name": "variance_partition_results",
        "source_path": "data/derived/variance_partition_results.parquet",
        "category": "model",
    },
    {
        "table_number": 5,
        "name": "barrier_topology_model_results",
        "source_path": "data/derived/barrier_topology_model_results.parquet",
        "category": "model",
    },
    {
        "table_number": 6,
        "name": "clone_structure_test_results",
        "source_path": "data/derived/clone_structure_test_results.parquet",
        "category": "model",
    },
    {
        "table_number": 7,
        "name": "claim_evidence_matrix",
        "source_path": "results/claim_evidence_matrix.tsv",
        "category": "validation",
    },
    {
        "table_number": 8,
        "name": "source_paper_comparison",
        "source_path": "reports/validation/source_paper_comparison.tsv",
        "category": "validation",
    },
    {
        "table_number": 9,
        "name": "validation_plan",
        "source_path": "governance/validation_plan.tsv",
        "category": "validation",
    },
    {
        "table_number": 10,
        "name": "hpv_claim_strength",
        "source_path": "results/hpv_claim_strength.tsv",
        "category": "validation",
    },
]


def build_results_tables(project_root: Path) -> dict:
    tables_dir = project_root / "tables"

    missing = [
        entry["source_path"]
        for entry in TABLE_MANIFEST
        if not (project_root / entry["source_path"]).is_file()
    ]
    if missing:
        raise PipelineError(f"Missing source table(s), cannot assemble results tables: {missing}")

    tables_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for entry in TABLE_MANIFEST:
        source = project_root / entry["source_path"]
        dest_name = f"Table_{entry['table_number']}_{entry['name']}.tsv"
        dest = tables_dir / dest_name
        if source.suffix == ".parquet":
            pd.read_parquet(source).to_csv(dest, sep="\t", index=False)
        else:
            dest.write_bytes(source.read_bytes())
        manifest_rows.append(
            {
                "table_number": entry["table_number"],
                "filename": dest_name,
                "category": entry["category"],
                "source_path": entry["source_path"],
            }
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = tables_dir / "MANIFEST.tsv"
    manifest_df.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_tables": len(manifest_rows),
        "tables_dir": str(tables_dir),
        "manifest_path": str(manifest_path),
    }
