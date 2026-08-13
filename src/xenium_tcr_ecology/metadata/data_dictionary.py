"""Compile config/metadata/data_dictionary_input.yaml into
metadata/data_dictionary.xlsx (`01_project_setup_and_governance/02_generate_data_dictionary.py`), one worksheet per table, and
cross-check it against the sample_manifest.tsv columns so the
dictionary cannot silently drift from what `01_project_setup_and_governance/01_build_sample_manifest.py` actually produced.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.validation import validate_records

DICTIONARY_FIELDS = ["field", "unit", "allowed_values", "missingness_code", "derivation_rule"]


def load_dictionary_input(path: Path) -> dict:
    if not path.is_file():
        raise PipelineError(f"Required data-dictionary input not found: '{path}'.")
    return yaml.safe_load(path.read_text())


def _cross_check_against_real_table(
    table_name: str, fields: list[dict], project_root: Path
) -> None:
    """If the table this dictionary block describes already exists on disk,
    confirm every dictionary field is an actual column and every actual
    column is documented -- catching the specific failure mode of a
    dictionary that quietly stops matching its table."""
    candidate_paths = list(project_root.glob(f"metadata/{table_name}")) + list(
        project_root.glob(f"**/{table_name}")
    )
    real_files = [p for p in candidate_paths if p.is_file() and "data/raw" not in str(p)]
    if not real_files:
        return  # table not produced yet -- nothing to cross-check

    real_columns = set(pd.read_csv(real_files[0], sep="\t", nrows=0).columns)
    documented_fields = {f["field"] for f in fields}

    missing_from_dict = real_columns - documented_fields
    missing_from_table = documented_fields - real_columns
    problems = []
    if missing_from_dict:
        problems.append(
            f"columns in '{real_files[0].name}' not documented: {sorted(missing_from_dict)}"
        )
    if missing_from_table:
        problems.append(
            f"documented fields not present in '{real_files[0].name}': {sorted(missing_from_table)}"
        )
    if problems:
        raise PipelineError(
            f"Data dictionary for '{table_name}' is out of sync: " + "; ".join(problems)
        )


def compile_data_dictionary(input_path: Path, output_path: Path, project_root: Path) -> dict:
    data = load_dictionary_input(input_path)
    tables = data["tables"]

    if not tables:
        raise PipelineError(f"'{input_path}' declares no tables.")

    # Validate every table fully before opening the ExcelWriter: a
    # PipelineError raised mid-write would leave an empty workbook whose own
    # __exit__ then raises a second, unrelated openpyxl error ("at least one
    # sheet must be visible"), masking the underlying failure. Fail fast,
    # before any file is touched.
    for table in tables:
        validate_records(
            table["fields"],
            required_fields=DICTIONARY_FIELDS,
            source_path=input_path,
            min_records=1,
            no_placeholder_fields=["derivation_rule"],
        )
        _cross_check_against_real_table(table["table_name"], table["fields"], project_root)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        for table in tables:
            df = pd.DataFrame(table["fields"])[DICTIONARY_FIELDS]
            sheet_name = table["table_name"].replace(".tsv", "")[:31]  # Excel sheet name limit
            df.to_excel(writer, sheet_name=sheet_name, index=False)

    return {"tables_documented": len(tables), "total_fields": sum(len(t["fields"]) for t in tables)}
