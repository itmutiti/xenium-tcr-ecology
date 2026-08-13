"""Unit tests for xenium_tcr_ecology.external_checkpoint.bulk_reference (`12_external_checkpoint_validation/00_project_provisional_signatures_to_bulk_reference.py`)."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from xenium_tcr_ecology.external_checkpoint.bulk_reference import (
    assign_t_cell_substate,
    compare_state_proportions,
    compute_module_score,
    parse_gse103322_full_matrix_t_cells,
    parse_gse103322_t_cells,
    required_gse103322_genes,
)


class TestRequiredGse103322Genes:
    def test_includes_expected_marker_categories(self):
        genes = required_gse103322_genes()
        assert "GZMB" in genes  # cytotoxicity
        assert "PDCD1" in genes  # exhaustion
        assert "MKI67" in genes  # proliferation
        assert "FOXP3" in genes  # treg
        assert "CD4" in genes and "CD8A" in genes  # lineage


class TestComputeModuleScore:
    def test_high_expression_gives_positive_score(self):
        expr = pd.DataFrame({"GENE_A": [10.0, 0.0, 5.0], "GENE_B": [10.0, 0.0, 5.0]})
        result = compute_module_score(expr, ["GENE_A", "GENE_B"])
        assert result.iloc[0] > result.iloc[1]

    def test_zero_variance_gene_does_not_crash(self):
        expr = pd.DataFrame({"GENE_A": [5.0, 5.0, 5.0], "GENE_B": [1.0, 2.0, 3.0]})
        result = compute_module_score(expr, ["GENE_A", "GENE_B"])
        assert len(result) == 3


class TestAssignTCellSubstate:
    def test_treg_takes_highest_priority(self):
        # Treg-positive AND cycling-positive -- Treg must win.
        state = assign_t_cell_substate(
            treg_score=np.array([1.0]),
            cycling_score=np.array([1.0]),
            exhaustion_score=np.array([0.0]),
            cytotoxicity_score=np.array([0.0]),
            cd4_expr=np.array([5.0]),
            cd8a_expr=np.array([1.0]),
        )
        assert state[0] == "Treg"

    def test_cycling_beats_exhausted_and_cytotoxic(self):
        state = assign_t_cell_substate(
            treg_score=np.array([-1.0]),
            cycling_score=np.array([1.0]),
            exhaustion_score=np.array([1.0]),
            cytotoxicity_score=np.array([1.0]),
            cd4_expr=np.array([0.0]),
            cd8a_expr=np.array([0.0]),
        )
        assert state[0] == "Cycling"

    def test_cd8_dominance_without_other_programs(self):
        state = assign_t_cell_substate(
            treg_score=np.array([-1.0]),
            cycling_score=np.array([-1.0]),
            exhaustion_score=np.array([-1.0]),
            cytotoxicity_score=np.array([-1.0]),
            cd4_expr=np.array([1.0]),
            cd8a_expr=np.array([5.0]),
        )
        assert state[0] == "CD8"

    def test_no_positive_program_gives_ambiguous(self):
        state = assign_t_cell_substate(
            treg_score=np.array([-1.0]),
            cycling_score=np.array([-1.0]),
            exhaustion_score=np.array([-1.0]),
            cytotoxicity_score=np.array([-1.0]),
            cd4_expr=np.array([0.0]),
            cd8a_expr=np.array([0.0]),
        )
        assert state[0] == "Ambiguous"


class TestCompareStateProportions:
    def test_real_side_by_side_fractions(self):
        reference = pd.Series(["Cycling", "Cycling", "Exhausted"])
        project = pd.Series(["Cycling", "Exhausted", "Exhausted", "Exhausted"])
        result = compare_state_proportions(reference, project).set_index("state")
        assert result.loc["Cycling", "reference_fraction"] == pytest.approx(2 / 3)
        assert result.loc["Exhausted", "project_fraction"] == pytest.approx(3 / 4)

    def test_state_only_in_one_side_gets_zero_on_the_other(self):
        reference = pd.Series(["Treg"])
        project = pd.Series(["Cycling"])
        result = compare_state_proportions(reference, project).set_index("state")
        assert result.loc["Treg", "project_fraction"] == 0.0
        assert result.loc["Cycling", "reference_fraction"] == 0.0


class TestParseGse103322TCells:
    def test_parses_a_synthetic_file_matching_the_real_format(self, tmp_path):
        genes = sorted(required_gse103322_genes())
        header = "\t" + "\t".join(["cellA", "cellB", "cellC"]) + "\n"
        lines = [header]
        lines.append("processed by Maxima enzyme\t1\t1\t1\n")
        lines.append("Lymph node\t0\t1\t0\n")
        lines.append("classified  as cancer cell\t0\t0\t1\n")
        lines.append("classified as non-cancer cells\t1\t1\t0\n")
        lines.append("non-cancer cell type\tT cell\tFibroblast\t0\n")
        for gene in genes:
            values = "\t".join(["1.5", "2.5", "3.5"])
            lines.append(f"'{gene}'\t{values}\n")

        path = tmp_path / "synthetic.txt.gz"
        with gzip.open(path, "wt") as f:
            f.writelines(lines)

        result = parse_gse103322_t_cells(path)
        assert list(result.index) == ["cellA"]  # only the real T cell
        assert result.loc["cellA", "GZMB"] == 1.5

    def test_missing_required_gene_raises(self, tmp_path):
        header = "\t" + "\t".join(["cellA"]) + "\n"
        lines = [header]
        lines.append("processed by Maxima enzyme\t1\n")
        lines.append("Lymph node\t0\n")
        lines.append("classified  as cancer cell\t0\n")
        lines.append("classified as non-cancer cells\t1\n")
        lines.append("non-cancer cell type\tT cell\n")
        lines.append("'GZMB'\t1.0\n")  # only one of many required genes

        path = tmp_path / "incomplete.txt.gz"
        with gzip.open(path, "wt") as f:
            f.writelines(lines)

        with pytest.raises(Exception):
            parse_gse103322_t_cells(path)


class TestParseGse103322FullMatrixTCells:
    def test_parses_every_gene_row_not_just_markers(self, tmp_path):
        header = "\t" + "\t".join(["cellA", "cellB", "cellC"]) + "\n"
        lines = [header]
        lines.append("processed by Maxima enzyme\t1\t1\t1\n")
        lines.append("Lymph node\t0\t1\t0\n")
        lines.append("classified  as cancer cell\t0\t0\t1\n")
        lines.append("classified as non-cancer cells\t1\t1\t0\n")
        lines.append("non-cancer cell type\tT cell\tFibroblast\t0\n")
        lines.append("'GZMB'\t1.5\t2.5\t3.5\n")
        lines.append("'SOME_UNRELATED_GENE'\t9.0\t8.0\t7.0\n")  # not a required marker

        path = tmp_path / "synthetic_full.txt.gz"
        with gzip.open(path, "wt") as f:
            f.writelines(lines)

        result = parse_gse103322_full_matrix_t_cells(path)
        assert list(result.index) == ["cellA"]  # only the real T cell
        assert result.loc["cellA", "GZMB"] == 1.5
        assert result.loc["cellA", "SOME_UNRELATED_GENE"] == 9.0  # non-marker gene retained

    def test_missing_metadata_row_raises(self, tmp_path):
        header = "\t" + "\t".join(["cellA"]) + "\n"
        lines = [header]
        lines.append("processed by Maxima enzyme\t1\n")
        lines.append("Lymph node\t0\n")
        lines.append("classified  as cancer cell\t0\n")
        lines.append("classified as non-cancer cells\t1\n")
        lines.append("some other row, not the expected metadata row\t1\n")
        lines.append("'GZMB'\t1.0\n")

        path = tmp_path / "malformed_full.txt.gz"
        with gzip.open(path, "wt") as f:
            f.writelines(lines)

        with pytest.raises(Exception):
            parse_gse103322_full_matrix_t_cells(path)
