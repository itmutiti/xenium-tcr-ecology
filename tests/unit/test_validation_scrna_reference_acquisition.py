"""Unit tests for xenium_tcr_ecology.validation.scrna_reference_acquisition (`16_external_validation_and_generalisation/02_acquire_hnscc_scrna_references.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.validation.scrna_reference_acquisition import find_til_sample_names


class TestFindTilSampleNames:
    def test_real_only_til_samples_are_returned(self):
        filenames = [
            "GSM1_HNSCC_1_PBMC_barcodes.tsv.gz",
            "GSM2_HNSCC_1_TIL_barcodes.tsv.gz",
            "GSM3_HNSCC_1_TIL_genes.tsv.gz",
            "GSM3_HNSCC_1_TIL_matrix.mtx.gz",
            "GSM4_HD_1_PBMC_barcodes.tsv.gz",
            "GSM5_HD_1_Tonsil_barcodes.tsv.gz",
            "GSM6_HNSCC_2_TIL_barcodes.tsv.gz",
        ]
        result = find_til_sample_names(filenames)
        assert result == {"HNSCC_1_TIL", "HNSCC_2_TIL"}

    def test_real_no_til_samples_gives_empty_set(self):
        filenames = ["GSM1_HD_1_PBMC_barcodes.tsv.gz", "GSM2_HD_2_Tonsil_barcodes.tsv.gz"]
        result = find_til_sample_names(filenames)
        assert result == set()
