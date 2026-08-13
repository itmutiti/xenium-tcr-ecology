"""Unit tests for xenium_tcr_ecology.clone_ecology.taxonomy_release (`13_clone_ecology_confirmatory_models/00_load_frozen_taxonomy_version.py`)."""

from __future__ import annotations

from xenium_tcr_ecology.clone_ecology.taxonomy_release import (
    check_hash_consistency,
    compute_file_hash,
)


class TestComputeFileHash:
    def test_real_file_hash_is_deterministic(self, tmp_path):
        path = tmp_path / "a.txt"
        path.write_text("hello world")
        first = compute_file_hash(path)
        second = compute_file_hash(path)
        assert first == second
        assert len(first) == 64  # real SHA256 hex digest length

    def test_different_content_gives_different_hash(self, tmp_path):
        path_a = tmp_path / "a.txt"
        path_b = tmp_path / "b.txt"
        path_a.write_text("hello")
        path_b.write_text("world")
        assert compute_file_hash(path_a) != compute_file_hash(path_b)


class TestCheckHashConsistency:
    def test_no_changes_gives_empty_list(self):
        previous = {"a.parquet": "hash1", "b.parquet": "hash2"}
        current = {"a.parquet": "hash1", "b.parquet": "hash2"}
        assert check_hash_consistency(previous, current) == []

    def test_changed_file_is_detected(self):
        previous = {"a.parquet": "hash1", "b.parquet": "hash2"}
        current = {"a.parquet": "hash1", "b.parquet": "DIFFERENT"}
        assert check_hash_consistency(previous, current) == ["b.parquet"]

    def test_multiple_changed_files_all_reported_sorted(self):
        previous = {"z.parquet": "h1", "a.parquet": "h2"}
        current = {"z.parquet": "CHANGED_Z", "a.parquet": "CHANGED_A"}
        assert check_hash_consistency(previous, current) == ["a.parquet", "z.parquet"]

    def test_file_missing_from_current_is_not_flagged_as_changed(self):
        # A file present before but absent now is a different real
        # condition (file removed) from a real content change -- this
        # function only reports genuine hash mismatches for files
        # present in both.
        previous = {"a.parquet": "hash1"}
        current = {}
        assert check_hash_consistency(previous, current) == []
