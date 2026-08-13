"""Unit tests for xenium_tcr_ecology.ingest.extract (`02_raw_data_ingestion/03_extract_archive_safely.py`), including
the path-traversal defence -- the single most important property of this
module to actually verify, not just assert exists in a docstring."""

from __future__ import annotations

import tarfile

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.ingest.extract import safe_extract


def _make_tar_with_raw_names(tmp_path, files: dict[str, bytes]):
    """Build a tar where member names are controlled exactly as given
    (unlike tarfile.add(), which derives the name from a real path) -- lets
    us construct a malicious ../ member name for the traversal test."""
    tar_path = tmp_path / "archive.tar"
    with tarfile.open(tar_path, "w") as tf:
        for name, content in files.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(content)
            import io

            tf.addfile(info, io.BytesIO(content))
    return tar_path


class TestSafeExtract:
    def test_extracts_grouped_by_gsm(self, tmp_path):
        tar_path = _make_tar_with_raw_names(
            tmp_path,
            {
                "GSM123_foo_cells.parquet.gz": b"data1",
                "GSM123_foo_cell_feature_matrix.h5": b"data2",
                "GSM456_bar_cells.parquet.gz": b"data3",
            },
        )
        dest_root = tmp_path / "staged"
        summary = safe_extract(tar_path, dest_root, project_root=tmp_path)

        assert summary["files_extracted"] == 3
        assert summary["samples"] == 2
        assert (dest_root / "GSM123" / "GSM123_foo_cells.parquet.gz").read_bytes() == b"data1"
        assert (dest_root / "GSM456" / "GSM456_bar_cells.parquet.gz").read_bytes() == b"data3"

    def test_sample_count_ignores_non_directory_entries(self, tmp_path):
        """Regression test: a stray non-directory file already present in
        dest_root (e.g. the scaffold's .gitkeep placeholder) must not
        inflate the reported sample count -- caught for real during the
        Raw Data Ingestion run against the actual 18-sample GEO archive, which
        initially (incorrectly) reported 19 samples because of exactly
        this."""
        tar_path = _make_tar_with_raw_names(tmp_path, {"GSM123_foo_cells.parquet.gz": b"data1"})
        dest_root = tmp_path / "staged"
        dest_root.mkdir()
        (dest_root / ".gitkeep").write_bytes(b"")

        summary = safe_extract(tar_path, dest_root, project_root=tmp_path)
        assert summary["samples"] == 1

    def test_rejects_path_traversal_member(self, tmp_path):
        tar_path = _make_tar_with_raw_names(
            tmp_path,
            {
                "../../etc/GSM999_evil.txt": b"malicious",
            },
        )
        with pytest.raises(PipelineError, match="path-traversal"):
            safe_extract(tar_path, tmp_path / "staged", project_root=tmp_path)

    def test_rejects_absolute_path_member(self, tmp_path):
        tar_path = _make_tar_with_raw_names(
            tmp_path,
            {
                "/etc/GSM999_evil.txt": b"malicious",
            },
        )
        with pytest.raises(PipelineError, match="absolute path"):
            safe_extract(tar_path, tmp_path / "staged", project_root=tmp_path)

    def test_traversal_member_never_reaches_disk_outside_dest(self, tmp_path):
        """Belt-and-braces: even if the raised-exception guard were somehow
        removed, destination paths are built from the basename only, so a
        traversal member could not escape dest_root. Verify this directly by
        checking no file was written outside dest_root before the guard
        raised."""
        tar_path = _make_tar_with_raw_names(
            tmp_path,
            {
                "../escaped.txt": b"malicious",
            },
        )
        canary = tmp_path / "escaped.txt"
        with pytest.raises(PipelineError):
            safe_extract(tar_path, tmp_path / "staged", project_root=tmp_path)
        assert not canary.exists()

    def test_raises_on_unrecognised_filename_pattern(self, tmp_path):
        tar_path = _make_tar_with_raw_names(tmp_path, {"not_a_gsm_file.txt": b"x"})
        with pytest.raises(PipelineError, match="Could not determine GSM"):
            safe_extract(tar_path, tmp_path / "staged", project_root=tmp_path)

    def test_raises_on_missing_archive(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            safe_extract(tmp_path / "nope.tar", tmp_path / "staged", project_root=tmp_path)
