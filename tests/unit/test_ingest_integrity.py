"""Unit tests for xenium_tcr_ecology.ingest.integrity (`02_raw_data_ingestion/02_verify_archive_checksums.py`), using
small synthetic tar fixtures -- doesn't require the real 52 GB download."""

from __future__ import annotations

import json
import tarfile

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.ingest.integrity import sha256_of_file, verify_archive_and_checksum


def _make_tar(tmp_path, files: dict[str, bytes]) -> "Path":  # noqa: F821
    tar_path = tmp_path / "test.tar"
    with tarfile.open(tar_path, "w") as tf:
        for name, content in files.items():
            member_path = tmp_path / name
            member_path.write_bytes(content)
            tf.add(member_path, arcname=name)
    return tar_path


def _write_snapshot(tmp_path, entries: list[dict]):
    path = tmp_path / "snapshot.json"
    path.write_text(json.dumps({"file_entries": entries}))
    return path


class TestSha256OfFile:
    def test_known_hash(self, tmp_path):
        p = tmp_path / "f.txt"
        p.write_bytes(b"hello world")
        assert (
            sha256_of_file(p)
            == "b94d27b9934d3e08a52e52d7da7dabfac484efe37a5380ee9088f7ace2efcde9"  # pragma: allowlist secret
        )


class TestVerifyArchiveAndChecksum:
    def test_passes_on_matching_archive(self, tmp_path):
        tar_path = _make_tar(tmp_path, {"a.txt": b"hello world", "b.txt": b"goodbye"})
        snapshot = _write_snapshot(
            tmp_path,
            [
                {"name": "a.txt", "size_bytes": 11},
                {"name": "b.txt", "size_bytes": 7},
            ],
        )
        summary = verify_archive_and_checksum(
            tar_path, snapshot, tmp_path / "SHA256SUMS", tmp_path / "integrity.tsv", tmp_path
        )
        assert summary["files_verified"] == 2
        assert summary["files_expected"] == 2
        assert (tmp_path / "SHA256SUMS").is_file()
        assert sha256_of_file(tar_path) == summary["archive_sha256"]

    def test_raises_on_size_mismatch(self, tmp_path):
        tar_path = _make_tar(tmp_path, {"a.txt": b"hello world"})
        snapshot = _write_snapshot(tmp_path, [{"name": "a.txt", "size_bytes": 999}])
        with pytest.raises(PipelineError, match="Archive integrity check failed"):
            verify_archive_and_checksum(
                tar_path, snapshot, tmp_path / "SHA256SUMS", tmp_path / "integrity.tsv", tmp_path
            )
        report = (tmp_path / "integrity.tsv").read_text()
        assert "size_mismatch" in report

    def test_raises_on_missing_expected_file(self, tmp_path):
        tar_path = _make_tar(tmp_path, {"a.txt": b"hello world"})
        snapshot = _write_snapshot(
            tmp_path,
            [
                {"name": "a.txt", "size_bytes": 11},
                {"name": "missing.txt", "size_bytes": 5},
            ],
        )
        with pytest.raises(PipelineError, match="Archive integrity check failed"):
            verify_archive_and_checksum(
                tar_path, snapshot, tmp_path / "SHA256SUMS", tmp_path / "integrity.tsv", tmp_path
            )
        report = (tmp_path / "integrity.tsv").read_text()
        assert "missing_from_archive" in report

    def test_raises_on_unexpected_file_in_archive(self, tmp_path):
        tar_path = _make_tar(tmp_path, {"a.txt": b"hello world", "extra.txt": b"surprise"})
        snapshot = _write_snapshot(tmp_path, [{"name": "a.txt", "size_bytes": 11}])
        with pytest.raises(PipelineError, match="Archive integrity check failed"):
            verify_archive_and_checksum(
                tar_path, snapshot, tmp_path / "SHA256SUMS", tmp_path / "integrity.tsv", tmp_path
            )
        report = (tmp_path / "integrity.tsv").read_text()
        assert "unexpected_in_archive" in report

    def test_raises_on_missing_archive(self, tmp_path):
        snapshot = _write_snapshot(tmp_path, [{"name": "a.txt", "size_bytes": 11}])
        with pytest.raises(PipelineError, match="not found"):
            verify_archive_and_checksum(
                tmp_path / "nope.tar",
                snapshot,
                tmp_path / "SHA256SUMS",
                tmp_path / "integrity.tsv",
                tmp_path,
            )
