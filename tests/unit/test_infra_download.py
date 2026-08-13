"""Unit tests for xenium_tcr_ecology.infra.download."""

from __future__ import annotations

import hashlib
from unittest.mock import MagicMock, patch

import pytest

from xenium_tcr_ecology.infra.download import (
    DEFAULT_HEADERS,
    download_file,
    verify_checksums,
)
from xenium_tcr_ecology.infra.exceptions import PipelineError


def _mock_response(content: bytes, status_code: int = 200):
    resp = MagicMock()
    resp.status_code = status_code
    resp.iter_content = lambda chunk_size: [content]
    if status_code >= 400:
        import requests

        resp.raise_for_status.side_effect = requests.HTTPError(f"{status_code} error")
    else:
        resp.raise_for_status.return_value = None
    resp.__enter__ = lambda self: resp
    resp.__exit__ = lambda self, *a: None
    return resp


class TestDownloadFile:
    def test_real_skips_download_when_size_already_matches(self, tmp_path):
        dest = tmp_path / "already_here.txt"
        dest.write_bytes(b"exact content")
        with patch("xenium_tcr_ecology.infra.download.requests.get") as mock_get:
            download_file(
                "https://example.invalid/f", dest, expected_size_bytes=len(b"exact content")
            )
        mock_get.assert_not_called()

    def test_real_downloads_when_missing_and_sends_default_headers(self, tmp_path):
        dest = tmp_path / "new_file.txt"
        content = b"downloaded content"
        with patch("xenium_tcr_ecology.infra.download.requests.get") as mock_get:
            mock_get.return_value = _mock_response(content)
            download_file("https://example.invalid/f", dest, expected_size_bytes=len(content))
        assert dest.read_bytes() == content
        _, kwargs = mock_get.call_args
        assert kwargs["headers"]["User-Agent"] == DEFAULT_HEADERS["User-Agent"]
        assert kwargs["headers"]["Accept-Encoding"] == "identity"

    def test_real_size_mismatch_after_download_raises_and_removes_file(self, tmp_path):
        dest = tmp_path / "truncated.txt"
        with patch("xenium_tcr_ecology.infra.download.requests.get") as mock_get:
            mock_get.return_value = _mock_response(b"short")
            with pytest.raises(PipelineError, match="does not match expected"):
                download_file("https://example.invalid/f", dest, expected_size_bytes=9999)
        assert not dest.exists()

    def test_real_retries_then_raises_pipeline_error_on_persistent_failure(self, tmp_path):
        import requests

        dest = tmp_path / "never.txt"
        with (
            patch("xenium_tcr_ecology.infra.download.requests.get") as mock_get,
            patch("xenium_tcr_ecology.infra.download.time.sleep"),
        ):
            mock_get.side_effect = requests.ConnectionError("real network failure")
            with pytest.raises(PipelineError, match="Failed to download"):
                download_file("https://example.invalid/f", dest)
        assert mock_get.call_count == 5

    def test_real_no_expected_size_still_downloads_and_writes_content(self, tmp_path):
        dest = tmp_path / "no_size_check.txt"
        content = b"whatever length"
        with patch("xenium_tcr_ecology.infra.download.requests.get") as mock_get:
            mock_get.return_value = _mock_response(content)
            download_file("https://example.invalid/f", dest)
        assert dest.read_bytes() == content


class TestVerifyChecksums:
    def test_real_matching_checksum_passes(self, tmp_path):
        content = b"real content"
        (tmp_path / "f.txt").write_bytes(content)
        (tmp_path / "checksums.sha256").write_text(
            f"{hashlib.sha256(content).hexdigest()}  f.txt\n"
        )
        assert verify_checksums(tmp_path) == {"f.txt": True}

    def test_real_mismatched_checksum_fails(self, tmp_path):
        (tmp_path / "f.txt").write_bytes(b"corrupted")
        (tmp_path / "checksums.sha256").write_text(f"{'0' * 64}  f.txt\n")
        assert verify_checksums(tmp_path) == {"f.txt": False}
