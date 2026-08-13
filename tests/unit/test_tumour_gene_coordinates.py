"""Unit tests for xenium_tcr_ecology.tumour.gene_coordinates."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.tumour.gene_coordinates import fetch_gene_coordinates


class _FakeResponse:
    def __init__(self, status_code, payload):
        self.status_code = status_code
        self._payload = payload
        self.text = str(payload)

    def json(self):
        return self._payload


class TestFetchGeneCoordinates:
    def test_parses_resolved_genes(self):
        fake_payload = {
            "EPCAM": {"seq_region_name": "2", "start": 47345158, "end": 47387601, "strand": 1},
            "BRCA1": {"seq_region_name": "17", "start": 43044292, "end": 43170245, "strand": -1},
        }
        with patch("requests.post", return_value=_FakeResponse(200, fake_payload)):
            result = fetch_gene_coordinates(["EPCAM", "BRCA1"])
        assert len(result) == 2
        assert set(result["gene"]) == {"EPCAM", "BRCA1"}
        assert result.set_index("gene").loc["EPCAM", "chromosome"] == "2"

    def test_unresolved_gene_is_dropped_not_fabricated(self):
        fake_payload = {"EPCAM": {"seq_region_name": "2", "start": 1, "end": 2, "strand": 1}}
        with patch("requests.post", return_value=_FakeResponse(200, fake_payload)):
            result = fetch_gene_coordinates(["EPCAM", "NOT_A_REAL_GENE"])
        assert len(result) == 1
        assert "NOT_A_REAL_GENE" not in result["gene"].to_list()

    def test_batches_requests(self):
        symbols = [f"GENE{i}" for i in range(5)]
        fake_payload = {
            s: {"seq_region_name": "1", "start": 1, "end": 2, "strand": 1} for s in symbols
        }
        with patch("requests.post", return_value=_FakeResponse(200, fake_payload)) as mock_post:
            fetch_gene_coordinates(symbols, batch_size=2)
        assert mock_post.call_count == 3  # ceil(5/2)

    def test_raises_on_non_200_response(self):
        with patch("requests.post", return_value=_FakeResponse(500, {})):
            with pytest.raises(PipelineError, match="status 500"):
                fetch_gene_coordinates(["EPCAM"])
