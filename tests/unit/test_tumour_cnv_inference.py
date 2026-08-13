"""Unit tests for xenium_tcr_ecology.tumour.cnv_inference."""

from __future__ import annotations

import anndata as ad
import numpy as np
from scipy import sparse

from xenium_tcr_ecology.tumour.cnv_inference import _to_dense_frame


def _make_adata(sparse_x=False):
    X = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], dtype=np.float32)
    if sparse_x:
        X = sparse.csr_matrix(X)
    adata = ad.AnnData(X=X)
    adata.var_names = ["A", "B", "C"]
    adata.obs_names = ["c1", "c2"]
    adata.layers["lognorm"] = X
    return adata


class TestToDenseFrame:
    def test_dense_input(self):
        adata = _make_adata(sparse_x=False)
        result = _to_dense_frame(adata, ["A", "C"], "lognorm")
        assert list(result.columns) == ["A", "C"]
        assert result.loc["c1", "A"] == 1.0
        assert result.loc["c2", "C"] == 6.0

    def test_sparse_input_is_densified(self):
        adata = _make_adata(sparse_x=True)
        result = _to_dense_frame(adata, ["A", "B"], "lognorm")
        assert result.shape == (2, 2)
        assert result.loc["c1", "B"] == 2.0
