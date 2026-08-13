"""Unit tests for xenium_tcr_ecology.validation.companion_reference_acquisition
(`06_cell_type_annotation/08_acquire_companion_scrna_and_vdj_reference.py`)."""

from __future__ import annotations

import errno
import tarfile
from unittest.mock import patch

from xenium_tcr_ecology.validation.companion_reference_acquisition import (
    VDJ_MEMBERS,
    VDJ_SAMPLES,
    _gsm_series_group,
    ensure_gse287301_vdj_acquired,
)


class TestGsmSeriesGroup:
    def test_real_last_three_digits_become_nnn(self):
        assert _gsm_series_group("GSM8743474") == "GSM8743nnn"


class TestVdjSamples:
    def test_real_exactly_sixteen_pools_registered(self):
        assert len(VDJ_SAMPLES) == 16

    def test_real_pool_names_are_unique(self):
        assert len(set(VDJ_SAMPLES.values())) == len(VDJ_SAMPLES)

    def test_real_gsm_ids_are_unique_and_well_formed(self):
        assert len(set(VDJ_SAMPLES)) == len(VDJ_SAMPLES)
        assert all(gsm.startswith("GSM") for gsm in VDJ_SAMPLES)

    def test_real_known_mapping_matches_geo_record(self):
        # Confirmed against GEO's own per-sample record
        # (Sample_supplementary_file_1), not inferred from numbering --
        # spot-check a few, including the one non-sequential pool name
        # (chip2pool16, not chip2pool8).
        assert VDJ_SAMPLES["GSM8743474"] == "chip1pool1"
        assert VDJ_SAMPLES["GSM8743481"] == "chip1pool8"
        assert VDJ_SAMPLES["GSM8743482"] == "chip2pool1"
        assert VDJ_SAMPLES["GSM8743489"] == "chip2pool16"


class TestEnsureGse287301VdjAcquiredCrossDeviceMove:
    """Regression test for a real defect: a Docker clean-room run (never
    native, where both paths are typically the same filesystem) hit
    `OSError: [Errno 18] Invalid cross-device link` because the original
    code used `Path.replace()` (os.rename(), same-filesystem only) to move
    an extracted file from a system tempfile.TemporaryDirectory() into
    data/external/GSE287301/vdj/ -- a different filesystem from Docker's
    internal /tmp under the /workspace bind mount. Fixed with shutil.move,
    which falls back to copy+remove on exactly this OSError."""

    def _fake_download_file(self, url, dest, **kwargs):
        # Stands in for a real GEO download: writes a real, minimal tar.gz
        # with the two members the real function extracts, so the rest of
        # the pipeline (tarfile.open/extract, then the move under test)
        # runs unmodified against real files, not a mock.
        pool_name = dest.name.split("_", 1)[1].removesuffix(".tar.gz")
        src_dir = dest.parent / "_src" / pool_name
        src_dir.mkdir(parents=True, exist_ok=True)
        for member in VDJ_MEMBERS:
            (src_dir / member).write_text(f"dummy {member} for {pool_name}\n")
        dest.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(dest, "w:gz") as tf:
            tf.add(src_dir, arcname=pool_name)

    def test_move_survives_cross_device_link_error(self, tmp_path):
        def always_raises_exdev(src, dst, *a, **kw):
            raise OSError(errno.EXDEV, "Invalid cross-device link")

        # Patch *both* os.rename (what shutil.move tries first, then
        # falls back from) and os.replace (what Path.replace() -- the
        # original, buggy code -- calls, with no fallback at all). If
        # this ever regresses back to Path.replace()/os.replace(), the
        # forced EXDEV has no fallback to catch it and this test fails
        # with a real, unhandled OSError -- it isn't just asserting
        # shutil.move's own behaviour in the abstract.
        with (
            patch(
                "xenium_tcr_ecology.validation.companion_reference_acquisition.download_file",
                side_effect=self._fake_download_file,
            ),
            patch("os.rename", side_effect=always_raises_exdev),
            patch("os.replace", side_effect=always_raises_exdev),
        ):
            vdj_dir = ensure_gse287301_vdj_acquired(tmp_path)

        for pool_name in VDJ_SAMPLES.values():
            for member in VDJ_MEMBERS:
                moved = vdj_dir / pool_name / member
                assert moved.is_file(), f"{moved} was not created despite forced EXDEV"
                assert moved.read_text() == f"dummy {member} for {pool_name}\n"
