"""Unit tests for xenium_tcr_ecology.release.reproduction_test (`17_statistical_closure_and_release/08_run_end_to_end_reproduction_test.sh`)."""

from __future__ import annotations

from xenium_tcr_ecology.release.reproduction_test import classify_run_result


class TestClassifyRunResult:
    def test_real_zero_exit_code_is_pass(self):
        assert classify_run_result(returncode=0) == "PASS"

    def test_real_nonzero_exit_code_is_fail(self):
        assert classify_run_result(returncode=1) == "FAIL"

    def test_real_negative_exit_code_from_a_real_signal_is_fail(self):
        assert classify_run_result(returncode=-9) == "FAIL"
