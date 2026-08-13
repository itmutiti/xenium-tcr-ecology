# Simulation / null-model calibration tests

This directory holds the calibration suite for Spatial Graph Construction and Calibration
(`scripts/09_spatial_graph_construction_and_calibration/07_generate_synthetic_ground_truth_patterns.py` and
`08_run_calibration_suite_on_synthetic_data.py`): tests that generate
synthetic spatial point patterns with known, controllable effect sizes and
confirm that every spatial null model used elsewhere in the pipeline
(constrained permutation, degree-preserving, graph-preserving) achieves
nominal Type I error and adequate power at the actual n≈10-patient scale.

**No null model in `src/xenium_tcr_ecology/nullmodels/` may be used outside
this test suite until it has a passing calibration test here** - see
the module docstring of `src/xenium_tcr_ecology/nullmodels/__init__.py`.

Empty at scaffold time; populated when Spatial Graph Construction and Calibration is implemented.
