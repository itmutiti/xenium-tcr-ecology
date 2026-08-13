"""Statistical closure, figures, manuscript and reproducible release,
the final stage (scripts/17_statistical_closure_and_release/).

`17_statistical_closure_and_release/00_freeze_primary_results.py`-17.09
(freeze, gatekeeping, figures, tables, manuscript, public data release,
software package, reproduction test, calibration regression) are all
implemented. 17.10 (archive/DOI) is not automated: DOI archival is a
one-time, irreversible action, and `10_archive_code_and_create_doi.sh`
prints the manual procedure and exits rather than performing it (see its
own header comment). Two reproducibility issues were found and fixed during
17.08/17.09: a stale environment lock missing r-arrow/r-testthat, and
non-deterministic RNG seeding in xenium_tcr_ecology.graphs.null_models
(now xenium_tcr_ecology.graphs.null_models.deterministic_seed).
"""
