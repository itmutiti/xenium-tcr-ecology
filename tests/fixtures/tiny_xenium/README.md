# Tiny synthetic/subsampled Xenium fixture - scope specification

**Status: not yet created.** This is a concrete, de-risked specification
for a planned CI fixture, so the decision to defer construction is
documented and actionable rather than an open-ended placeholder.
`snakemake -n` (a dry run) validates that the DAG's declared dependencies
are structurally well-formed, but cannot detect a rule that races ahead
of an upstream dependency it never declared - that class of defect is
only observable by actually executing the DAG under parallel
scheduling.

## Purpose (and explicit non-purpose)

This fixture exists to let CI **execute** a bounded, representative slice
of the Snakemake DAG end-to-end, in minutes, to catch **wiring
breakage**: missing/incorrect `input:` dependencies, file-discovery and
schema-handling bugs, and broken Python/R invocation.

**It must never be described as, or used as evidence for, scientific
validity.** It does not and cannot test whether the pipeline's
statistical conclusions are correct, whether thresholds are well-chosen,
or whether the full 28-patient cohort's results are reproducible - those
are covered by `tests/simulation/` (null-model calibration on synthetic
data with known ground truth) and by an actual full-cohort run
respectively. A tiny, heavily subsampled or synthetic dataset has no
statistical power and no biological realism to validate against.

## Proposed scope: phases 1–4 only (project setup through core QC)

Deliberately bounded to the phases where exercising the DAG with
representative rules is achievable **without inventing new scientific
assumptions**:

- **Phase 1** (project setup and governance) - no data dependency, already
  fully exercised by the existing `pytest`/R suites and by direct script
  invocation.
- **Phase 2** (raw data ingestion) - needs only a tiny tar archive with a
  plausible internal structure; a download-then-verify pair of rules like
  `phase02_01`/`phase02_02` is exactly the shape of dependency this
  fixture needs to exercise under real execution.
- **Phase 3** (SpatialData import) - needs one non-synthetic tiny
  Xenium-format section; see "Data source" below.
- **Phase 4, rules 4.00–4.03 and 4.06–4.09 only** (core QC - the same 9
  rules `qc_and_integrity_release` now gates on, see
  `workflow/rules/_checkpoints.smk`). Exercises Python (4.00–4.03) *and*
  R (4.06, 4.08) in the same run. **Deliberately excludes 4.04/4.05**
  (need Cell Type Annotation output -- exactly the
  kind of phase06+ biological dependency this fixture should not need)
  and everything from phase05 onward (normalisation, annotation, TCR
  calling, and beyond all require either biological signal to be
  meaningful or external comparison data - neither of which a
  tiny fixture can provide).

This scope reaches the `canonical_objects` and (a narrowed)
`qc_and_integrity_release` checkpoints - two of the fifteen checkpoints
in the full chain - which is sufficient to prove the DAG executes
correctly in dependency order, without extending into phases whose
correctness depends on real biology.

## Data source: a licensed, already-present public dataset - not fabricated synthetic biology

**Do not fabricate synthetic Xenium-format files from scratch.** The
6 mandatory per-section files (`MANDATORY_ROLES` in
`src/xenium_tcr_ecology/ingest/xenium_inventory.py`:
`cell_boundaries`, `cell_feature_matrix`, `cells`, `morphology`,
`nucleus_boundaries`, `transcripts`) constitute a semi-standardised
10x Genomics multi-file format (Parquet, HDF5, and - the hard part - a
pyramidal OME-TIFF for `morphology`). Recreating this format correctly
from nothing risks subtly-wrong files that pass a naive check but don't
exercise what `spatialdata-io`'s reader actually does.

Instead: **spatially crop a small bounding-box region from one of the
already-present, already-CC-BY-4.0-licensed public Xenium datasets** this
repository already uses for external validation (`16_external_validation_and_generalisation`)
- `data/external/spatial/Xenium_Janesick_BreastCancer_Rep1/` or
`Xenium_Oliveira_ColorectalCancer_P1/`. This produces
correctly-formatted Xenium output (transcripts/cells/boundaries filtered
by a spatial bounding box; `cell_feature_matrix.h5` subset to the
retained cell barcodes) with no fabricated biology - only existing
tissue data, downsampled by extent. No licensing concern:
redistributing a small crop of an already-CC-BY-4.0 dataset is
unambiguously covered by the same licence as the whole.

**The hard part**: `morphology.ome.tif.gz` is a pyramidal
(multi-resolution) OME-TIFF; cropping it to match the same bounding box
without corrupting the pyramid structure `spatialdata-io` expects needs
care (region-of-interest extraction per pyramid level, not just a naive
crop of the full-resolution layer) and its own dedicated validation
(load the cropped result back through `spatialdata-io` and confirm it
aligns with the cropped transcripts/cells before trusting it in CI). This
is the reason construction remains deferred - it warrants focused effort
and independent validation of its own, not a small addition alongside
unrelated work.

## What "done" looks like

- `tests/fixtures/tiny_xenium/<GSM-like-id>/` containing the 6 mandatory
  files for 1–2 cropped sections, target total size a few MB (not the
  138 MB+ of a full section).
- A CI job (extending `.github/workflows/ci.yml`'s `smoke-tier`, which
  currently only dry-runs) that stages this fixture where
  `03_spatialdata_import` expects `data/standardised/<section_id>/`,
  runs `snakemake --cores 2 qc_and_integrity_release` (not `-n`),
  and asserts the run succeeds - this is the actual DAG-execution check
  the dry-run alone cannot provide.
- A short comment at the top of the CI job's own YAML stating plainly
  what it does and does not validate, mirroring this file's "Purpose"
  section, so a CI failure (or a green run) is never misread as a
  scientific-validity signal in either direction.
- A unit or integration test asserting the fixture's own files pass
  `spatialdata-io`'s reader in isolation, so a future accidental edit to
  the fixture fails fast and locally, not only inside a full CI DAG run.

## Why this remains deferred, not built yet

This fixture is worth building only if it can be done without fabricating
new scientific assumptions, and the scope above achieves that - but the
OME-TIFF pyramid cropping is non-trivial, independently error-prone
engineering work that deserves focused effort and iterative validation of
its own. Until it exists, exercising the DAG's execution order
(rather than only its dry-run structure) requires either running the
pipeline against real data directly, or an equivalent full clean-room
execution.
