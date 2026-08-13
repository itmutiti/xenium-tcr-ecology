# Changelog

All notable changes to this project are documented here. Format loosely
follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased]

### Added
- `tools/build_markdown_pdf.py`: forces line breaks into long unbroken
  identifiers/paths in rendered tables and body text - `xhtml2pdf`
  does not enforce a fixed table-cell width or honour `word-break`/
  `<wbr>`/an explicit `<colgroup>` width against a long token with no
  whitespace, confirmed by rendering several candidate fixes before
  finding one that actually works with this renderer; previously such
  cells overflowed into, or clipped past, neighbouring columns.
- PNG (300 DPI) and SVG export alongside every existing figure PDF
  (`xenium_tcr_ecology.infra.figure_export`, via `pdftocairo`).
- A general-purpose Markdown→PDF build tool
  (`tools/build_markdown_pdf.py`); `tools/build_execution_manual_pdf.py`
  now delegates to it.
- `tools/run_pipeline.sh`: a reviewer-facing launcher that probes Docker,
  then Apptainer, then native conda/Snakemake, in that priority order,
  and runs the complete workflow through the first usable
  route - see `docs/execution_manual/EXECUTION_MANUAL.md`
  "Reproducibility" and "Running the pipeline". Regression coverage in
  `tests/unit/test_run_pipeline_launcher.py`, including a direct
  cross-route DAG-identity assertion.
- `containers/docker-entrypoint.sh`: self-registers an arbitrary running
  UID/GID (with a writable `$HOME`) so the Docker image can be run
  non-root, keeping files it writes into the bind-mounted repository
  owned by the invoking user rather than root.
- "Computational requirements" section in `docs/execution_manual/
  EXECUTION_MANUAL.md` (CPU/RAM/disk/network/runtime, each labelled
  measured or recommended) and a concise summary in `README.md`.
- Implemented `17_statistical_closure_and_release/
  06_build_public_data_release.py` (`xenium_tcr_ecology.release.
  public_data_release`), previously a permanent `NotImplementedError`
  stub. Builds `release/data/` from three sources - the frozen primary
  confirmatory results, the 10 manuscript result tables, and the data
  dictionary - under a CC BY 4.0 licence, with a hash-manifested `MANIFEST.json`/
  `checksums.sha256` that refuses to silently rebuild if any source
  input changed since the last build. Deliberately excludes anything
  cell-level or per-patient beyond the de-identified codes already
  public in the manuscript's own tables. Does not itself substitute for
  the manual privacy/licensing review the release procedure still
  requires - see `docs/execution_manual/EXECUTION_MANUAL.md` "Release
  actions". Regression coverage in `tests/unit/
  test_release_public_data_release.py`.

### Changed
- `docs/execution_manual/EXECUTION_MANUAL.md`: restructured (purpose,
  quick start, and a concise verification-status summary now open the
  document; detailed per-route explanations follow) and audited against
  actual repository state rather than trusted as previously written.
  Resolved two real factual contradictions: `environment.lock` was
  described in one place as committed/tracked and in another as
  generated-locally-only (it is committed; the "Building locally"
  Docker instructions no longer imply a lock-generation step is needed
  first). Disk-usage figures (52 GiB primary cohort vs. ~110 GB vs. ~60
  GB, previously easy to mistake for disagreeing) are now stated
  separately and re-measured (primary cohort ~53 GiB, six external
  datasets ~7 GB, all seven combined ~60 GB, datasets plus every
  generated output ~100 GB). Updated Docker/Apptainer verified status:
  both have now independently reached `computational_analysis_complete`
  on a real clean-room instance (2026-08-07, part of the
  cross-environment reproducibility investigation below), superseding
  the document's previous "neither has completed a full run" and
  "Docker was not exercised by any clean-room run" statements. Also
  corrected a Python test-suite nuance carried over inaccurately from an
  earlier pass: on a public checkout missing the two optional private
  governance files, one dependent test fails (no skip guard, by design)
  and a different one skips (explicit guard) - previously both were
  described as skipping. Re-verified the
  documented DAG job count (152, confirmed still current against
  `tests/unit/test_workflow_dag.py::TestDocumentedJobCountMatchesLiveDag`,
  which compares it to a live `snakemake --forceall -n` dry run) and
  updated the Python test count to its current measured value (727, not
  709).
- `17_statistical_closure_and_release/10_archive_code_and_create_doi.sh`:
  rewritten to document both DOIs the release architecture actually
  requires (previously covered only the software DOI). The software DOI
  fires automatically from the GitHub-Zenodo webhook on a published
  GitHub Release; the data DOI is a separate, manual Zenodo upload of
  `release/data/`, which is gitignored and never swept into the tagged
  software archive regardless of the release-packaging option used.
  Still deliberately non-automated and exits 1 on every run - DOI
  archival is a one-time, irreversible, human-gated action, not a gap
  awaiting further implementation.
- Set `governance/analysis_registry.tsv`'s `registered_by` column,
  `config/geo/sample_manifest_input.yaml`'s `compiled_by` field, and
  `06_create_analysis_registry.py`'s `--registered-by` CLI default to
  "Irvine Tatenda Mutiti" ahead of the first public release.
- Reordered the documented execution-route priority to Docker
  (preferred), Apptainer (fallback), native conda/Snakemake (final
  fallback, still fully supported) - previously presented as three
  equally-optional routes with native framed as canonical. See
  `docs/execution_manual/EXECUTION_MANUAL.md` "Reproducibility".
- `environment/conda/environment.lock` is now committed (previously
  gitignored as "machine-generated"). It's a genuine dependency lock,
  platform-pinned to linux-64 -- the one platform `containers/Dockerfile`
  and `containers/Apptainer.def` ever build on, so committing it carries
  none of the cross-platform-drift risk the native route's own fresh
  `main.yml` resolve has on an arbitrary reviewer host. This closes a
  real architectural gap: backend selection previously skipped Docker
  and Apptainer outright on a fresh checkout because the lock they build
  from didn't exist yet, forcing a manual native bootstrap first just to
  generate it. `environment/conda/environment.resolved.yml` and
  `environment/R_sessionInfo.txt` (machine-specific: hostname,
  kernel, capture timestamp) remain gitignored, generated provenance.
- `tools/run_pipeline.sh`: added `--backend auto|docker|apptainer|native`
  (env var `XENIUM_FORCE_ROUTE` remains a legacy synonym). An explicit
  `--backend` never silently falls back -- it explains why the requested
  backend is unusable, states that no automatic fallback occurred, and
  lists available alternatives with the exact command for
  each. Backend probing no longer references
  `environment/conda/environment.lock`'s presence at all -- Docker's
  usability is a pure capability probe, unconditional on any generated
  file. Every rejected backend's reason, the selected backend, the image
  tag/digest (or `.sif` path/checksum), and the lock checksum are now
  saved to `results/logs/run_pipeline/provenance.json` in addition to
  the existing `route_selection.log`. New regression coverage in
  `tests/unit/test_run_pipeline_fresh_checkout.py`, exercising the real
  script against fake backend stubs on a fresh checkout (no
  lock, no image, no `.sif`).
- `containers/Dockerfile`: builds from the committed
  `environment/conda/environment.lock` when present (the normal case)
  and falls back to a fresh `mamba env create -f main.yml` resolve
  automatically if it's ever absent, rather than requiring the
  lock to exist beforehand.
- Renamed the 6 main figures from internal `analysis_id` shorthand to
  descriptive, reviewer-facing names; `analysis_id` retained as a
  `MANIFEST.tsv` column for traceability.
- Centralised random-seed configuration (`xenium_tcr_ecology.infra.
  seeding`, `config/config.yaml`'s `default_seed`/`annotation_seed`, and
  an equivalent R helper) across scripts that previously duplicated the
  same hardcoded seed locally. Verified value-preserving (byte-identical
  output on re-run) across the full Python and R test suites.
- Corrected the release-gatekeeping statistic for the barrier-topology
  confirmatory claim (`q3_barrier_topology_confirmatory`) from an ad hoc
  Wald test to a properly nested single-coefficient likelihood-ratio
  test, matching the claim actually being gated. Conclusion (PASS)
  unchanged. See `docs/analysis_amendments.md`.
- Renamed `src/xenium_tcr_ecology/checkpoint/` →
  `src/xenium_tcr_ecology/external_checkpoint/` (avoids a naming
  collision with "immune checkpoint" biology terminology used elsewhere
  in the same codebase).
- Separated the two external-facing publication actions (public data
  export, code archival/DOI creation) from the default Snakemake target;
  `rule all` now stops at a `computational_analysis_complete` checkpoint,
  and `publication_release` is a separate, explicitly-named target run by
  hand.
- Renamed all 17 computational workflow stage folders (`scripts/`,
  `workflow/rules/`) and their governance decision logs from
  abbreviated phase-number names to descriptive workflow names
  throughout the repository (code, tests, docs, governance prose).
- Populated `manifests/dataset_registry.yaml` with all datasets actually
  used by the pipeline.

### Fixed
- `ensure_gse287301_vdj_acquired()` (`06_cell_type_annotation/08_acquire_
  companion_scrna_and_vdj_reference.py`) used `Path.replace()` to move a
  downloaded/extracted file from a `tempfile.TemporaryDirectory()` into
  `data/external/GSE287301/vdj/` -- `OSError: [Errno 18] Invalid
  cross-device link` under Docker specifically, where the container's
  internal `/tmp` and the `/workspace` bind mount are different
  filesystems (never hit natively, where both are typically the same
  filesystem). Found by this project's first real, full-dataset Docker
  clean-room execution. Fixed with `shutil.move`, which falls back to
  copy+remove on exactly this error; regression test in
  `tests/unit/test_validation_companion_reference_acquisition.py` forces
  the failure via a mocked `os.rename`/`os.replace` and asserts the move
  still completes.
- `find_project_root()` depended on `commandArgs()`'s `--file=`, which is
  absent when an R script is merely `source()`-d rather than run
  directly; replaced with a self-contained `getwd()`-based lookup.
- A title-clipping rendering bug affecting 3 figures.
- A stale citation in the manuscript's `.qmd` source (a glob pattern
  referencing a governance file path that no longer existed after the
  workflow-stage rename above).
- Documented, without retrofitting, a reproducibility gap: some
  stochastic analyses were not RNG-seeded prior to the seed
  centralisation above, so exact CI bounds were not byte-reproducible
  across re-runs before that fix (qualitative conclusions were stable
  throughout).
- `data/standardised/`'s per-section symlinks used absolute host paths
  (`link_path.symlink_to(source.resolve())`), silently broken under any
  bind mount at a different path (Docker's `/workspace`, and equally
  Apptainer's whole-pipeline mode) even though the target file was
  present - found because two tests failed only under Docker,
  never natively, on the identical mounted filesystem. Now relative,
  correct regardless of mount point.
- `.dockerignore` had no entry for `manuscript/`, `docs/Journals/`, or
  the private prompt file - Docker doesn't consult `.gitignore` on its
  own, so a build would have copied these into the image's build
  context (confirmed: context dropped from would-be 300+ MB to 10 MB
  after the fix).
- Docker's `ENTRYPOINT` ran as root, leaving root-owned files in the
  bind-mounted repository that then broke a subsequent native run with
  a `PermissionError` - see `containers/docker-entrypoint.sh` above.
- `tests/unit/test_workflow_dag.py` hardcoded `.venv/bin/snakemake`,
  whose shebang bakes in this machine's absolute path and breaks under
  any bind mount at a different path; now resolves via `PATH` first,
  falling back to the native `.venv`.
- `containers/Apptainer.def`'s `%post`/`%test` assumed a bash shell
  (`pipefail`, `source`), but Apptainer runs section scriptlets under
  `/bin/sh` by default; `%files` copied to `/tmp`, which Apptainer's own
  default configuration bind-mounts from the host during a build,
  silently shadowing the copied file. Neither was ever caught before,
  because the only prior Apptainer evidence was the `docker-daemon://`
  conversion route, which exercises neither `%post` nor `%test`. Found
  and fixed by this project's first genuine `--fakeroot` build; see
  `docs/execution_manual/EXECUTION_MANUAL.md` "Reproducibility" →
  Apptainer for the full account.
- `containers/Dockerfile` and `containers/Apptainer.def` built `FROM`/
  `From:` an unpinned `condaforge/miniforge3:latest` tag. A Docker image
  and an Apptainer image built from the identical committed
  `environment.lock`, minutes apart on the same clean-room VM, produced a
  different result for one bootstrap-based confirmatory model
  (`suppressive_myeloid_barrier_fraction`: -0.368259 in Docker vs.
  -0.343317 -- the manuscript value -- in Apptainer and native). Both
  files now pin the identical immutable base-image digest
  (`sha256:8edc74dee5b29511a9c1a8941c77efbf438b47edffb6320ab8675dc06bd618c2`);
  rebuilding Docker from scratch against the pinned digest eliminated the
  discrepancy exactly (confirmed both in an isolated single-script
  re-test and in a subsequent full 152-job clean-room pipeline run). Full
  experimental account, including every ruled-out alternative hypothesis,
  in `docs/reproducibility_investigation.md`.
- `xenium_tcr_ecology.tumour.gene_coordinates.fetch_gene_coordinates()`
  called Ensembl's REST API with no retry logic (a bare `requests.post`,
  raising on any non-200 response or exception), unlike this project's
  other external-dependency helper (`infra.download.download_file`).
  Surfaced by a genuine cache-miss during a fresh clean-room run (the
  panel's gene set was not a subset of the previously cached lookup
  table), which triggered a real API call that hit two transient
  Ensembl-side 503/500 responses under batch load -- confirmed not a
  local connectivity problem (`/info/ping` returned 200 throughout from
  the host, the VM, and inside the container). Added a bounded
  retry-with-backoff (`_post_batch_with_retries`, 5 attempts, 5s delay).

## [1.0.0] - 2026-07-11

### Added
- Full implementation of all 17 computational workflow stages (literature/
  novelty planning through statistical closure and manuscript rendering),
  each with test-first implementation (Python `pytest` + R `testthat`)
  and a governance decision log.
- External validation against five independent public datasets beyond the
  primary GSE300147 cohort: GSE287301 (companion scRNA-seq), GSE103322
  (Puram et al. 2017), GSE139324 (Cillo et al. 2020), TCGA-HNSC bulk
  RNA-seq, and an independent Xenium breast-cancer dataset (Janesick et
  al. 2023) used to confirm the calibrated null-model framework
  generalises beyond the primary cohort.
- Three prespecified confirmatory analyses (Q1 framework generalisation,
  Q2 clone ecological-structure variance partition, Q3 barrier-topology
  model) plus a capped, prospectively-power-analysed HPV-stratified
  analysis, all frozen and gatekept in the statistical-closure stage.
- Installable, versioned software package (`xenium_tcr_ecology` v1.0.0),
  Docker container image, Snakemake workflow orchestration, and a
  reproduction/calibration-regression test suite.
- `tools/scaffold_repository.py`, kept as a permanent idempotent scaffold
  and policy-verification tool.
