# Architecture

## Package vs. pipeline separation

Analysis logic lives in `src/xenium_tcr_ecology/` rather than directly in
`scripts/<phase>/`, so it is importable, unit-testable, and reusable across
stages. `scripts/<phase>/` stays a thin CLI wrapper per stage: each script
parses arguments, resolves the project root, calls one or more functions
from the package, and reports a one-line summary. This split means the
package can be installed and used independently of the Snakemake pipeline
(`pip install -e .`), and every non-trivial function has a corresponding
unit test in `tests/unit/` or `tests/r_unit/` without needing to run the
full pipeline.

## Repository layout

```
src/xenium_tcr_ecology/   installable package: reusable, tested analysis logic
scripts/<NN_stage>/       thin CLI entry points, one subfolder per computational workflow stage
workflow/rules/           Snakemake rule definitions, one file per stage
governance/               decision logs, analysis registry, taxonomy version log,
                          freeze decisions
config/                   thresholds, seeds, graph/interaction parameters, policy
manifests/                dataset registry, project-root marker, phase registry
environment/, containers/ conda spec + lock, Docker + Apptainer image definitions
                          (see tools/run_pipeline.sh: Docker preferred, Apptainer
                          fallback, native conda/Snakemake final fallback)
tests/                    unit/, simulation/ (null-model calibration), r_unit/, fixtures/
data/, results/, reports/, figures/, tables/   gitignored except .gitkeep
docs/                     architecture + phase reference (mkdocs)
release/                  installable software package, public data bundle (once built)
```

## Governance-as-code

Two scientific decisions are gated mechanically rather than by convention
alone, so a reviewer does not have to trust that a rule was followed by
hand:

- **Taxonomy freeze/revise gate**: the T-cell-state taxonomy produced by
  `11_clone_spatial_descriptors` is provisional until External Checkpoint
  Validation's module-transfer and directional-consistency checks pass a
  predeclared rule (`12_external_checkpoint_validation/03_decide_freeze_or_revise.py`,
  recorded in `governance/freeze_decision.tsv`). No script downstream of
  that gate may use the taxonomy until the decision is recorded.
- **HPV-contrast cap**: at most one primary HPV-positive-vs-negative
  contrast may be registered for confirmatory use in the entire project
  (`15_hpv_stratified_analysis/01_prespecify_primary_hpv_contrasts.py`,
  recorded in `governance/hpv_primary_contrasts.yaml`); this is enforced
  by `tools/scaffold_repository.py`'s policy check
  (`config/reproducibility_policy.yaml`), not left as an unenforced
  convention.

## Prespecification

`governance/analysis_registry.tsv` records every hypothesis-bearing or
formal validation analysis before it is run: its hypothesis, unit of
analysis, exclusion criteria, primary endpoint, and multiplicity family.
`docs/analysis_amendments.md` records every analysis added after the
initial registry was frozen, with its rationale and relationship to the
prespecified set.

See `EXECUTION_MANUAL.md` for how to run the pipeline, and `phases.md` for
the full per-stage script reference.
