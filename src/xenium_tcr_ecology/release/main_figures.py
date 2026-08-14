"""Assembles the manuscript main figures from already-generated,
already-validated PDF outputs (`17_statistical_closure_and_release/02_
generate_all_main_figures.py`). No analysis is regenerated here: every
source PDF was already produced, tested, and governance-documented by
the phase that owns it. Four of the six figures are retrospectively-built
composites combining multiple already-validated data sources into one
figure (`xenium_tcr_ecology.release.redesigned_main_figures`,
`17_statistical_closure_and_release/11_build_redesigned_manuscript_
figures.py`).

The figure set was redesigned (2026-07-13; see `docs/analysis_
amendments.md`) after the original six were finalised before later
findings existed: a second, independent framework-generalisation
dataset, the Q2 feature-set sensitivity check, the Q3 covariate-ablation
check, and the TCR probe-vs-VDJ validation (which had no figure at all).
The sequence below folds each into the figure it most directly bears on,
and consolidates the two HPV figures into one.

Filenames are descriptive (e.g. `Figure_4_Clone_Ecological_Structure_
Variance_Partition_With_Sensitivity.pdf`), not the internal
`governance/analysis_registry.tsv` `analysis_id` shorthand
(`q2_variance_partition_confirmatory`); `analysis_id` is retained as its
own `MANIFEST.tsv` column for cross-referencing back to the registry.

Every figure is rendered to PDF (primary), a 300 DPI PNG, and an SVG
(`xenium_tcr_ecology.infra.figure_export`) -- no JPEG, since lossy
compression is inappropriate for statistical plots with sharp text and
thin lines.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.figure_export import export_pdf_to_png_svg_and_tiff

MAIN_FIGURE_MANIFEST: list[dict] = [
    {
        "figure_number": 1,
        "analysis_id": "q1_framework_generalisation",
        "descriptive_name": "Calibrated_Null_Models_Control_False_Positive_Associations",
        "source_path": "reports/manuscript_figures/framework_generalisation_three_tumour_types.pdf",
        "caption": "Calibration of spatial null models across synthetic and independent tumour datasets. Empirical rejection rate of three spatial null-model constructions -- constrained permutation, degree-preserving permutation, and graph-preserving (graph-hop-distance) permutation -- across seven injected effect sizes (0 = true null). In all panels, points show the empirical rejection rate at p < 0.05; vertical lines show the 95% Clopper-Pearson confidence interval (CI); the dashed line marks the nominal alpha = 0.05; marker shape and colour identify the null model (shared legend below). (A) Synthetic tissue data: 10 synthetic patients per effect size. Each patient has a spatially contiguous tumour region and a clone population sampled with a known, injected distance-dependent bias toward the tumour. (B, C) The same three null models and test statistic were applied, unmodified, to two independent, single-cell-resolution Xenium datasets not used in framework development: breast cancer (B; Janesick et al. 2023) and colorectal cancer (C; de Oliveira et al. 2025). Each of these datasets comprises one tissue section rather than multiple patients, so each point summarises 10 spatially contiguous resamples of that section rather than 10 independent patients as in (A). A known clone effect of the same design as in (A) was injected onto each section's own spatial topology.",
    },
    {
        "figure_number": 2,
        "analysis_id": "tcr_probe_vdj_ground_truth_validation",
        "descriptive_name": "TCR_Probe_Assignments_vs_Paired_scTCR_seq",
        "source_path": "reports/tcr/probe_vdj_ground_truth_validation.pdf",
        "caption": "Concordance between Xenium CDR3-probe assignments and paired single-cell TCR sequencing. (A) Per-patient count of patient-specific complementarity-determining region 3 (CDR3) probes statistically mapped to an intended patient (n = 105 probes across 11 patients, including the additional ameloblastoma-specimen patient P10; Methods). Bars are stacked by whether the probe's CDR3 amino-acid sequence was also detected in that same patient's independently generated single-cell TCR-sequencing (scTCR-seq) VDJ data (GSE287301): confirmed (blue) or not confirmed (grey). Patients are ordered by total mapped-probe count. (B) For the 80 confirmed probes, Xenium spatial detection rate (fraction of that patient's T cells with the probe detected) is plotted against the same probe's clonal abundance rank within its patient and T-cell receptor chain in the independently sequenced VDJ data (rank 1 = most abundant clonotype); both axes are log-scaled. The annotated statistic (two-sided Spearman rank correlation, rho = -0.19, p = 0.099) indicates the association between Xenium detection rate and independently measured clonal abundance.",
    },
    {
        "figure_number": 3,
        "analysis_id": "q2_discrete_vs_continuous_structure_test",
        "descriptive_name": "Discrete_vs_Continuous_Clone_Structure_Test",
        "source_path": "reports/clone_ecology/structure_test.pdf",
        "caption": "Discrete-versus-continuous clone ecological-structure test. Two model-selection criteria were evaluated across candidate discrete cluster counts (K = 1-6) for two independent descriptor sets. Blue points and lines show each criterion's value at each K; the vermillion point marks the K selected by that criterion. (A, C) Bayesian information criterion (BIC; lower is better) of the best-fitting K-component discrete mixture model at each K. The dashed line marks the BIC of a one-factor continuous-latent-axis model fitted to the same data, for reference; the vermillion point marks the K with the lowest discrete-model BIC. (B, D) Tibshirani gap statistic at each K; error bars show +/-1 SD from the statistic's reference-distribution resampling. The vermillion point marks the K with the highest gap statistic. (A, B) Primary descriptor set (n = 261 clone-sections). (C, D) Sensitivity descriptor set (17 features, n = 120 clone-sections), meeting a stricter completeness threshold.",
    },
    {
        "figure_number": 4,
        "analysis_id": "q2_variance_partition_confirmatory",
        "descriptive_name": "Variance_Partition_Feature_Set_Sensitivity",
        "source_path": "reports/manuscript_figures/variance_partition_with_sensitivity.pdf",
        "caption": "Variance-component partition of the clone ecological-structure score. Three-component variance partition (spatial context, patient, clonal identity) of the clone ecological-structure score, estimated by a three-level nested random-intercept model (clone-sections nested within clones, nested within patients). The model was fit to 261 clone-sections, 158 clones, and 10 patients, using the prespecified 11-feature descriptor set (blue) and a sensitivity feature set excluding cycling_fraction (vermillion, 10 features). Points show the estimated proportion of total variance attributed to each component; horizontal lines show the 95% bootstrap CI (500 replicates, lme4::bootMer). Grey segments connect each component's point estimate between the two feature sets. Components are ordered top-to-bottom by the primary (11-feature) estimate.",
    },
    {
        "figure_number": 5,
        "analysis_id": "q3_barrier_topology_confirmatory",
        "descriptive_name": "Suppressive_Myeloid_Barrier_Association_With_Engagement",
        "source_path": "reports/manuscript_figures/barrier_topology_with_ablation.pdf",
        "caption": "Barrier-topology model of clone-tumour engagement. Fibroblast and suppressive-myeloid barrier fraction is the composition of cell types encountered along the shortest spatial-graph path between a clone and its nearest tumour cell. (A) Fixed-effect estimates for the two barrier covariates from a nested random-intercept mixed-effects model of clone-tumour engagement ratio, adjusted for cell-intrinsic state composition and niche composition. The model was fit to 152 clone-sections (105 clones, 10 patients); clone-sections with no intermediate barrier cells were excluded by definition (Methods). Points show the fixed-effect estimate; horizontal lines show the 95% bootstrap CI (500 replicates). The dashed line marks zero effect. (B) Covariate-ablation stress test. The suppressive-myeloid barrier estimate was re-fitted using the barrier terms alone, each adjustment block alone (state composition; niche composition), the single largest individually added covariate (a niche-archetype fraction), and the full joint-adjustment model (state and niche composition together, vermillion). Points show the fixed-effect estimate; horizontal lines show the 95% Wald CI (estimate +/- 1.96 x SE). This differs from the bootstrap CI used in (A). (C) Raw (unadjusted) Pearson correlation between each barrier covariate and engagement in this cohort (circles), shown alongside the one available published cross-cancer raw correlation for the suppressive-myeloid analogue (diamond: stromal density versus T-cell exclusion, lung cancer, Grout et al. 2022). The dashed line marks zero correlation. The two studies differ in cohort, spatial unit, and analytical scale, and no formal test of agreement between them is shown (Methods).",
    },
    {
        "figure_number": 6,
        "analysis_id": "hpv_primary_contrast_consolidated",
        "descriptive_name": "HPV_Underpowered_Composition_Structure_and_Discordance",
        "source_path": "reports/manuscript_figures/hpv_composition_structure_and_discordance_qc.pdf",
        "caption": "HPV-status comparison and clinical-molecular concordance. (A, B) Significance (-log10 Benjamini-Hochberg [BH] q-value) of patient-level comparisons between HPV-positive and HPV-negative patients (n = 4 versus n = 4, exact two-sided Wilcoxon tests). (A) Cellular-composition lineages (n = 12 lineages). (B) Ecosystem/clone-structure metrics (n = 13: abundance and mixing index for each of six niche archetypes, plus the clone ecological-structure score). Points mark each comparison's -log10(q); the dashed line marks q = 0.05. Comparisons are ordered by significance within each panel. (C) Per-patient Xenium HPV16 E6/E7 probe-positive cell fraction (log scale) for the six patients with probe coverage. Points are coloured by concordance with clinical p16 immunohistochemistry status: confirmed positive (probe-positive and clinically p16-positive), clinically positive/probe-negative, or probe-positive/clinically untested.",
    },
]


def build_main_figures(project_root: Path) -> dict:
    figures_dir = project_root / "figures" / "main"

    missing = [
        entry["source_path"]
        for entry in MAIN_FIGURE_MANIFEST
        if not (project_root / entry["source_path"]).is_file()
    ]
    if missing:
        raise PipelineError(f"Missing source PDF(s), cannot assemble main figures: {missing}")

    figures_dir.mkdir(parents=True, exist_ok=True)
    current_stems = {
        f"Figure_{entry['figure_number']}_{entry['descriptive_name']}"
        for entry in MAIN_FIGURE_MANIFEST
    }
    for stale in figures_dir.glob("Figure_*"):
        if stale.stem not in current_stems:
            stale.unlink()

    manifest_rows = []
    for entry in MAIN_FIGURE_MANIFEST:
        source = project_root / entry["source_path"]
        stem = f"Figure_{entry['figure_number']}_{entry['descriptive_name']}"
        dest = figures_dir / f"{stem}.pdf"
        dest.write_bytes(source.read_bytes())
        export_pdf_to_png_svg_and_tiff(dest, figures_dir, stem)
        manifest_rows.append(
            {
                "figure_number": entry["figure_number"],
                "filename": dest.name,
                "analysis_id": entry["analysis_id"],
                "source_path": entry["source_path"],
                "caption": entry["caption"],
            }
        )

    manifest_df = pd.DataFrame(manifest_rows)
    manifest_path = figures_dir / "MANIFEST.tsv"
    manifest_df.to_csv(manifest_path, sep="\t", index=False)

    return {
        "n_figures": len(manifest_rows),
        "figures_dir": str(figures_dir),
        "manifest_path": str(manifest_path),
    }
