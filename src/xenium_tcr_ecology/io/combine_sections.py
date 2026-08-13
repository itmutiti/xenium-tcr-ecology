"""Concatenate all per-section AnnData objects into one combined analysis
object (`03_spatialdata_import/05_build_combined_analysis_object.py`), with globally unique cell identifiers (prefixed by
section_id, since raw Xenium cell_id values like 'aaadggoi-1' are only
unique within a single section and collide across sections) and an
explicit patient/run hierarchy already present from `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py`'s join.

Gene panels are not identical across sections: sections differ by a ~70-feature block whose names are
literal CDR3 amino-acid sequences (e.g. '230913_CASRPLSYNEQFF_TRB'),
confirming these are patient-specific TCR probes as described in the
specification. Panels vary by processing batch, not strictly per-patient: P09
and P13's sections carry an *identical* non-core probe set, meaning a
custom-probe batch's panel can include multiple patients' clones together
(exactly the "off-patient probe" structure TCR Clonal Analysis's false-positive-rate
estimation needs). Concatenation therefore uses an outer join (union of all
genes across all sections) rather than requiring an identical panel.

Because AnnData's outer-join concat fills a missing gene with 0 for any
section that never measured it, and an actual 0-count and "not in this
section's panel at all" are not the same thing, this module also writes a
per-section panel-membership record (results/tables/03_spatialdata_import/
gene_panel_membership.parquet) so TCR Clonal Analysis (which needs this distinction to
correctly identify a probe as off-patient/negative-control rather than
"detected zero") is not misled by the 0-fill.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

REQUIRED_OBS_COLUMNS = ["section_id", "patient_id", "cell_id"]


def build_combined_object(
    anndata_root: Path, output_path: Path, panel_membership_path: Path
) -> dict:
    h5ad_paths = sorted(anndata_root.glob("*.h5ad"))
    if not h5ad_paths:
        raise PipelineError(
            f"No .h5ad files found under '{anndata_root}'. Run `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` first."
        )

    adatas = {}
    panel_by_section: dict[str, set[str]] = {}
    for path in h5ad_paths:
        adata = ad.read_h5ad(path)
        missing = [c for c in REQUIRED_OBS_COLUMNS if c not in adata.obs.columns]
        if missing:
            raise PipelineError(
                f"'{path}': obs missing required column(s) {missing} -- run `03_spatialdata_import/04_attach_clinical_and_technical_metadata.py` first."
            )

        section_id = adata.obs["section_id"].iloc[0]
        panel_by_section[section_id] = set(adata.var_names)
        # Xenium cell_id values (e.g. 'aaadggoi-1') are only unique within a
        # single section -- prefixing with section_id before concatenation
        # is what makes the combined object's obs_names globally unique.
        adata.obs_names = [f"{section_id}_{cid}" for cid in adata.obs["cell_id"]]
        adatas[section_id] = adata

    core_genes = set.intersection(*panel_by_section.values())
    all_genes = sorted(set.union(*panel_by_section.values()))

    combined = ad.concat(adatas, join="outer", index_unique=None, label="_concat_section_id")
    if combined.obs_names.duplicated().any():
        n_dup = combined.obs_names.duplicated().sum()
        raise PipelineError(
            f"Combined object has {n_dup} duplicate obs_names after prefixing -- globally unique ID construction failed."
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output_path)

    membership = pd.DataFrame(
        {
            gene: {section: gene in genes for section, genes in panel_by_section.items()}
            for gene in all_genes
        }
    ).T
    panel_membership_path.parent.mkdir(parents=True, exist_ok=True)
    membership.to_parquet(panel_membership_path)

    return {
        "n_sections": len(adatas),
        "n_cells_total": combined.n_obs,
        "n_genes_union": combined.n_vars,
        "n_genes_core": len(core_genes),
        "n_patients": combined.obs["patient_id"].nunique(),
    }
