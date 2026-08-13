"""Blinded ecosystem annotation (`10_niche_and_ecosystem_discovery/04_annotate_ecosystems_with_blinded_rules.py`).

Assigns a descriptive, human-readable name to each of `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s
unsupervised archetypes -- after discovery, from a documented, mechanical
rubric applied identically to every archetype, not hand-picked per result
and not requiring a domain expert to view tissue images. This is
different in kind from `06_cell_type_annotation/07_blinded_annotation_review.py`'s blinded panel review
: that milestone
required a pathologist's genuine visual judgment of tissue morphology,
which this pipeline has no source for and does not fabricate. Here, the
"rubric" input is itself already a fully computed quantity (`10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R`'s
archetype centroids vs. the dataset-wide mean composition) -- naming
from it is deterministic and fully automatable, closer to the naming
convention used for published cellular-neighbourhood analyses (e.g.
Schurch et al. 2020's CN labels) than to expert histopathological
adjudication.

**Rubric:** for each archetype and each of `06_cell_type_annotation/06_integrate_annotation_evidence.py`'s 12 lineages,
compute the enrichment ratio = archetype's mean neighbour-composition
fraction for that lineage / the dataset-wide mean fraction for that
lineage (both from the same radius_30.0um scale as archetype discovery
itself). Raw composition fraction alone is a poor naming signal because
it is dominated by globally abundant lineages (Epithelial_Tumour is
~35.7% of the average neighbourhood everywhere) regardless of whether an
archetype is differentiated by that lineage; the enrichment
ratio corrects for this baseline. Any lineage with enrichment ratio >=
`ENRICHMENT_THRESHOLD = 2.0` (a standard "at least doubled relative to
background" convention) is included in the archetype's name, ordered by
descending enrichment, joined with "/". An archetype with no lineage
clearing the threshold is labelled "Mixed/non-specific niche" rather
than forced into a spurious specific name.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

ENRICHMENT_THRESHOLD = 2.0
SCALE_PREFIX = "radius_30.0um__"

LINEAGE_DISPLAY_NAMES = {
    "B_cell": "B-cell",
    "Dendritic_cell": "Dendritic-cell",
    "Endothelial": "Endothelial",
    "Epithelial_Tumour": "Tumour",
    "Erythroid": "Erythroid",
    "Fibroblast": "Fibroblast",
    "Mast_cell": "Mast-cell",
    "Myeloid": "Myeloid",
    "NK_cell": "NK-cell",
    "Perivascular_SmoothMuscle": "Perivascular",
    "Plasma_cell": "Plasma-cell",
    "T_cell": "T-cell",
}


def compute_enrichment_ratios(
    centroids: pd.DataFrame, global_mean: pd.Series, lineage_cols: list[str]
) -> pd.DataFrame:
    """Pure, testable: archetype centroid fraction / dataset-wide mean
    fraction, per lineage. `centroids` indexed by archetype id."""
    return centroids[lineage_cols].div(global_mean[lineage_cols], axis=1)


def label_ecosystem(
    enrichment_row: pd.Series,
    threshold: float = ENRICHMENT_THRESHOLD,
    display_names: dict[str, str] = LINEAGE_DISPLAY_NAMES,
) -> str:
    """Pure, testable: builds one archetype's descriptive name from its
    enrichment ratios, per the documented rubric (module docstring)."""
    enriched = enrichment_row[enrichment_row >= threshold].sort_values(ascending=False)
    if len(enriched) == 0:
        return "Mixed/non-specific niche"
    names = [display_names.get(lineage, lineage) for lineage in enriched.index]
    return "/".join(names) + " niche"


def build_ecosystem_annotation(project_root: Path) -> dict:
    compositions_path = project_root / "data" / "derived" / "local_compositions.parquet"
    centroids_path = project_root / "data" / "derived" / "neighbourhood_archetype_centroids.parquet"
    output_path = project_root / "metadata" / "ecosystem_annotation.tsv"

    if not compositions_path.is_file():
        raise PipelineError(
            f"'{compositions_path}' not found. Run `10_niche_and_ecosystem_discovery/01_compute_local_neighbourhood_compositions.py` first."
        )
    if not centroids_path.is_file():
        raise PipelineError(
            f"'{centroids_path}' not found. Run `10_niche_and_ecosystem_discovery/02_discover_neighbourhood_archetypes.R` first."
        )

    compositions = pd.read_parquet(compositions_path)
    scale_cols = [c for c in compositions.columns if c.startswith(SCALE_PREFIX)]
    lineage_cols = [c[len(SCALE_PREFIX) :] for c in scale_cols]
    global_mean = compositions[scale_cols].mean()
    global_mean.index = lineage_cols

    centroids = pd.read_parquet(centroids_path).set_index("archetype")
    if not set(lineage_cols).issubset(centroids.columns):
        raise PipelineError(
            f"'{centroids_path}' is missing expected lineage columns {sorted(set(lineage_cols) - set(centroids.columns))}."
        )

    enrichment = compute_enrichment_ratios(centroids, global_mean, lineage_cols)
    ecosystem_labels = enrichment.apply(label_ecosystem, axis=1)
    top_enriched = enrichment.apply(
        lambda row: ", ".join(
            f"{lineage}={ratio:.2f}x"
            for lineage, ratio in row.sort_values(ascending=False).head(3).items()
        ),
        axis=1,
    )

    result = pd.DataFrame(
        {
            "archetype": centroids.index,
            "ecosystem_label": ecosystem_labels.to_numpy(),
            "dominant_lineage": centroids["dominant_lineage"].to_numpy(),
            "top_enriched_lineages": top_enriched.to_numpy(),
            "n_cells": centroids["n_cells"].to_numpy(),
        }
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_archetypes": len(result),
        "n_mixed_non_specific": int(
            (result["ecosystem_label"] == "Mixed/non-specific niche").sum()
        ),
        "ecosystem_labels": result["ecosystem_label"].tolist(),
        "output_path": str(output_path),
    }
