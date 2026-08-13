"""Marker and reference registry compilation (`06_cell_type_annotation/00_compile_marker_and_reference_registry.py`).

Builds a versioned registry mapping the 399-gene `biological_gene`
panel to feasible cell identities, at the resolution the panel actually
supports -- not an aspirational literature taxonomy applied blindly.

Two things were checked against the panel before curating
this registry, not assumed:

1. An initial pass built candidate marker lists by browsing the panel and
   picking plausible-looking genes -- this is methodologically circular
   (every gene "found" was trivially "present" by construction) and was
   discarded in favour of starting from independent, genuine domain
   knowledge of each lineage's canonical markers and checking presence
   against the panel, the same discipline used for `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s program
   scores and the earlier lineage-coverage check in
   Quality Control.

2. The panel contains a complete set of pancreatic islet hormones (INS,
   GCG, PPY, SST, GHRL, PCSK2 -- 6/6 present), plus renal (AQP2, SLC4A1,
   UMOD, AQP8), hepatic (CYP3A4, CYP2B6, APOA5), and melanocyte (MLANA,
   SNCA, SNCG) markers -- genes with no plausible relevance to HNSCC
   oropharynx tumour tissue. This is strong internal evidence that the
   399-gene core is a standard commercial multi-tissue/pan-cancer Xenium
   panel (10x Genomics), not a bespoke HNSCC design, with this study's
   patient-specific TCR/CDR3 probes added on top, consistent with the
   probe-design finding from SpatialData Import. This is an evidence-based inference, not independently
   confirmed against 10x's published panel catalogue -- stated as such.
   A material fraction of the panel is therefore "off-target" for this
   project's actual biology and is excluded from the registry below
   rather than force-mapped to a cell identity.

Confidence tiers reflect genuine lineage specificity in HNSCC tissue, not
merely presence in the panel: e.g. no squamous-specific epithelial markers
(TP63, KRT5, KRT14, KRT17, KRT19, KRT8, KRT18, SFN) are in this panel at
all, so `Epithelial_Tumour` is callable only as a general epithelial/
carcinoma identity, not confidently as squamous-subtype-specific -- a
load-bearing limitation for this project's tumour-focused claims, recorded
here and expected to propagate as an explicit caveat through `06_cell_type_annotation/02_score_major_lineages.py`/
6.06's confidence scoring.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

REGISTRY_VERSION = "v1"

# hierarchy_level: "major_lineage" (`06_cell_type_annotation/02_score_major_lineages.py`'s target compartments) or
# "substate" (finer resolution feeding `06_cell_type_annotation/04_resolve_t_cell_substates.R`, `06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`).
CELL_TYPE_MARKER_REGISTRY: list[dict] = [
    {
        "cell_identity": "Epithelial_Tumour",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["EPCAM", "MET", "ERBB2", "EGFR", "KRT7"],
        "confidence_tier": "moderate",
        "rationale": (
            "General epithelial/carcinoma markers present; no squamous-specific markers "
            "(TP63, KRT5, KRT14, KRT17, KRT19, KRT8, KRT18, SFN) are in this panel -- callable "
            "as epithelial/tumour, not confidently as squamous-subtype-specific."
        ),
    },
    {
        "cell_identity": "T_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": [
            "CD3D",
            "CD3E",
            "CD2",
            "CD7",
            "TRAC",
            "TRBC1",
            "TRBC2",
            "TRDC",
            "CD247",
            "CD28",
            "CD27",
            "CD226",
            "IL7R",
        ],
        "confidence_tier": "high",
        "rationale": "Strong, near-canonical T-cell receptor complex and pan-T-cell marker coverage.",
    },
    {
        "cell_identity": "NK_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["NKG7", "GNLY", "KLRD1", "KLRB1", "KLRC1", "FGFBP2"],
        "confidence_tier": "moderate",
        "rationale": (
            "Cytotoxic/NK receptor markers present, but overlap substantially with cytotoxic "
            "CD8 T cells (GZMB/PRF1 deliberately excluded here as non-specific) -- NK identity "
            "requires CD3-negativity as an explicit exclusion criterion, not markers alone."
        ),
    },
    {
        "cell_identity": "B_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["CD19", "MS4A1", "CD79A", "BANK1", "TCL1A", "IGHD"],
        "confidence_tier": "high",
        "rationale": "Canonical pan-B-cell markers all present.",
    },
    {
        "cell_identity": "Plasma_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["MZB1", "DERL3", "TNFRSF17", "FKBP11", "TNFRSF13B"],
        "confidence_tier": "high",
        "rationale": "Canonical plasma-cell markers all present; distinguishable from B_cell (CD19/MS4A1-low).",
    },
    {
        "cell_identity": "Myeloid",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": [
            "CD68",
            "CD14",
            "CD163",
            "MRC1",
            "MARCO",
            "AIF1",
            "MPEG1",
            "FCGR3A",
            "FCGR1A",
            "MNDA",
            "FCN1",
            "TREM2",
            "CD5L",
            "ADGRE1",
            "VSIG4",
        ],
        "confidence_tier": "high",
        "rationale": "Broad, strong monocyte/macrophage marker coverage; substate resolution (`06_cell_type_annotation/05_resolve_myeloid_and_stromal_substates.R`) needed to separate macrophage polarisation states.",
    },
    {
        "cell_identity": "Dendritic_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["CD1A", "CD1C", "CD1E", "LAMP3", "CLEC10A", "FCER1A", "LILRA4", "IRF8", "SPIB"],
        "confidence_tier": "high",
        "rationale": "Covers conventional (CD1C/CLEC10A/FCER1A), mature/migratory (LAMP3), and plasmacytoid (LILRA4/IRF8/SPIB) DC axes.",
    },
    {
        "cell_identity": "Mast_cell",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["CPA3", "MS4A2", "KIT", "CTSG"],
        "confidence_tier": "high",
        "rationale": "Canonical mast-cell markers all present.",
    },
    {
        "cell_identity": "Fibroblast",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": [
            "PDGFRA",
            "FBLN1",
            "FBN1",
            "ASPN",
            "SFRP2",
            "SFRP4",
            "DPT",
            "COL5A2",
            "MFAP5",
            "OGN",
            "VCAN",
            "TNC",
        ],
        "confidence_tier": "high",
        "rationale": "Strong ECM/stromal marker coverage. PDGFRB and THY1 excluded here (shared with Perivascular_SmoothMuscle / general mesenchymal identity) to keep the fibroblast set specific.",
    },
    {
        "cell_identity": "Endothelial",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": [
            "PECAM1",
            "VWF",
            "EGFL7",
            "CLEC14A",
            "RAMP2",
            "SOX17",
            "SOX18",
            "ERG",
            "MMRN2",
            "ECSCR",
            "GNG11",
        ],
        "confidence_tier": "high",
        "rationale": "Canonical pan-endothelial markers all present.",
    },
    {
        "cell_identity": "Lymphatic_endothelial",
        "hierarchy_level": "substate",
        "parent_identity": "Endothelial",
        "markers": ["PROX1", "LYVE1", "MMRN1"],
        "confidence_tier": "moderate",
        "rationale": "Distinguishes lymphatic from blood endothelium within the Endothelial major lineage.",
    },
    {
        "cell_identity": "Perivascular_SmoothMuscle",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["ACTA2", "MYH11", "MYLK", "RERGL", "HIGD1B", "DES", "CNN1", "PDGFRB"],
        "confidence_tier": "moderate",
        "rationale": (
            "Covers pericyte/smooth-muscle identity; ACTA2 and PDGFRB are also myofibroblast/"
            "activated-stroma markers, so this identity can overlap with an activated Fibroblast "
            "state, not always cleanly separable in a targeted panel."
        ),
    },
    {
        "cell_identity": "Erythroid",
        "hierarchy_level": "major_lineage",
        "parent_identity": None,
        "markers": ["GYPA", "GYPB", "ALAS2", "HEMGN", "AHSP"],
        "confidence_tier": "high",
        "rationale": "Canonical erythroid markers -- expected as a background/blood-contamination signal in tissue sections, not a tumour-microenvironment cell type of biological interest.",
    },
]

# Panel genes with strong, independent evidence of being off-target for
# HNSCC tissue (see module docstring) -- excluded from the registry above
# rather than force-mapped. Recorded explicitly so this is a documented
# exclusion, not a silent gap.
OFF_TARGET_GENE_EXAMPLES = [
    "INS",
    "GCG",
    "PPY",
    "SST",
    "GHRL",
    "PCSK2",  # pancreatic islet
    "AQP2",
    "SLC4A1",
    "UMOD",
    "AQP8",  # renal
    "CYP3A4",
    "CYP2B6",
    "APOA5",  # hepatic
    "MLANA",
    "SNCA",
    "SNCG",  # melanocyte
]


def validate_registry(
    available_genes: set[str], registry: list[dict] = CELL_TYPE_MARKER_REGISTRY
) -> list[dict]:
    """Restricts each identity's markers to genes actually present, raising
    if an identity would lose ALL its markers (a registry/panel mismatch
    serious enough to require attention, not silent omission)."""
    validated = []
    for entry in registry:
        present = [g for g in entry["markers"] if g in available_genes]
        if len(present) == 0:
            raise PipelineError(
                f"Identity '{entry['cell_identity']}' has 0 markers present in the panel."
            )
        validated.append({**entry, "markers": present, "n_markers_in_panel": len(present)})
    return validated


def build_marker_registry_table(available_genes: set[str]) -> pd.DataFrame:
    validated = validate_registry(available_genes)
    rows = []
    for entry in validated:
        rows.append(
            {
                "cell_identity": entry["cell_identity"],
                "hierarchy_level": entry["hierarchy_level"],
                "parent_identity": entry["parent_identity"] or "",
                "markers": ";".join(entry["markers"]),
                "n_markers_in_panel": entry["n_markers_in_panel"],
                "confidence_tier": entry["confidence_tier"],
                "rationale": entry["rationale"],
                "registry_version": REGISTRY_VERSION,
            }
        )
    return pd.DataFrame(rows)


def build_marker_registry_report(project_root: Path) -> dict:
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "references" / "cell_type_marker_registry.tsv"

    if not feature_annotation_path.is_file():
        raise PipelineError(
            f"'{feature_annotation_path}' not found. Run `05_preprocessing_and_normalisation/00_separate_gene_and_control_features.py` first."
        )

    feat = pd.read_csv(feature_annotation_path, sep="\t")
    available_genes = set(feat.loc[feat["feature_class"] == "biological_gene", "feature_name"])

    table = build_marker_registry_table(available_genes)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_path, sep="\t", index=False)

    all_registry_genes = set()
    for entry in CELL_TYPE_MARKER_REGISTRY:
        all_registry_genes.update(entry["markers"])
    unmapped_genes = available_genes - all_registry_genes

    return {
        "n_identities": len(table),
        "n_major_lineages": int((table["hierarchy_level"] == "major_lineage").sum()),
        "n_substates": int((table["hierarchy_level"] == "substate").sum()),
        "confidence_tier_counts": table["confidence_tier"].value_counts().to_dict(),
        "n_panel_genes_total": len(available_genes),
        "n_panel_genes_mapped": len(all_registry_genes & available_genes),
        "n_panel_genes_unmapped": len(unmapped_genes),
        "registry_version": REGISTRY_VERSION,
    }
