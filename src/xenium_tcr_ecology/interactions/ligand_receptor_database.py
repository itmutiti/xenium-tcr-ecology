"""Ligand-receptor database filtered to the panel (`14_spatial_interactions_and_barriers/01_filter_ligand_receptor_database_to_panel.py`).

**Scope decision:** this project's analysis matrix is a targeted,
623-gene panel (immuno-oncology/TCR-focused, per `05_preprocessing_and_normalisation/03_calculate_program_scores.py`'s already-documented
coverage assessment), not a whole-transcriptome assay. A comprehensive
genome-wide ligand-receptor database (CellPhoneDB, CellChatDB,
connectomeDB2020 etc. typically list 2,000+ pairs) would have the vast
majority of its pairs immediately discarded as "not in panel" -- a
low-information operation for a panel this size and this specifically
immune-focused. No such generic database file exists in this project,
and fetching one (an external network resource) was judged not worth
the effort for a targeted panel where a smaller, directly relevant,
hand-curated list is more informative anyway.

Instead, this module hand-curates a standard, well-known set of
canonical immune checkpoint/costimulatory/homing ligand-receptor pairs
directly relevant to `14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`'s 4 sender-receiver pairs -- each a
standard pair from canonical immunology (PD-1/PD-L1, CTLA4/CD80-CD86,
TIGIT/PVR, LAG3/MHC-II, etc.), not sourced from a specific external
database file. This is a scoped alternative for a targeted panel,
stated directly rather than presented as a comprehensive database
filter.

Each pair is checked against `metadata/feature_
annotation.tsv`'s `in_analysis_matrix` genes -- a pair with both ligand
and receptor present is `pair_complete = True`; a pair with only one
side present is reportable (a one-sided expression signal can still be
informative) but explicitly labelled `pair_complete = False`, per this
milestone's stated purpose ("prevent overinterpretation" -- an
incomplete pair cannot support a ligand-receptor signalling claim, only
a one-sided expression observation).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

# Standard, canonical immune checkpoint/costimulatory/homing pairs.
# `programs` links each pair back to `14_spatial_interactions_and_barriers/00_define_sender_receiver_pairs.py`'s programme gene sets,
# for direct cross-reference.
CANDIDATE_LR_PAIRS = [
    {"pair_id": "PD1_PDL1", "ligand": "CD274", "receptor": "PDCD1", "programs": ["exhaustion"]},
    {"pair_id": "PD1_PDL2", "ligand": "PDCD1LG2", "receptor": "PDCD1", "programs": ["exhaustion"]},
    {"pair_id": "CTLA4_CD80", "ligand": "CD80", "receptor": "CTLA4", "programs": ["exhaustion"]},
    {"pair_id": "CTLA4_CD86", "ligand": "CD86", "receptor": "CTLA4", "programs": ["exhaustion"]},
    {"pair_id": "CD28_CD80", "ligand": "CD80", "receptor": "CD28", "programs": ["activation"]},
    {"pair_id": "CD28_CD86", "ligand": "CD86", "receptor": "CD28", "programs": ["activation"]},
    {"pair_id": "TIGIT_PVR", "ligand": "PVR", "receptor": "TIGIT", "programs": ["exhaustion"]},
    {"pair_id": "TIGIT_PVRL2", "ligand": "PVRL2", "receptor": "TIGIT", "programs": ["exhaustion"]},
    {
        "pair_id": "LAG3_MHCII",
        "ligand": "HLA-DRA",
        "receptor": "LAG3",
        "programs": ["exhaustion", "antigen_presentation"],
    },
    {"pair_id": "TIM3_GAL9", "ligand": "LGALS9", "receptor": "HAVCR2", "programs": ["exhaustion"]},
    {"pair_id": "CD27_CD70", "ligand": "CD70", "receptor": "CD27", "programs": ["activation"]},
    {
        "pair_id": "4-1BB_4-1BBL",
        "ligand": "TNFSF9",
        "receptor": "TNFRSF9",
        "programs": ["activation"],
    },
    {"pair_id": "ICOS_ICOSLG", "ligand": "ICOSLG", "receptor": "ICOS", "programs": ["activation"]},
    {"pair_id": "CD40_CD40LG", "ligand": "CD40LG", "receptor": "CD40", "programs": ["activation"]},
    {
        "pair_id": "CXCR3_CXCL9",
        "ligand": "CXCL9",
        "receptor": "CXCR3",
        "programs": ["interferon", "chemokine"],
    },
    {
        "pair_id": "CXCR3_CXCL10",
        "ligand": "CXCL10",
        "receptor": "CXCR3",
        "programs": ["interferon", "chemokine"],
    },
    {"pair_id": "CCR7_CCL19", "ligand": "CCL19", "receptor": "CCR7", "programs": ["chemokine"]},
    {"pair_id": "CXCR4_CXCL12", "ligand": "CXCL12", "receptor": "CXCR4", "programs": ["chemokine"]},
    {"pair_id": "CD47_SIRPA", "ligand": "CD47", "receptor": "SIRPA", "programs": []},
    {"pair_id": "CD39_CD73", "ligand": "ENTPD1", "receptor": "NT5E", "programs": ["exhaustion"]},
]


def check_pair_panel_support(pair: dict, panel_genes: set[str]) -> dict:
    """Pure, testable: per-pair panel-support labelling. A pair is
    `pair_complete` only if both ligand and receptor are present panel
    genes."""
    ligand_present = pair["ligand"] in panel_genes
    receptor_present = pair["receptor"] in panel_genes
    return {
        **pair,
        "ligand_present": ligand_present,
        "receptor_present": receptor_present,
        "pair_complete": ligand_present and receptor_present,
    }


def build_panel_supported_interactions(project_root: Path) -> dict:
    feature_annotation_path = project_root / "metadata" / "feature_annotation.tsv"
    output_path = project_root / "references" / "panel_supported_interactions.tsv"

    if not feature_annotation_path.is_file():
        raise PipelineError(f"'{feature_annotation_path}' not found.")

    feature_annotation = pd.read_csv(feature_annotation_path, sep="\t")
    panel_genes = set(
        feature_annotation.loc[feature_annotation["in_analysis_matrix"], "feature_name"]
    )

    rows = [check_pair_panel_support(pair, panel_genes) for pair in CANDIDATE_LR_PAIRS]
    result = pd.DataFrame(rows)
    result["programs"] = result["programs"].apply(lambda p: ";".join(p))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, sep="\t", index=False)

    return {
        "n_pairs_total": len(result),
        "n_pairs_complete": int(result["pair_complete"].sum()),
        "n_pairs_incomplete": int((~result["pair_complete"]).sum()),
        "output_path": str(output_path),
    }
