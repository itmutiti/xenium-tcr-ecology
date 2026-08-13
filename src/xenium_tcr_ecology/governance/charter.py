"""Project charter compilation (`01_project_setup_and_governance/00_validate_project_scope.py`).

Records the project's primary research questions and how they relate to
the source paper's own published findings. This is a factual snapshot,
not a gate: generating this charter, and every downstream computational
stage, proceeds unconditionally.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import yaml

PRIMARY_QUESTIONS = {
    "Q1": "Can clone-resolved spatial ecology be measured with statistically calibrated, "
    "generalisable methods, and does an open framework built on this cohort transfer to an "
    "independent spatial dataset?",
    "Q2": "How much of the variability in T-cell clonal spatial and phenotypic behaviour is "
    "attributable to intrinsic clonal/transcriptional identity, versus local microenvironmental "
    "context, versus patient identity -- and is that partition stable in an independent cohort?",
    "Q3": "Does fibroblast/myeloid barrier topology explain residual variance in clone-tumour "
    "engagement beyond cell-intrinsic state and niche composition, and how does this compare "
    "quantitatively to published CAF-exclusion mechanisms in other cancers?",
}


def build_project_charter(project_root: Path, output_path: Path) -> dict:
    """Build project_charter.yaml. Always succeeds."""
    charter = {
        "project_name": "xenium-tcr-ecology",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "primary_questions": PRIMARY_QUESTIONS,
        "source_dataset": "GSE300147 (McCord et al., Science Immunology 11, eaec3133, 2026)",
        "published_findings_vs_original_hypotheses": (
            "The source paper already shows spatial context shapes clonal T-cell state "
            "(qualitative, illustrative-clone-level). This project's Q1-Q3 are framed as "
            "quantitative, calibrated, and externally validated extensions of that finding, "
            "not a restatement of it."
        ),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.dump(charter, sort_keys=False, width=100))
    return charter
