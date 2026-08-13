"""CDR3 cross-patient similarity screen (`08_tcr_clonal_analysis/05_screen_cdr3_cross_patient_similarity.py`).

Screens all pairs of probed CDR3 amino acid sequences (`08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py`'s
registry) for high sequence similarity within the same TCR chain (TRA-TRA
or TRB-TRB only -- a cross-chain comparison is not biologically
meaningful), flagging pairs similar enough to plausibly cross-react at
the probe-hybridisation level. Cross-patient similar pairs are the
specific concern this milestone exists to surface: two patients'
independently-selected clonotypes sharing a near-identical CDR3 is
consistent with an expected phenomenon documented in McCord et al. 2026
(Sci Immunol 11:eaec3133, Figure 3: "VDJdb-matched microbial-reactive
(EBV, influenza) TCRs are used as bystander-clone markers") --
public/quasi-public viral-reactive
TCR motifs are recurrent across unrelated individuals, not a probe-design
error, so a cross-patient hit is not automatically "wrong," but it is a
risk factor for probe cross-reactivity that should be surfaced and
cross-referenced against `08_tcr_clonal_analysis/01_map_tcr_probes_to_patients.py`'s patient-specificity results and
`08_tcr_clonal_analysis/04_estimate_false_positive_tcr_calls.R`'s false-positive-rate estimates.

**Method:** plain Levenshtein (edit) distance, implemented directly
(pure Python DP) rather than adding a new dependency (`python-Levenshtein`/
`rapidfuzz` are not installed, and at this scale -- 216 probes, at most
C(118,2) = 6,903 same-chain pairs -- a hand-written O(n*m) implementation
is fast enough and avoids an unnecessary dependency for such a small,
well-understood algorithm). `MAX_SIMILARITY_EDIT_DISTANCE = 2`: checked
against the CDR3 length distribution (9-19 amino acids, median
13-14) before choosing -- an edit distance of <=2 on a ~13-14 residue
loop is >=85% sequence identity, a conservative, standard
"near-identical, plausible cross-reactivity risk" bar for CDR3 sequences,
not a loose one that would flag every pair.
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import pandas as pd

from xenium_tcr_ecology.infra.exceptions import PipelineError

MAX_SIMILARITY_EDIT_DISTANCE = 2


def levenshtein_distance(a: str, b: str) -> int:
    """Standard O(len(a) * len(b)) edit-distance dynamic program."""
    if a == b:
        return 0
    if len(a) == 0:
        return len(b)
    if len(b) == 0:
        return len(a)

    prev_row = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        curr_row = [i] + [0] * len(b)
        for j, cb in enumerate(b, start=1):
            cost = 0 if ca == cb else 1
            curr_row[j] = min(
                prev_row[j] + 1,  # deletion
                curr_row[j - 1] + 1,  # insertion
                prev_row[j - 1] + cost,  # substitution
            )
        prev_row = curr_row
    return prev_row[-1]


def screen_cdr3_pairwise_similarity(
    registry: pd.DataFrame, max_distance: int = MAX_SIMILARITY_EDIT_DISTANCE
) -> pd.DataFrame:
    """Pure, testable pairwise screen over a registry-shaped DataFrame
    with `probe_name`, `cdr3_amino_acid_sequence`, `tcr_chain`,
    `patients_with_probe` (semicolon-joined) columns. Only same-chain
    pairs are compared."""
    rows = []
    for chain, group in registry.groupby("tcr_chain"):
        records = group[
            ["probe_name", "cdr3_amino_acid_sequence", "patients_with_probe"]
        ].to_records(index=False)
        for (name_a, seq_a, patients_a), (name_b, seq_b, patients_b) in combinations(records, 2):
            distance = levenshtein_distance(seq_a, seq_b)
            if distance > max_distance:
                continue
            patients_a_set = set(str(patients_a).split(";"))
            patients_b_set = set(str(patients_b).split(";"))
            rows.append(
                {
                    "probe_a": name_a,
                    "probe_b": name_b,
                    "tcr_chain": chain,
                    "cdr3_a": seq_a,
                    "cdr3_b": seq_b,
                    "edit_distance": distance,
                    "is_cross_patient": patients_a_set.isdisjoint(patients_b_set),
                }
            )
    return pd.DataFrame(rows)


def build_cdr3_similarity_screen(project_root: Path) -> dict:
    registry_path = project_root / "metadata" / "tcr_probe_registry.tsv"
    output_path = project_root / "reports" / "tcr" / "cdr3_similarity_screen.tsv"

    if not registry_path.is_file():
        raise PipelineError(
            f"'{registry_path}' not found. Run `08_tcr_clonal_analysis/00_identify_tcr_cdr3_probe_features.py` first."
        )

    registry = pd.read_csv(registry_path, sep="\t")
    screen = screen_cdr3_pairwise_similarity(registry)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    screen.to_csv(output_path, sep="\t", index=False)

    n_pairs_total = sum(
        1 for _ in combinations(range(len(registry[registry["tcr_chain"] == "TRA"])), 2)
    ) + sum(1 for _ in combinations(range(len(registry[registry["tcr_chain"] == "TRB"])), 2))

    return {
        "n_probes": len(registry),
        "n_same_chain_pairs_screened": n_pairs_total,
        "n_similar_pairs_flagged": len(screen),
        "n_cross_patient_similar_pairs": (
            int(screen["is_cross_patient"].sum()) if len(screen) else 0
        ),
        "n_probes_involved_in_cross_patient_similarity": (
            len(
                set(screen.loc[screen["is_cross_patient"], "probe_a"])
                | set(screen.loc[screen["is_cross_patient"], "probe_b"])
            )
            if len(screen)
            else 0
        ),
        "max_similarity_edit_distance": MAX_SIMILARITY_EDIT_DISTANCE,
        "output_path": str(output_path),
    }
