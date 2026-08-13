"""Unit tests for xenium_tcr_ecology.ingest (`01_project_setup_and_governance/01_build_sample_manifest.py` sample manifest)."""

from __future__ import annotations

import pytest
import yaml

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.ingest.manifest import compile_sample_manifest


def _sample(patient_number, run_number, specimen_type="HNSCC tumor", p16="Positive"):
    return {
        "gsm_accession": f"GSM_{patient_number}_{run_number}",
        "patient_number": patient_number,
        "run_number": run_number,
        "specimen_type": specimen_type,
        "p16_ihc_status": p16,
        "recurrence_status": "Primary",
        "smoking_pack_years": 10,
        "tumour_resection_site": "Oropharynx",
    }


class TestCompileSampleManifest:
    def test_derives_replicate_and_cohort_flags_correctly(self, tmp_path):
        samples = [
            _sample(1, 1),  # single-run HNSCC patient
            _sample(2, 1, p16="Negative"),  # replicated HNSCC patient, run 1
            _sample(2, 2, p16="Negative"),  # replicated HNSCC patient, run 2
            _sample(3, 1, specimen_type="Ameloblastoma", p16="N/A"),  # non-HNSCC
        ]
        input_path = tmp_path / "input.yaml"
        input_path.write_text(yaml.dump({"samples": samples}))
        output_path = tmp_path / "manifest.tsv"

        summary = compile_sample_manifest(input_path, output_path, project_root=tmp_path)

        assert summary["total_samples"] == 4
        assert summary["hnscc_patients"] == 2
        assert summary["hnscc_sections"] == 3
        assert summary["replicated_patients"] == 1
        assert summary["hpv_positive_patients"] == 1
        assert summary["ameloblastoma_specimens"] == 1

        rows = output_path.read_text().splitlines()
        header = rows[0].split("\t")
        data_rows = [dict(zip(header, r.split("\t"))) for r in rows[1:]]

        p1 = next(r for r in data_rows if r["patient_id"] == "P01")
        assert p1["is_technical_replicate"] == "False"
        assert p1["included_in_primary_hnscc_cohort"] == "True"

        p2_run1 = next(r for r in data_rows if r["section_id"] == "P02_run1")
        assert p2_run1["is_technical_replicate"] == "True"
        assert p2_run1["hpv_p16_positive"] == "False"

        p3 = next(r for r in data_rows if r["patient_id"] == "P03")
        assert p3["included_in_primary_hnscc_cohort"] == "False"

    def test_raises_on_duplicate_gsm(self, tmp_path):
        samples = [_sample(1, 1), _sample(1, 1)]
        for s in samples:
            s["gsm_accession"] = "GSM_DUPLICATE"
        input_path = tmp_path / "input.yaml"
        input_path.write_text(yaml.dump({"samples": samples}))

        with pytest.raises(PipelineError, match="Duplicate GSM accession"):
            compile_sample_manifest(input_path, tmp_path / "out.tsv", project_root=tmp_path)

    def test_raises_on_missing_input(self, tmp_path):
        with pytest.raises(PipelineError, match="not found"):
            compile_sample_manifest(
                tmp_path / "missing.yaml", tmp_path / "out.tsv", project_root=tmp_path
            )


class TestRealProjectManifest:
    def test_real_geo_input_matches_known_cohort_structure(self):
        """Exercises the actual config/geo/sample_manifest_input.yaml checked
        into this project (fetched from real GSM records), not just a
        synthetic fixture -- confirms it independently reproduces the cohort
        structure already documented from the source paper (10 patients, 5
        HPV+, 7 replicated = 17 sections, 1 ameloblastoma), which is exactly
        the cross-check this script exists to provide."""
        project_root = find_project_root()
        input_path = project_root / "config" / "geo" / "sample_manifest_input.yaml"
        output_path = project_root / "metadata" / "_test_sample_manifest.tsv"

        summary = compile_sample_manifest(input_path, output_path, project_root=project_root)
        output_path.unlink()

        assert summary["total_samples"] == 18
        assert summary["hnscc_patients"] == 10
        assert summary["hnscc_sections"] == 17
        assert summary["replicated_patients"] == 7
        assert summary["hpv_positive_patients"] == 5
        assert summary["ameloblastoma_specimens"] == 1
