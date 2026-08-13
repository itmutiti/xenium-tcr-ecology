"""Unit tests for xenium_tcr_ecology.infra -- the adopted-from-legacy-audit
utilities. These are the only functions in the package with real
implementations at scaffold time; everything else raises
NotImplementedError deliberately."""

from __future__ import annotations

import json
import os

import pytest

from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.inventory import InventoryWriter
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import ENV_ROOT_VAR, PROJECT_ROOT_MARKER, find_project_root


class TestFindProjectRoot:
    def test_finds_root_via_marker_walkup(self, tmp_path):
        root = tmp_path / "repo"
        (root / "manifests").mkdir(parents=True)
        (root / PROJECT_ROOT_MARKER).write_text("data: {}\n")
        nested = root / "scripts" / "00_literature"
        nested.mkdir(parents=True)

        found = find_project_root(start=nested / "some_script.py")
        assert found == root

    def test_cli_arg_takes_precedence(self, tmp_path, monkeypatch):
        root = tmp_path / "repo"
        (root / "manifests").mkdir(parents=True)
        (root / PROJECT_ROOT_MARKER).write_text("data: {}\n")
        monkeypatch.delenv(ENV_ROOT_VAR, raising=False)

        found = find_project_root(cli_arg=str(root))
        assert found == root

    def test_cli_arg_without_marker_raises(self, tmp_path):
        with pytest.raises(PipelineError, match="does not contain"):
            find_project_root(cli_arg=str(tmp_path))

    def test_env_var_used_when_no_cli_arg(self, tmp_path, monkeypatch):
        root = tmp_path / "repo"
        (root / "manifests").mkdir(parents=True)
        (root / PROJECT_ROOT_MARKER).write_text("data: {}\n")
        monkeypatch.setenv(ENV_ROOT_VAR, str(root))

        found = find_project_root()
        assert found == root

    def test_no_marker_anywhere_raises(self, tmp_path, monkeypatch):
        monkeypatch.delenv(ENV_ROOT_VAR, raising=False)
        isolated = tmp_path / "no_repo_here"
        isolated.mkdir()
        with pytest.raises(PipelineError, match="Could not locate"):
            find_project_root(start=isolated / "script.py")


class TestInventoryWriter:
    def test_writes_header_once_and_appends(self, tmp_path):
        inv_path = tmp_path / "inventory.tsv"
        writer = InventoryWriter(inv_path, project_root=tmp_path, fields=["a", "b"])
        writer.write_row(a="1", b="2")
        writer.write_row(a="3", b="4")

        lines = inv_path.read_text().splitlines()
        assert lines[0] == "a\tb"
        assert lines[1] == "1\t2"
        assert lines[2] == "3\t4"

    def test_relativizes_paths_to_project_root(self, tmp_path):
        inv_path = tmp_path / "inventory.tsv"
        writer = InventoryWriter(inv_path, project_root=tmp_path, fields=["path"])
        writer.write_row(path=tmp_path / "data" / "raw" / "file.tar")

        line = inv_path.read_text().splitlines()[1]
        assert line == os.path.join("data", "raw", "file.tar")

    def test_unknown_field_raises(self, tmp_path):
        writer = InventoryWriter(tmp_path / "inventory.tsv", project_root=tmp_path, fields=["a"])
        with pytest.raises(PipelineError, match="Unknown inventory field"):
            writer.write_row(a="1", unexpected="2")


class TestJsonRunLogger:
    def test_write_produces_valid_json_with_expected_keys(self, tmp_path):
        logger = JsonRunLogger(
            logs_dir=tmp_path / "logs",
            script_name="test_script",
            project_root=tmp_path,
            phase="00_literature",
        )
        logger.log_event(detail="did a thing")
        logger.log_error("something went wrong")
        out_path = logger.write(status="failed")

        record = json.loads(out_path.read_text())
        assert record["script"] == "test_script"
        assert record["phase"] == "00_literature"
        assert record["status"] == "failed"
        assert record["errors"] == ["something went wrong"]
        assert len(record["events"]) == 1
        assert record["events"][0]["detail"] == "did a thing"
        assert "environment" in record and "python_version" in record["environment"]
