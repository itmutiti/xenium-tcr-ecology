#!/usr/bin/env python3
"""
`03_spatialdata_import/01_import_each_section_to_spatialdata.py`

Runs the custom Xenium reader (`03_spatialdata_import/00_detect_xenium_format_version.py` established the standard
spatialdata_io.xenium() reader cannot be used on this dataset) across every
standardised section, writing one SpatialData Zarr store per section.

Idempotent: an existing, non-empty <section>.zarr is skipped unless --force
is passed, since each section takes substantial time (large morphology image
decompression + polygon construction over tens of thousands of cells).

Primary output: data/objects/spatialdata/<section>.zarr
"""

from __future__ import annotations

import shutil
import sys
import time

from xenium_tcr_ecology.cli import base_parser
from xenium_tcr_ecology.infra.exceptions import PipelineError
from xenium_tcr_ecology.infra.logging import JsonRunLogger
from xenium_tcr_ecology.infra.paths import find_project_root
from xenium_tcr_ecology.io.xenium_reader import import_section


def main() -> int:
    parser = base_parser(__doc__)
    parser.add_argument(
        "--force", action="store_true", help="Re-import sections that already have a .zarr store."
    )
    args = parser.parse_args()

    try:
        project_root = find_project_root(args.project_root)
    except PipelineError as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1

    standardised_root = project_root / "data" / "standardised"
    output_root = project_root / "data" / "objects" / "spatialdata"
    output_root.mkdir(parents=True, exist_ok=True)

    section_dirs = sorted(p for p in standardised_root.iterdir() if p.is_dir())
    if not section_dirs:
        print(
            f"[ERROR] No section directories found under '{standardised_root}'. Run `02_raw_data_ingestion/05_standardise_sample_directory_layout.py` first.",
            file=sys.stderr,
        )
        return 1

    logger = JsonRunLogger(
        logs_dir=project_root / "results" / "logs" / "03_spatialdata_import",
        script_name="01_import_each_section_to_spatialdata",
        project_root=project_root,
        phase="03_spatialdata_import",
    )

    imported, skipped, failed = 0, 0, []
    for section_dir in section_dirs:
        section_id = section_dir.name
        zarr_path = output_root / f"{section_id}.zarr"

        if zarr_path.exists() and not args.force:
            print(f"[SKIP]  {section_id}: {zarr_path} already exists.")
            skipped += 1
            continue

        t0 = time.time()
        try:
            sdata = import_section(section_dir, section_id)
            if zarr_path.exists():
                shutil.rmtree(zarr_path)
            sdata.write(zarr_path)
        except (
            Exception
        ) as exc:  # noqa: BLE001 -- surfaced per-section, doesn't abort the whole run
            print(f"[ERROR] {section_id}: {exc}", file=sys.stderr)
            logger.log_error(f"{section_id}: {exc}")
            failed.append(section_id)
            continue

        elapsed = time.time() - t0
        n_cells = sdata["table"].n_obs
        print(f"[OK]   {section_id}: {n_cells} cells, wrote {zarr_path} ({elapsed:.1f}s)")
        logger.log_event(section_id=section_id, n_cells=n_cells, elapsed_s=round(elapsed, 1))
        imported += 1

    status = "ok" if not failed else "failed"
    logger.log_event(imported=imported, skipped=skipped, failed=failed)
    logger.write(status=status)

    print(f"\n[SUMMARY] {imported} imported, {skipped} skipped, {len(failed)} failed.")
    if failed:
        print(f"[ERROR] Failed sections: {failed}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
