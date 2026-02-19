#!/usr/bin/env python3
"""Step 0: Configuration validation and environment setup.

- Parses config.yaml
- Validates all input file paths exist
- Creates output directory structure
- Logs environment info and config snapshot
"""

import json
import platform
import shutil
import sys
from datetime import datetime
from pathlib import Path

import yaml

# Add Scripts dir to path for utils imports
SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS_DIR))

from utils.parsers import load_config, resolve_path


def validate_file(path: Path, description: str) -> None:
    """Raise FileNotFoundError if path does not exist."""
    if not path.exists():
        raise FileNotFoundError(f"{description}: {path}")


def validate_sheet(filepath: Path, sheet_name: str, description: str) -> None:
    """Verify an Excel sheet exists within a workbook (without loading data)."""
    import openpyxl
    wb = openpyxl.load_workbook(filepath, read_only=True, data_only=True)
    available = wb.sheetnames
    wb.close()
    if sheet_name not in available:
        raise ValueError(
            f"{description}: sheet '{sheet_name}' not found in {filepath.name}. "
            f"Available sheets: {available}"
        )


def get_package_versions() -> dict:
    """Get versions of key dependencies."""
    versions = {}
    packages = [
        "pandas", "openpyxl", "numpy", "scipy", "matplotlib", "seaborn",
        "yaml", "requests", "gprofiler", "mygene", "pyteomics", "Bio",
        "gseapy", "matplotlib_venn", "upsetplot",
    ]
    for pkg in packages:
        try:
            mod = __import__(pkg)
            versions[pkg] = getattr(mod, "__version__", "installed (no version)")
        except ImportError:
            versions[pkg] = "NOT INSTALLED"
    return versions


def main():
    config_path = SCRIPTS_DIR / "config.yaml"
    config = load_config(config_path)
    raw_dir = (SCRIPTS_DIR / config["paths"]["raw_dir"]).resolve()

    print(f"[00_setup] Config loaded from {config_path}")
    print(f"[00_setup] Raw data directory: {raw_dir}")

    # ------------------------------------------------------------------
    # 1. Validate all dataset file paths
    # ------------------------------------------------------------------
    print("\n--- Validating dataset files ---")
    for ds_id, ds_cfg in config["datasets"].items():
        filepath = resolve_path(raw_dir, ds_cfg["file"])
        validate_file(filepath, f"Dataset {ds_id}")

        # For Excel files with sheet spec, validate sheet exists
        sheet = ds_cfg.get("sheet")
        if sheet and filepath.suffix in (".xlsx", ".xls"):
            validate_sheet(filepath, sheet, f"Dataset {ds_id}")

        print(f"  {ds_id}: OK ({filepath.name})")

    # EV dataset has multiple sheets
    ev_cfg = config["datasets"].get("EV", {})
    if "sheets" in ev_cfg:
        ev_path = resolve_path(raw_dir, ev_cfg["file"])
        for sheet_key, sheet_name in ev_cfg["sheets"].items():
            validate_sheet(ev_path, sheet_name, f"EV sheet '{sheet_key}'")

    # ------------------------------------------------------------------
    # 2. Validate reference files
    # ------------------------------------------------------------------
    print("\n--- Validating reference files ---")
    for ref_id, ref_cfg in config["references"].items():
        filepath = resolve_path(raw_dir, ref_cfg["file"])
        validate_file(filepath, f"Reference {ref_id}")

        for sheet_key in ["deduplicated_sheet", "full_sheet", "sheet"]:
            sheet = ref_cfg.get(sheet_key)
            if sheet and filepath.suffix in (".xlsx", ".xls"):
                validate_sheet(filepath, sheet, f"Reference {ref_id} ({sheet_key})")

        print(f"  {ref_id}: OK ({filepath.name})")

    # ------------------------------------------------------------------
    # 3. Create output directory structure
    # ------------------------------------------------------------------
    print("\n--- Creating directory structure ---")
    base_dir = SCRIPTS_DIR.parent  # DataProc/
    dirs_to_create = [
        "DerivedData/standardised",
        "DerivedData/orthology_cache",
        "DerivedData/uniprot_cache",
        "DerivedData/evidence_scores",
        "DerivedData/autophagy_panel",
        "DerivedData/peptide_feasibility",
        "DerivedData/modules",
        "QC",
        "Outputs/figures",
        "Log",
    ]
    for d in dirs_to_create:
        dirpath = base_dir / d
        dirpath.mkdir(parents=True, exist_ok=True)
        print(f"  {d}/")

    # ------------------------------------------------------------------
    # 4. Log environment info
    # ------------------------------------------------------------------
    print("\n--- Logging environment info ---")
    log_dir = (SCRIPTS_DIR / config["paths"]["log_dir"]).resolve()
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    log_entry = {
        "timestamp": datetime.now().isoformat(),
        "python_version": sys.version,
        "platform": platform.platform(),
        "package_versions": get_package_versions(),
        "config": config,
    }

    log_path = log_dir / f"setup_{timestamp}.json"
    with open(log_path, "w") as f:
        json.dump(log_entry, f, indent=2, default=str)
    print(f"  Environment log: {log_path}")

    # Config snapshot
    snapshot_path = log_dir / "config_snapshot.yaml"
    shutil.copy2(config_path, snapshot_path)
    print(f"  Config snapshot: {snapshot_path}")

    # ------------------------------------------------------------------
    # 5. Summary
    # ------------------------------------------------------------------
    n_datasets = len(config["datasets"])
    n_refs = len(config["references"])
    print(f"\n[00_setup] DONE: {n_datasets} datasets, {n_refs} references validated.")
    print(f"[00_setup] Directory structure ready.")

    # Check optional external files
    plasma_ext = config.get("plasma_proteins", {}).get("external_list")
    if plasma_ext:
        ext_path = resolve_path(raw_dir, plasma_ext)
        if not ext_path.exists():
            print(f"  WARNING: Optional external plasma list not found: {ext_path}")

    return config


if __name__ == "__main__":
    main()
