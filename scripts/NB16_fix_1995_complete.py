#!/usr/bin/env python3
"""
NB16_fix_1995_complete.py
=========================
Fix the corrupted 1995 SST file and recompute all Phase 1 results
with the complete 33-year record.

This script:
  1. Re-downloads 1995 in half-year chunks via copernicusmarine CLI
  2. Merges the two halves with cdo mergetime
  3. Regrids the merged 1995 file to 0.125 deg monthly
  4. Re-runs NB14 (full Phase 1 analysis) with all 33 years

IMPORTANT: This script must be run from a Gadi LOGIN NODE (not a PBS
compute node) because copernicusmarine requires interactive network
access for authentication. Run inside a screen session:

  screen -S fix_1995
  module purge && unset PYTHONPATH
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  python NB16_fix_1995_complete.py

After starting, detach with Ctrl+A then D.

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import subprocess
import sys
import os
import time as _time

# =============================================================================
# CONFIGURATION
# =============================================================================
SST_DIR = Path("/g/data/gv90/xl1657/cmems_sst")
NOTEBOOKS_DIR = Path("/g/data/gv90/xl1657/cmems_adt/notebooks")
REGRID_DIR = Path("/g/data/gv90/xl1657/cmems_adt/meander_sst_project/sst_regridded")

YEAR = 1995
FINAL_FILE = SST_DIR / f"cmems_ostia_sst_so30S_{YEAR}.nc"
H1_FILE = SST_DIR / f"cmems_ostia_sst_so30S_{YEAR}_h1.nc"
H2_FILE = SST_DIR / f"cmems_ostia_sst_so30S_{YEAR}_h2.nc"
REGRIDDED_FILE = REGRID_DIR / f"sst_monthly_0125deg_{YEAR}.nc"


def run_cmd(cmd, description):
    """Run a shell command and check for success."""
    print(f"\n  {description}...")
    print(f"    $ {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"    FAILED (exit code {result.returncode})")
        if result.stderr:
            print(f"    stderr: {result.stderr[:500]}")
        return False
    if result.stdout:
        print(f"    stdout: {result.stdout[:300]}")
    return True


# =============================================================================
# STEP 1: RE-DOWNLOAD 1995 IN HALF-YEAR CHUNKS
# =============================================================================
def step1_download():
    """Download 1995 in two half-year chunks to avoid login node memory limits."""
    print("\n" + "=" * 70)
    print("Step 1: Re-download 1995 SST in half-year chunks")
    print("=" * 70)

    # Remove any existing corrupted or partial files
    for fp in [FINAL_FILE, H1_FILE, H2_FILE]:
        if fp.exists():
            print(f"  Removing existing file: {fp.name}")
            fp.unlink()

    # Remove orphaned temp files
    for fp in SST_DIR.glob(f"cmems_ostia_sst_so30S_{YEAR}.nc.*"):
        print(f"  Removing orphan: {fp.name}")
        fp.unlink()

    # Download first half (Jan-Jun)
    cmd_h1 = (
        f"copernicusmarine subset "
        f"--dataset-id METOFFICE-GLO-SST-L4-REP-OBS-SST "
        f"--variable analysed_sst "
        f"--start-datetime '{YEAR}-01-01T00:00:00' "
        f"--end-datetime '{YEAR}-06-30T23:59:59' "
        f"--minimum-latitude -90.0 --maximum-latitude -30.0 "
        f"--minimum-longitude -180.0 --maximum-longitude 179.95 "
        f"--output-directory {SST_DIR} "
        f"--output-filename 'cmems_ostia_sst_so30S_{YEAR}_h1.nc'"
    )
    if not run_cmd(cmd_h1, f"Downloading {YEAR} first half (Jan-Jun)"):
        return False

    if not H1_FILE.exists():
        print(f"    ERROR: {H1_FILE.name} not created")
        return False
    print(f"    H1 size: {H1_FILE.stat().st_size / 1024**3:.1f} GB")

    # Download second half (Jul-Dec)
    cmd_h2 = (
        f"copernicusmarine subset "
        f"--dataset-id METOFFICE-GLO-SST-L4-REP-OBS-SST "
        f"--variable analysed_sst "
        f"--start-datetime '{YEAR}-07-01T00:00:00' "
        f"--end-datetime '{YEAR}-12-31T23:59:59' "
        f"--minimum-latitude -90.0 --maximum-latitude -30.0 "
        f"--minimum-longitude -180.0 --maximum-longitude 179.95 "
        f"--output-directory {SST_DIR} "
        f"--output-filename 'cmems_ostia_sst_so30S_{YEAR}_h2.nc'"
    )
    if not run_cmd(cmd_h2, f"Downloading {YEAR} second half (Jul-Dec)"):
        return False

    if not H2_FILE.exists():
        print(f"    ERROR: {H2_FILE.name} not created")
        return False
    print(f"    H2 size: {H2_FILE.stat().st_size / 1024**3:.1f} GB")

    return True


# =============================================================================
# STEP 2: MERGE HALF-YEAR FILES WITH CDO
# =============================================================================
def step2_merge():
    """Merge the two half-year files into a single annual file."""
    print("\n" + "=" * 70)
    print("Step 2: Merge half-year files with cdo")
    print("=" * 70)

    if not H1_FILE.exists() or not H2_FILE.exists():
        print("  ERROR: Half-year files not found")
        return False

    # Load cdo module and merge
    cmd = (
        f"module load cdo && "
        f"cdo mergetime {H1_FILE} {H2_FILE} {FINAL_FILE}"
    )
    if not run_cmd(cmd, "Merging with cdo mergetime"):
        # cdo may not be available via module in the venv; try direct path
        cmd_alt = f"cdo mergetime {H1_FILE} {H2_FILE} {FINAL_FILE}"
        if not run_cmd(cmd_alt, "Merging with cdo (direct)"):
            return False

    if not FINAL_FILE.exists():
        print("  ERROR: Merged file not created")
        return False

    size_gb = FINAL_FILE.stat().st_size / 1024**3
    print(f"  Merged file size: {size_gb:.1f} GB")

    # Verify time steps
    import xarray as xr
    ds = xr.open_dataset(FINAL_FILE)
    nt = ds.sizes.get("time", 0)
    ds.close()
    expected = 365  # 1995 is not a leap year
    print(f"  Time steps: {nt} (expected {expected})")

    if abs(nt - expected) > 1:
        print(f"  WARNING: Unexpected number of time steps")
        return False

    # Clean up half-year files
    print("  Removing half-year files...")
    H1_FILE.unlink()
    H2_FILE.unlink()
    print("  Cleanup done.")

    return True


# =============================================================================
# STEP 3: REGRID 1995 TO 0.125 DEG MONTHLY
# =============================================================================
def step3_regrid():
    """Regrid the merged 1995 file using NB12."""
    print("\n" + "=" * 70)
    print("Step 3: Regrid 1995 to 0.125 deg monthly")
    print("=" * 70)

    # Remove existing regridded file if present
    if REGRIDDED_FILE.exists():
        REGRIDDED_FILE.unlink()
        print(f"  Removed existing {REGRIDDED_FILE.name}")

    # Run NB12 for 1995 only
    cmd = (
        f"cd {NOTEBOOKS_DIR} && "
        f"python NB12_regrid_sst.py 1995"
    )
    if not run_cmd(cmd, "Running NB12_regrid_sst.py for 1995"):
        return False

    if not REGRIDDED_FILE.exists():
        print("  ERROR: Regridded file not created")
        return False

    size_mb = REGRIDDED_FILE.stat().st_size / 1024**2
    print(f"  Regridded file size: {size_mb:.0f} MB")

    # Verify
    import xarray as xr
    ds = xr.open_dataset(REGRIDDED_FILE)
    nt = ds.sizes.get("time", 0)
    ds.close()
    print(f"  Time steps: {nt} months (expected 12)")

    return nt == 12


# =============================================================================
# STEP 4: RE-RUN FULL PHASE 1 ANALYSIS WITH 33 YEARS
# =============================================================================
def step4_rerun_analysis():
    """Re-run NB14 to recompute all trends with the complete 33-year record."""
    print("\n" + "=" * 70)
    print("Step 4: Re-run Phase 1 analysis with all 33 years")
    print("=" * 70)

    # Verify all 33 regridded files exist
    n_files = len(list(REGRID_DIR.glob("sst_monthly_0125deg_*.nc")))
    print(f"  Regridded files available: {n_files}")

    if n_files < 33:
        print(f"  WARNING: Expected 33 files, found {n_files}")
        missing = []
        for y in range(1993, 2026):
            fp = REGRID_DIR / f"sst_monthly_0125deg_{y}.nc"
            if not fp.exists():
                missing.append(y)
        if missing:
            print(f"  Missing years: {missing}")
        if n_files < 32:
            print("  ERROR: Too few files to proceed")
            return False

    # Submit NB14 as a PBS job since it requires 2.5+ hours of compute
    pbs_script = NOTEBOOKS_DIR / "run_nb14_33yr.pbs"
    pbs_content = """#!/bin/bash
#PBS -N nb14_33yr
#PBS -P gv90
#PBS -q normal
#PBS -l walltime=06:00:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -l storage=gdata/gv90
#PBS -o /g/data/gv90/xl1657/cmems_adt/notebooks/nb14_33yr.out
#PBS -e /g/data/gv90/xl1657/cmems_adt/notebooks/nb14_33yr.err

module purge
unset PYTHONPATH
source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
cd /g/data/gv90/xl1657/cmems_adt/notebooks/
python NB14_fix_1995_and_par_transect.py
"""
    with open(pbs_script, "w") as f:
        f.write(pbs_content)

    result = subprocess.run(f"qsub {pbs_script}", shell=True,
                           capture_output=True, text=True)
    if result.returncode == 0:
        job_id = result.stdout.strip()
        print(f"  Submitted NB14 PBS job: {job_id}")
        print(f"  Monitor: tail -30 {NOTEBOOKS_DIR}/nb14_33yr.out")
        return True
    else:
        print(f"  Failed to submit PBS job: {result.stderr}")
        return False


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB16: Fix 1995 SST + Recompute Complete 33-Year Phase 1 Analysis")
    print("=" * 70)

    t0 = _time.time()

    # Steps 1-3 run on the login node (require network access)
    ok = step1_download()
    if not ok:
        print("\nStep 1 FAILED. Aborting.")
        sys.exit(1)

    ok = step2_merge()
    if not ok:
        print("\nStep 2 FAILED. Aborting.")
        sys.exit(1)

    ok = step3_regrid()
    if not ok:
        print("\nStep 3 FAILED. Check NB12 output for errors.")
        sys.exit(1)

    # Step 4 submits a PBS job (runs on compute node)
    ok = step4_rerun_analysis()
    if not ok:
        print("\nStep 4 FAILED. Submit NB14 manually.")
        sys.exit(1)

    elapsed = (_time.time() - t0) / 60
    print(f"\nSteps 1-3 completed in {elapsed:.0f} minutes.")
    print("Step 4 (NB14 reanalysis) is running as a PBS job.")
    print("Check nb14_33yr.out when the job completes.")
    print("=" * 70)
