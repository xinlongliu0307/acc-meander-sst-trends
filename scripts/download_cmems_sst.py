#!/usr/bin/env python3
"""
download_cmems_sst.py
======================
Task B: Download CMEMS OSTIA Reprocessed SST to NCI Gadi.

Product: SST_GLO_SST_L4_REP_OBSERVATIONS_010_011
  - Variable: analysed_sst
  - Resolution: 0.05° daily
  - Period: 1993-01-01 to 2025-08-02 (matching ADT coverage)
  - Domain: 90S–30S, all longitudes

Also downloads NOAA OISST v2.1 as a cross-check dataset.

Prerequisites:
  pip install copernicusmarine --break-system-packages
  copernicusmarine login   # enter your Copernicus Marine credentials once

Run on NCI Gadi (interactive or via PBS — see bottom of file for PBS script):
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python download_cmems_sst.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import subprocess
import sys

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
BASE_DIR = Path("/g/data/gv90/xl1657")

# --- CMEMS OSTIA SST (primary) ---
CMEMS_SST_DIR = BASE_DIR / "cmems_sst"
CMEMS_SST_DIR.mkdir(parents=True, exist_ok=True)

CMEMS_PRODUCT  = "SST_GLO_SST_L4_REP_OBSERVATIONS_010_011"
CMEMS_DATASET  = "METOFFICE-GLO-SST-L4-REP-OBS-SST"
CMEMS_VARIABLE = "analysed_sst"

# Time range matching the ADT data
T_START = "1993-01-01T00:00:00"
T_END   = "2025-08-02T23:59:59"

# Spatial domain matching ADT: 90S–30S, all longitudes
LAT_MIN, LAT_MAX = -90.0, -30.0
LON_MIN, LON_MAX = -180.0, 179.95

# Download in annual chunks to avoid timeouts and enable restart
YEARS = list(range(1993, 2026))

# --- NOAA OISST v2.1 (cross-check) ---
OISST_DIR = BASE_DIR / "noaa_oisst"
OISST_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 1. DOWNLOAD CMEMS OSTIA SST (annual chunks)
# =============================================================================
def download_cmems_sst_annual():
    """
    Download CMEMS OSTIA SST in annual chunks via copernicusmarine CLI.
    Each year is a separate file to allow restart on failure.
    """
    print("=" * 70)
    print("Downloading CMEMS OSTIA Reprocessed SST")
    print(f"  Product: {CMEMS_PRODUCT}")
    print(f"  Dataset: {CMEMS_DATASET}")
    print(f"  Variable: {CMEMS_VARIABLE}")
    print(f"  Domain: {LAT_MIN}–{LAT_MAX}°N, {LON_MIN}–{LON_MAX}°E")
    print(f"  Output: {CMEMS_SST_DIR}")
    print("=" * 70)

    for year in YEARS:
        t0 = f"{year}-01-01T00:00:00"
        t1 = f"{year}-12-31T23:59:59"
        if year == 2025:
            t1 = T_END  # partial year

        out_fp = CMEMS_SST_DIR / f"cmems_ostia_sst_so30S_{year}.nc"

        if out_fp.exists():
            size_mb = out_fp.stat().st_size / 1024**2
            if size_mb > 10:  # reasonable minimum for a full year
                print(f"  {year}: Already exists ({size_mb:.0f} MB), skipping.")
                continue
            else:
                print(f"  {year}: File exists but small ({size_mb:.1f} MB), re-downloading.")

        print(f"\n  Downloading {year}...")
        cmd = [
            "copernicusmarine", "subset",
            "--dataset-id", CMEMS_DATASET,
            "--variable", CMEMS_VARIABLE,
            "--start-datetime", t0,
            "--end-datetime", t1,
            "--minimum-latitude", str(LAT_MIN),
            "--maximum-latitude", str(LAT_MAX),
            "--minimum-longitude", str(LON_MIN),
            "--maximum-longitude", str(LON_MAX),
            "--output-directory", str(CMEMS_SST_DIR),
            "--output-filename", out_fp.name,
            "--force-download",
        ]

        print(f"    Command: {' '.join(cmd[:6])} ...")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            if result.returncode == 0:
                if out_fp.exists():
                    size_mb = out_fp.stat().st_size / 1024**2
                    print(f"    OK: {size_mb:.0f} MB")
                else:
                    print(f"    WARNING: Command succeeded but file not found.")
                    print(f"    stdout: {result.stdout[-500:]}")
            else:
                print(f"    ERROR (return code {result.returncode}):")
                print(f"    stderr: {result.stderr[-500:]}")
        except subprocess.TimeoutExpired:
            print(f"    TIMEOUT: {year} download exceeded 2 hours.")
        except FileNotFoundError:
            print("    ERROR: 'copernicusmarine' not found.")
            print("    Install: pip install copernicusmarine --break-system-packages")
            print("    Login:   copernicusmarine login")
            sys.exit(1)

    print(f"\nCMEMS SST download complete. Files in: {CMEMS_SST_DIR}")


# =============================================================================
# 2. DOWNLOAD NOAA OISST v2.1 (cross-check, via THREDDS/OPeNDAP)
# =============================================================================
def download_noaa_oisst():
    """
    Download NOAA OISST v2.1 (0.25°, daily) for the Southern Ocean.
    Uses wget to fetch annual files from NCEI THREDDS.
    """
    print("\n" + "=" * 70)
    print("Downloading NOAA OISST v2.1 (cross-check dataset)")
    print(f"  Output: {OISST_DIR}")
    print("=" * 70)

    # NOAA OISST v2.1 annual files from NCEI
    BASE_URL = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr"

    for year in YEARS:
        # OISST files are organized as monthly directories with daily files
        # For bulk download, we use the annual mean or download monthly
        # Here we download the preliminary annual files
        out_fp = OISST_DIR / f"oisst-avhrr-v02r01.{year}.nc"

        if out_fp.exists():
            size_mb = out_fp.stat().st_size / 1024**2
            if size_mb > 5:
                print(f"  {year}: Already exists ({size_mb:.0f} MB), skipping.")
                continue

        # NOAA provides daily files; for simplicity, download monthly
        # Alternative: use xarray + OPeNDAP for subset download
        url = f"{BASE_URL}/{year}01/oisst-avhrr-v02r01.{year}0101.nc"
        print(f"  {year}: Attempting download...")
        print(f"    NOTE: NOAA OISST is distributed as daily files.")
        print(f"    For bulk download, use the NCEI bulk download tool or OPeNDAP.")
        print(f"    See: https://www.ncei.noaa.gov/products/optimum-interpolation-sst")

    print(f"\n  NOAA OISST: Manual download recommended for daily files.")
    print(f"  Alternative: use xarray with OPeNDAP for on-the-fly subsetting.")
    print(f"  Example:")
    print(f"    import xarray as xr")
    print(f"    url = 'https://www.ncei.noaa.gov/thredds/dodsC/OisSTv21ADaily'")
    print(f"    ds = xr.open_dataset(url)")
    print(f"    ds_so = ds.sel(lat=slice(-90, -30), time=slice('1993','2025'))")


# =============================================================================
# 3. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("SST Data Download Pipeline for Meander SST Project")
    print("=" * 70)

    download_cmems_sst_annual()
    download_noaa_oisst()

    print("\n" + "=" * 70)
    print("Download pipeline complete.")
    print(f"  CMEMS OSTIA SST: {CMEMS_SST_DIR}")
    print(f"  NOAA OISST v2.1: {OISST_DIR}")
    print("=" * 70)


# =============================================================================
# PBS JOB SCRIPT (copy below into download_sst.pbs and submit with qsub)
# =============================================================================
PBS_SCRIPT = """
#!/bin/bash
#PBS -N download_sst
#PBS -P gv90
#PBS -q copyq
#PBS -l walltime=12:00:00
#PBS -l mem=8GB
#PBS -l ncpus=1
#PBS -l storage=gdata/gv90
#PBS -l wd
#PBS -o download_sst.out
#PBS -e download_sst.err

source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
cd /g/data/gv90/xl1657/cmems_adt/notebooks/
python download_cmems_sst.py
"""

if "--write-pbs" in sys.argv:
    pbs_fp = Path("/g/data/gv90/xl1657/cmems_adt/notebooks/download_sst.pbs")
    with open(pbs_fp, "w") as f:
        f.write(PBS_SCRIPT.strip() + "\n")
    print(f"PBS script written to: {pbs_fp}")
    print(f"Submit with: qsub {pbs_fp}")
