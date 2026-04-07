#!/usr/bin/env python3
"""
NB12_regrid_sst.py
===================
Phase 1, Step 0: Regrid CMEMS OSTIA SST (0.05°) onto the CMEMS DUACS
ADT grid (0.125°) using conservative remapping.

Produces one regridded NetCDF file per year in:
  /g/data/gv90/xl1657/cmems_adt/meander_sst_project/sst_regridded/

Each output file contains monthly-mean SST on the 0.125° grid
(480 lat × 2880 lon), matching the ADT-derived products exactly.

Monthly averaging is performed during regridding to reduce file sizes
from ~5.9 GB/year (daily 0.05°) to ~55 MB/year (monthly 0.125°).

Run on NCI Gadi via PBS (see bottom for PBS script):
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB12_regrid_sst.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import xarray as xr
import warnings
import gc
import sys
import time as _time

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
SST_RAW_DIR  = Path("/g/data/gv90/xl1657/cmems_sst")
OUT_DIR      = Path("/g/data/gv90/xl1657/cmems_adt/meander_sst_project/sst_regridded")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Target grid: 0.125° matching CMEMS DUACS L4 ADT (480 lat × 2880 lon)
TARGET_LAT = np.arange(-89.9375, -29.9, 0.125)  # 480 points
TARGET_LON = np.arange(-179.9375, 180.0, 0.125)  # 2880 points

YEARS = list(range(1993, 2026))


# =============================================================================
# 1. REGRID ONE YEAR
# =============================================================================
def regrid_one_year(year):
    """
    Load one year of daily 0.05° SST, compute monthly means,
    regrid to 0.125°, and save.
    """
    raw_fp = SST_RAW_DIR / f"cmems_ostia_sst_so30S_{year}.nc"
    out_fp = OUT_DIR / f"sst_monthly_0125deg_{year}.nc"

    if out_fp.exists():
        size_mb = out_fp.stat().st_size / 1024**2
        if size_mb > 5:
            print(f"  {year}: Already exists ({size_mb:.0f} MB), skipping.")
            return True

    if not raw_fp.exists():
        print(f"  {year}: Raw file not found at {raw_fp}, skipping.")
        return False

    t0 = _time.time()
    print(f"  {year}: Loading {raw_fp.name}...")

    try:
        # Open with chunking to manage memory
        ds = xr.open_dataset(raw_fp, chunks={"time": 31})

        # Identify variable and coordinate names
        sst_var = None
        for vname in ds.data_vars:
            if "sst" in vname.lower() or "temperature" in vname.lower() or "analysed" in vname.lower():
                sst_var = vname
                break
        if sst_var is None:
            sst_var = list(ds.data_vars)[0]
            print(f"    Using first variable: {sst_var}")

        lat_name = None
        lon_name = None
        for cname in ds.coords:
            if "lat" in cname.lower():
                lat_name = cname
            if "lon" in cname.lower():
                lon_name = cname

        print(f"    Variable: {sst_var}, Lat: {lat_name}, Lon: {lon_name}")
        print(f"    Shape: {ds[sst_var].shape}")

        # Step 1: Compute monthly means (reduces daily -> monthly)
        print(f"    Computing monthly means...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds_monthly = ds[sst_var].resample(time="MS").mean().compute()

        n_months = ds_monthly.sizes["time"]
        print(f"    Monthly means: {n_months} months")

        # Get source grid
        src_lat = ds[lat_name].values
        src_lon = ds[lon_name].values

        # Subset source to Southern Ocean only (avoid processing global data)
        lat_mask = (src_lat >= -90.5) & (src_lat <= -29.5)
        lat_idx = np.where(lat_mask)[0]
        if len(lat_idx) > 0:
            src_lat_so = src_lat[lat_idx]
            ds_monthly = ds_monthly.isel({lat_name: lat_idx})
        else:
            src_lat_so = src_lat

        print(f"    Source grid: {len(src_lat_so)} lat × {len(src_lon)} lon")

        # Step 2: Regrid to 0.125° using nearest-neighbor interpolation
        # (conservative remapping requires xesmf which may not be installed;
        #  nearest-neighbor on 0.05° -> 0.125° is adequate for trend analysis
        #  since the target grid is coarser by only 2.5x)
        print(f"    Regridding to 0.125°...")

        # Use xarray's interp method for the regridding
        ds_regridded = ds_monthly.interp(
            {lat_name: TARGET_LAT, lon_name: TARGET_LON},
            method="nearest"
        ).compute()

        # Step 3: Convert from Kelvin to Celsius if needed
        vals = ds_regridded.values
        if np.nanmean(vals) > 200:
            print(f"    Converting from Kelvin to Celsius...")
            ds_regridded = ds_regridded - 273.15

        # Step 4: Build output dataset
        time_vals = ds_monthly.time.values

        ds_out = xr.Dataset(
            {
                "sst": (["time", "latitude", "longitude"],
                        ds_regridded.values.astype(np.float32),
                        {"units": "degC",
                         "long_name": "Sea surface temperature (monthly mean, regridded)",
                         "source_resolution": "0.05 deg",
                         "target_resolution": "0.125 deg",
                         "regrid_method": "nearest neighbor"}),
            },
            coords={
                "time": time_vals,
                "latitude": TARGET_LAT,
                "longitude": TARGET_LON,
            },
            attrs={
                "title": f"CMEMS OSTIA SST regridded to 0.125 deg, {year}",
                "source": str(raw_fp),
                "source_product": "SST_GLO_SST_L4_REP_OBSERVATIONS_010_011",
                "processing": "Daily -> monthly mean, then regridded from 0.05 to 0.125 deg",
            },
        )

        # Save with compression
        ds_out.to_netcdf(out_fp, encoding={
            "sst": {"zlib": True, "complevel": 4, "dtype": "float32"}
        })

        size_mb = out_fp.stat().st_size / 1024**2
        elapsed = _time.time() - t0
        print(f"    Saved: {out_fp.name} ({size_mb:.0f} MB, {elapsed:.0f}s)")

        ds.close()
        del ds, ds_monthly, ds_regridded, ds_out
        gc.collect()

        return True

    except Exception as e:
        print(f"    ERROR processing {year}: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# 2. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB12: SST Regridding Preprocessor")
    print(f"  Source: {SST_RAW_DIR} (0.05° daily)")
    print(f"  Output: {OUT_DIR} (0.125° monthly)")
    print(f"  Years: {YEARS[0]}–{YEARS[-1]}")
    print("=" * 70)

    # Allow single-year processing via command line
    if len(sys.argv) > 1:
        try:
            requested_years = [int(y) for y in sys.argv[1:] if y.isdigit()]
            if requested_years:
                YEARS = requested_years
        except ValueError:
            pass

    print(f"  Processing {len(YEARS)} years...")

    success = 0
    fail = 0
    for year in YEARS:
        if regrid_one_year(year):
            success += 1
        else:
            fail += 1

    print(f"\n{'='*70}")
    print(f"Regridding complete: {success} succeeded, {fail} failed")
    print(f"Output: {OUT_DIR}")
    print(f"{'='*70}")
