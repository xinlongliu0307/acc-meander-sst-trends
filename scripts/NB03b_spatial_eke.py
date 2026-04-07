#!/usr/bin/env python3
"""
NB03b_spatial_eke.py
=====================
Task A: Compute spatially resolved EKE anomaly fields and per-grid-point
EKE trend maps for the meander SST project.

This script extends NB03_speed_eke_trends.py by RETAINING the full
(n_months × nlat × nlon) EKE anomaly cube rather than collapsing it to
a scalar. It then computes Sen's slope and Mann-Kendall significance at
each grid point, producing a spatial trend map.

Outputs per site (saved to meander_sst_project/eke_spatial/):
  - spatial_eke_{SITE}_monthly.nc: monthly EKE field (n_months, nlat, nlon)
  - spatial_eke_{SITE}_trend.nc: per-grid-point trend map (nlat, nlon)

Memory note: Each site's EKE cube is ~390 months × ~90 lat × ~560 lon
× 8 bytes ≈ 150 MB. Fits comfortably in a normal Gadi job (4–8 GB RAM).

Run on NCI Gadi (recommend PBS due to ~2–4 hours per site):
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB03b_spatial_eke.py          # all 4 sites
  python NB03b_spatial_eke.py CP       # single site

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
from netCDF4 import Dataset
import pymannkendall as mk
from statsmodels.tsa.stattools import acf
import warnings, gc, sys, time as _time

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
BASE_DIR    = Path("/g/data/gv90/xl1657/cmems_adt")
ADT_FP      = BASE_DIR / "cmems_so30S_19930101_20250802_adt_sla_ugos_vgos_batch.nc"
OUT_DIR     = BASE_DIR / "meander_sst_project" / "eke_spatial"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SITES = {
    "CP":   {"name": "Campbell Plateau",
             "inner_lon": (150, -150), "inner_lat": (-57, -46), "wraps": True},
    "PAR":  {"name": "Pacific-Antarctic Ridge",
             "inner_lon": (-150, -80), "inner_lat": (-60, -48), "wraps": False},
    "SEIR": {"name": "Southeast Indian Ridge",
             "inner_lon": (130, 152),  "inner_lat": (-56, -44), "wraps": False},
    "SWIR": {"name": "Southwest Indian Ridge",
             "inner_lon": (15, 45),    "inner_lat": (-58, -44), "wraps": False},
}

TIME_CHUNK = 30  # days per read chunk


# =============================================================================
# 1. COMPUTE SPATIAL EKE FIELDS
# =============================================================================
def compute_spatial_eke(site_key):
    """
    For one site, compute the full monthly EKE field at each grid point
    and save as NetCDF. Logic matches NB03 exactly except the spatial
    field is retained instead of collapsed.
    """
    site = SITES[site_key]
    print(f"\n{'='*70}")
    print(f"Computing spatial EKE: {site['name']} ({site_key})")
    print(f"{'='*70}")
    t0 = _time.time()

    with Dataset(str(ADT_FP), "r") as src:
        lat_name = "latitude" if "latitude" in src.variables else "lat"
        lon_name = "longitude" if "longitude" in src.variables else "lon"
        lat = np.array(src.variables[lat_name][:], dtype=np.float64)
        lon = np.array(src.variables[lon_name][:], dtype=np.float64)
        time_raw = src.variables["time"]
        time_all = pd.DatetimeIndex(
            xr.coding.times.decode_cf_datetime(
                time_raw[:], time_raw.units,
                calendar=getattr(time_raw, "calendar", "standard"),
            )
        )
        nt = len(time_all)

        # --- Select site grid points ---
        ilon = site["inner_lon"]
        ilat_range = site["inner_lat"]
        wraps = site.get("wraps", False)
        jlat = np.where((lat >= ilat_range[0]) & (lat <= ilat_range[1]))[0]
        if not wraps:
            jlon = np.where((lon >= ilon[0]) & (lon <= ilon[1]))[0]
        else:
            jlon = np.where((lon >= ilon[0]) | (lon <= ilon[1]))[0]

        lat_s = slice(int(jlat[0]), int(jlat[-1]) + 1)
        nlat_site = int(jlat[-1]) - int(jlat[0]) + 1
        lat_site = lat[jlat[0]:jlat[-1]+1]

        if not wraps or len(jlon) == 0:
            lon_s = slice(int(jlon[0]), int(jlon[-1]) + 1)
            lon_contiguous = True
            nlon_site = int(jlon[-1]) - int(jlon[0]) + 1
            lon_site = lon[jlon[0]:jlon[-1]+1]
        else:
            gap = np.where(np.diff(jlon) > 1)[0]
            if len(gap) > 0:
                split = gap[0] + 1
                lon_s_1 = slice(int(jlon[:split][0]), int(jlon[:split][-1]) + 1)
                lon_s_2 = slice(int(jlon[split:][0]), int(jlon[split:][-1]) + 1)
                lon_contiguous = False
                nlon_site = len(jlon)
                lon_site = lon[jlon]
            else:
                lon_s = slice(int(jlon[0]), int(jlon[-1]) + 1)
                lon_contiguous = True
                nlon_site = int(jlon[-1]) - int(jlon[0]) + 1
                lon_site = lon[jlon[0]:jlon[-1]+1]

        print(f"  Grid: {nlat_site} lat × {nlon_site} lon = {nlat_site * nlon_site} points")

        # --- Build month index ---
        day_months = time_all.to_period("M")
        unique_months = day_months.unique()
        n_months = len(unique_months)
        month_dates = unique_months.to_timestamp() + pd.Timedelta(days=14)
        month_lookup = {m: i for i, m in enumerate(unique_months)}
        day_month_idx = np.array([month_lookup[m] for m in day_months])

        # --- Pass 1: Accumulate monthly-mean u, v at each grid point ---
        print(f"  Pass 1: Accumulating monthly means ({nt} days)...")
        u_month_sum = np.zeros((n_months, nlat_site, nlon_site), dtype=np.float64)
        v_month_sum = np.zeros((n_months, nlat_site, nlon_site), dtype=np.float64)
        month_count = np.zeros(n_months, dtype=np.int32)

        for cs in range(0, nt, TIME_CHUNK):
            ce = min(cs + TIME_CHUNK, nt)
            if lon_contiguous:
                u = np.array(src.variables["ugos"][cs:ce, lat_s, lon_s], dtype=np.float32)
                v = np.array(src.variables["vgos"][cs:ce, lat_s, lon_s], dtype=np.float32)
            else:
                u = np.concatenate([
                    np.array(src.variables["ugos"][cs:ce, lat_s, lon_s_2], dtype=np.float32),
                    np.array(src.variables["ugos"][cs:ce, lat_s, lon_s_1], dtype=np.float32),
                ], axis=2)
                v = np.concatenate([
                    np.array(src.variables["vgos"][cs:ce, lat_s, lon_s_2], dtype=np.float32),
                    np.array(src.variables["vgos"][cs:ce, lat_s, lon_s_1], dtype=np.float32),
                ], axis=2)

            u = np.where(np.abs(u) > 50, np.nan, u)
            v = np.where(np.abs(v) > 50, np.nan, v)

            for k in range(ce - cs):
                mi = day_month_idx[cs + k]
                u_day = np.where(np.isfinite(u[k]), u[k], 0.0)
                v_day = np.where(np.isfinite(v[k]), v[k], 0.0)
                u_month_sum[mi] += u_day
                v_month_sum[mi] += v_day
                month_count[mi] += 1

            if (cs // TIME_CHUNK) % 100 == 0 and cs > 0:
                print(f"    Day {cs}/{nt} ({_time.time()-t0:.0f}s)")

        # Monthly means
        for mi in range(n_months):
            if month_count[mi] > 0:
                u_month_sum[mi] /= month_count[mi]
                v_month_sum[mi] /= month_count[mi]
            else:
                u_month_sum[mi] = np.nan
                v_month_sum[mi] = np.nan

        u_monthly = u_month_sum
        v_monthly = v_month_sum

        # --- Pass 2: Monthly climatology at each grid point ---
        print("  Pass 2: Computing monthly climatology...")
        cal_months = np.array([m.month for m in unique_months])
        u_clim = np.zeros((12, nlat_site, nlon_site), dtype=np.float64)
        v_clim = np.zeros((12, nlat_site, nlon_site), dtype=np.float64)
        clim_count = np.zeros(12, dtype=np.int32)

        for mi in range(n_months):
            cm = cal_months[mi] - 1
            u_clim[cm] += u_monthly[mi]
            v_clim[cm] += v_monthly[mi]
            clim_count[cm] += 1
        for cm in range(12):
            if clim_count[cm] > 0:
                u_clim[cm] /= clim_count[cm]
                v_clim[cm] /= clim_count[cm]

        # --- Pass 3: EKE at each grid point (RETAIN SPATIAL FIELD) ---
        print("  Pass 3: Computing spatial EKE fields...")
        eke_spatial = np.full((n_months, nlat_site, nlon_site), np.nan, dtype=np.float32)

        for mi in range(n_months):
            cm = cal_months[mi] - 1
            u_anom = u_monthly[mi] - u_clim[cm]
            v_anom = v_monthly[mi] - v_clim[cm]
            eke_spatial[mi] = 0.5 * (u_anom**2 + v_anom**2)

    # Free large arrays
    del u_monthly, v_monthly, u_clim, v_clim, u_month_sum, v_month_sum
    gc.collect()

    elapsed = _time.time() - t0
    print(f"  EKE computed in {elapsed/60:.1f} min")
    print(f"  EKE cube shape: {eke_spatial.shape}")
    print(f"  EKE range: {np.nanmin(eke_spatial):.6f}–{np.nanmax(eke_spatial):.6f} m²/s²")

    # --- Save monthly EKE spatial field ---
    ds_eke = xr.Dataset(
        {
            "eke": (["month", "latitude", "longitude"], eke_spatial,
                    {"units": "m2 s-2",
                     "long_name": "Eddy kinetic energy (per grid point)",
                     "description": "EKE = 0.5*(u_anom^2 + v_anom^2), anomalies relative to monthly climatology"}),
        },
        coords={
            "month": month_dates.values,
            "latitude": lat_site,
            "longitude": lon_site,
        },
        attrs={
            "site": site_key,
            "site_name": site["name"],
            "source": str(ADT_FP),
            "method": "Per-grid-point EKE from geostrophic velocity anomalies (CORRECTED pipeline, matching NB03)",
        },
    )
    fp_monthly = OUT_DIR / f"spatial_eke_{site_key}_monthly.nc"
    ds_eke.to_netcdf(fp_monthly, encoding={"eke": {"zlib": True, "complevel": 4}})
    print(f"  Saved monthly EKE: {fp_monthly} ({fp_monthly.stat().st_size/1024**2:.0f} MB)")

    # --- Compute per-grid-point trend ---
    print("  Computing per-grid-point EKE trends (Sen's slope + MK)...")
    t1 = _time.time()

    slope_map = np.full((nlat_site, nlon_site), np.nan, dtype=np.float32)
    pvalue_map = np.full((nlat_site, nlon_site), np.nan, dtype=np.float32)
    sig_map = np.zeros((nlat_site, nlon_site), dtype=np.int8)
    n_computed = 0

    for i in range(nlat_site):
        for j in range(nlon_site):
            ts = eke_spatial[:, i, j]
            valid = np.isfinite(ts)
            if valid.sum() < 24:
                continue

            # Remove seasonal cycle
            ts_clean = ts[valid]
            idx_valid = np.where(valid)[0]
            cal_m = np.array([cal_months[k] for k in idx_valid])
            clim_12 = np.zeros(12)
            clim_n = np.zeros(12)
            for k, cm in enumerate(cal_m):
                clim_12[cm - 1] += ts_clean[k]
                clim_n[cm - 1] += 1
            for cm in range(12):
                if clim_n[cm] > 0:
                    clim_12[cm] /= clim_n[cm]
            ts_anom = ts_clean - np.array([clim_12[cm - 1] for cm in cal_m])

            # Autocorrelation check
            try:
                acf_vals = acf(ts_anom, nlags=5, fft=True)
                lag1 = acf_vals[1]
            except Exception:
                lag1 = 0.0

            # Mann-Kendall test
            try:
                if abs(lag1) > 0.1:
                    result = mk.hamed_rao_modification_test(ts_anom)
                else:
                    result = mk.original_test(ts_anom)
                slope_map[i, j] = result.slope * 12 * 10  # per decade
                pvalue_map[i, j] = result.p
                sig_map[i, j] = 1 if result.p < 0.05 else 0
                n_computed += 1
            except Exception:
                continue

        if (i + 1) % 20 == 0:
            print(f"    Row {i+1}/{nlat_site} ({_time.time()-t1:.0f}s)")

    print(f"  Trends computed at {n_computed} / {nlat_site * nlon_site} grid points")
    print(f"  Significant (p<0.05): {sig_map.sum()} grid points")

    # --- Save trend map ---
    ds_trend = xr.Dataset(
        {
            "eke_trend": (["latitude", "longitude"], slope_map,
                         {"units": "m2 s-2 per decade",
                          "long_name": "EKE Sen slope per decade"}),
            "eke_pvalue": (["latitude", "longitude"], pvalue_map,
                          {"long_name": "Mann-Kendall p-value"}),
            "eke_significant": (["latitude", "longitude"], sig_map,
                               {"long_name": "Significant at p<0.05 (1=yes, 0=no)"}),
        },
        coords={
            "latitude": lat_site,
            "longitude": lon_site,
        },
        attrs={
            "site": site_key,
            "site_name": site["name"],
            "method": "Sen slope + Mann-Kendall (Hamed-Rao modified when |ACF lag-1| > 0.1)",
            "period": f"{month_dates[0].date()} to {month_dates[-1].date()}",
        },
    )
    fp_trend = OUT_DIR / f"spatial_eke_{site_key}_trend.nc"
    ds_trend.to_netcdf(fp_trend)
    print(f"  Saved trend map: {fp_trend}")

    total = _time.time() - t0
    print(f"\n  {site_key} complete in {total/60:.1f} minutes.")

    del eke_spatial
    gc.collect()


# =============================================================================
# 2. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB03b: Spatial EKE Computation for Meander SST Project")
    print(f"  Source: {ADT_FP}")
    print(f"  Output: {OUT_DIR}")
    print("=" * 70)

    if not ADT_FP.exists():
        raise FileNotFoundError(f"ADT file not found: {ADT_FP}")

    # Allow single-site execution via command line
    if len(sys.argv) > 1:
        requested = [s.upper() for s in sys.argv[1:] if s.upper() in SITES]
        if not requested:
            print(f"Unknown site(s). Available: {list(SITES.keys())}")
            sys.exit(1)
    else:
        requested = list(SITES.keys())

    print(f"  Sites to process: {requested}")

    for site_key in requested:
        compute_spatial_eke(site_key)

    print("\n" + "=" * 70)
    print("All sites complete.")
    print("=" * 70)


# =============================================================================
# PBS JOB SCRIPT (copy below into run_spatial_eke.pbs and submit with qsub)
# =============================================================================
PBS_SCRIPT = """
#!/bin/bash
#PBS -N spatial_eke
#PBS -P gv90
#PBS -q normal
#PBS -l walltime=16:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -l storage=gdata/gv90
#PBS -l wd
#PBS -o spatial_eke.out
#PBS -e spatial_eke.err

source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
cd /g/data/gv90/xl1657/cmems_adt/notebooks/
python NB03b_spatial_eke.py
"""

if "--write-pbs" in sys.argv:
    pbs_fp = Path("/g/data/gv90/xl1657/cmems_adt/notebooks/run_spatial_eke.pbs")
    with open(pbs_fp, "w") as f:
        f.write(PBS_SCRIPT.strip() + "\n")
    print(f"PBS script written to: {pbs_fp}")
    print(f"Submit with: qsub {pbs_fp}")
