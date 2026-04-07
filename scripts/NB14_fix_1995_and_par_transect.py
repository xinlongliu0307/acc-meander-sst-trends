#!/usr/bin/env python3
"""
NB14_fix_1995_and_par_transect.py
==================================
Fix two issues from Phase 1:

1. Regrid the missing 1995 SST file
2. Recompute the full SST trend map with all 33 years (including 1995)
3. Recompute the along-ACC transect with corrected longitude convention
   for CP (crosses dateline) and PAR (Western Hemisphere)
4. Recompute all regional statistics and tests with the complete dataset

This script replaces the outputs of NB12 (for 1995 only) and NB13 (fully).

Run on NCI Gadi:
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB14_fix_1995_and_par_transect.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import pymannkendall as mk
from statsmodels.tsa.stattools import acf
from scipy import stats as sp_stats
import warnings
import gc
import time as _time

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
BASE_DIR     = Path("/g/data/gv90/xl1657/cmems_adt")
PROJ_DIR     = BASE_DIR / "meander_sst_project"
SST_RAW_DIR  = Path("/g/data/gv90/xl1657/cmems_sst")
SST_REG_DIR  = PROJ_DIR / "sst_regridded"
MASK_DIR     = PROJ_DIR / "masks"
EKE_DIR      = PROJ_DIR / "eke_spatial"
OUT_DIR      = PROJ_DIR / "sst_trends"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_LAT = np.arange(-89.9375, -29.9, 0.125)
TARGET_LON = np.arange(-179.9375, 180.0, 0.125)

SITES = ["CP", "PAR", "SEIR", "SWIR"]
CONTROLS = ["CTRL_SE_PAC", "CTRL_S_IND", "CTRL_S_ATL", "CTRL_CP_SEIR"]
YEARS = list(range(1993, 2026))


# =============================================================================
# 1. FIX 1: REGRID MISSING 1995
# =============================================================================
def regrid_1995():
    """Regrid the 1995 SST file that failed in NB12."""
    year = 1995
    raw_fp = SST_RAW_DIR / f"cmems_ostia_sst_so30S_{year}.nc"
    out_fp = SST_REG_DIR / f"sst_monthly_0125deg_{year}.nc"

    if out_fp.exists():
        size_mb = out_fp.stat().st_size / 1024**2
        if size_mb > 5:
            print(f"  1995: Already regridded ({size_mb:.0f} MB), skipping.")
            return True

    if not raw_fp.exists():
        print(f"  ERROR: Raw 1995 file not found at {raw_fp}")
        return False

    print(f"  Loading {raw_fp.name}...")
    t0 = _time.time()

    try:
        # Load with smaller chunks to reduce peak memory
        ds = xr.open_dataset(raw_fp, chunks={"time": 15})

        # Identify variables
        sst_var = None
        for vname in ds.data_vars:
            if "sst" in vname.lower() or "analysed" in vname.lower():
                sst_var = vname
                break
        if sst_var is None:
            sst_var = list(ds.data_vars)[0]

        lat_name = lon_name = None
        for cname in ds.coords:
            if "lat" in cname.lower():
                lat_name = cname
            if "lon" in cname.lower():
                lon_name = cname

        print(f"    Variable: {sst_var}, shape: {ds[sst_var].shape}")

        # Monthly means
        print(f"    Computing monthly means...")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ds_monthly = ds[sst_var].resample(time="MS").mean().compute()

        n_months = ds_monthly.sizes["time"]
        print(f"    Monthly means: {n_months} months")

        # Subset to Southern Ocean
        src_lat = ds[lat_name].values
        lat_mask = (src_lat >= -90.5) & (src_lat <= -29.5)
        lat_idx = np.where(lat_mask)[0]
        if len(lat_idx) > 0:
            ds_monthly = ds_monthly.isel({lat_name: lat_idx})

        # Regrid
        print(f"    Regridding to 0.125°...")
        ds_regridded = ds_monthly.interp(
            {lat_name: TARGET_LAT, lon_name: TARGET_LON},
            method="nearest"
        ).compute()

        # Kelvin to Celsius
        if np.nanmean(ds_regridded.values) > 200:
            print(f"    Converting K to °C...")
            ds_regridded = ds_regridded - 273.15

        # Save
        time_vals = ds_monthly.time.values
        ds_out = xr.Dataset(
            {"sst": (["time", "latitude", "longitude"],
                     ds_regridded.values.astype(np.float32),
                     {"units": "degC", "long_name": "SST monthly mean regridded"})},
            coords={"time": time_vals, "latitude": TARGET_LAT, "longitude": TARGET_LON},
            attrs={"title": f"CMEMS OSTIA SST regridded, {year}",
                   "source": str(raw_fp)},
        )
        ds_out.to_netcdf(out_fp, encoding={"sst": {"zlib": True, "complevel": 4}})

        size_mb = out_fp.stat().st_size / 1024**2
        elapsed = _time.time() - t0
        print(f"    Saved: {out_fp.name} ({size_mb:.0f} MB, {elapsed:.0f}s)")

        ds.close()
        del ds, ds_monthly, ds_regridded, ds_out
        gc.collect()
        return True

    except Exception as e:
        print(f"    ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False


# =============================================================================
# 2. LOAD ALL REGRIDDED SST
# =============================================================================
def load_sst_timeseries():
    """Load all 33 regridded monthly SST files."""
    print("\nLoading regridded SST data...")
    datasets = []
    for year in YEARS:
        fp = SST_REG_DIR / f"sst_monthly_0125deg_{year}.nc"
        if not fp.exists():
            print(f"  WARNING: {fp.name} not found!")
            continue
        ds = xr.open_dataset(fp)
        datasets.append(ds["sst"])

    sst_all = xr.concat(datasets, dim="time")
    print(f"  Loaded: {sst_all.sizes['time']} months from {len(datasets)} years")

    for ds in datasets:
        ds.close()
    return sst_all


# =============================================================================
# 3. COMPUTE PER-GRID-POINT SST TRENDS (same as NB13)
# =============================================================================
def compute_sst_trends(sst_all):
    """Compute Sen slope + MK at each grid point."""
    print("\nComputing per-grid-point SST trends...")
    t0 = _time.time()

    lat = sst_all.latitude.values
    lon = sst_all.longitude.values
    nlat, nlon = len(lat), len(lon)

    sst_data = sst_all.values
    time_idx = pd.DatetimeIndex(sst_all.time.values)
    cal_months = time_idx.month

    slope_map = np.full((nlat, nlon), np.nan, dtype=np.float32)
    pvalue_map = np.full((nlat, nlon), np.nan, dtype=np.float32)
    sig_map = np.zeros((nlat, nlon), dtype=np.int8)
    mean_sst_map = np.full((nlat, nlon), np.nan, dtype=np.float32)

    n_computed = 0
    n_skipped = 0

    for i in range(nlat):
        for j in range(nlon):
            ts = sst_data[:, i, j]
            valid = np.isfinite(ts)
            if valid.sum() < 60:
                n_skipped += 1
                continue

            mean_sst_map[i, j] = np.nanmean(ts)
            ts_clean = ts[valid]
            cal_m = cal_months[valid]

            clim_12 = np.zeros(12)
            clim_n = np.zeros(12)
            for k in range(len(ts_clean)):
                clim_12[cal_m[k] - 1] += ts_clean[k]
                clim_n[cal_m[k] - 1] += 1
            for cm in range(12):
                if clim_n[cm] > 0:
                    clim_12[cm] /= clim_n[cm]

            ts_anom = ts_clean - np.array([clim_12[cm - 1] for cm in cal_m])

            try:
                acf_vals = acf(ts_anom, nlags=5, fft=True)
                lag1 = acf_vals[1]
            except Exception:
                lag1 = 0.0

            try:
                if abs(lag1) > 0.1:
                    result = mk.hamed_rao_modification_test(ts_anom)
                else:
                    result = mk.original_test(ts_anom)
                slope_map[i, j] = result.slope * 12 * 10
                pvalue_map[i, j] = result.p
                sig_map[i, j] = 1 if result.p < 0.05 else 0
                n_computed += 1
            except Exception:
                n_skipped += 1

        if (i + 1) % 40 == 0:
            print(f"    Row {i+1}/{nlat} ({_time.time()-t0:.0f}s, {n_computed} computed)")

    print(f"  Computed: {n_computed}, Skipped: {n_skipped}")
    print(f"  Significant: {sig_map.sum()}")
    print(f"  Time: {(_time.time()-t0)/60:.1f} min")

    ds_trend = xr.Dataset(
        {
            "sst_trend": (["latitude", "longitude"], slope_map,
                         {"units": "degC per decade"}),
            "sst_pvalue": (["latitude", "longitude"], pvalue_map),
            "sst_significant": (["latitude", "longitude"], sig_map),
            "sst_mean": (["latitude", "longitude"], mean_sst_map,
                        {"units": "degC"}),
        },
        coords={"latitude": lat, "longitude": lon},
        attrs={"title": "Southern Ocean SST trends (1993-2025), all 33 years",
               "method": "Sen slope + MK (Hamed-Rao when |ACF lag-1|>0.1)"},
    )
    fp_out = OUT_DIR / "sst_trend_map.nc"
    ds_trend.to_netcdf(fp_out)
    print(f"  Saved: {fp_out}")

    return slope_map, pvalue_map, sig_map, mean_sst_map, lat, lon


# =============================================================================
# 4. REGIONAL STATISTICS (same as NB13)
# =============================================================================
def extract_regional_stats(slope_map, pvalue_map, sig_map, lat, lon):
    """Extract meander and control regional SST trend statistics."""
    print("\nExtracting regional statistics...")
    rows = []

    for region_key, region_type, mask_prefix in \
        [(s, "meander", "meander_envelope_mask") for s in SITES] + \
        [(c, "control", "control_mask") for c in CONTROLS]:

        fp = MASK_DIR / f"{mask_prefix}_{region_key}.nc"
        if not fp.exists():
            continue

        ds_mask = xr.open_dataset(fp)
        mask_var = "meander_mask" if region_type == "meander" else "control_mask"
        mask = ds_mask[mask_var].values.astype(bool)
        ds_mask.close()

        trends = slope_map[mask]
        valid = np.isfinite(trends)
        if valid.sum() == 0:
            continue

        tv = trends[valid]
        sv = sig_map[mask][valid]

        lat_2d = np.broadcast_to(lat[:, np.newaxis], slope_map.shape)
        wt = np.cos(np.deg2rad(lat_2d))[mask][valid]
        wmean = np.average(tv, weights=wt)

        stats = {
            "region": region_key, "type": region_type,
            "n_gridpoints": int(valid.sum()),
            "mean_trend_degC_dec": float(wmean),
            "median_trend_degC_dec": float(np.median(tv)),
            "std_trend_degC_dec": float(np.std(tv)),
            "frac_significant": float(sv.sum() / valid.sum()),
            "frac_cooling": float((tv < 0).sum() / len(tv)),
            "frac_warming": float((tv > 0).sum() / len(tv)),
        }
        rows.append(stats)
        print(f"  {region_key}: mean={wmean:+.4f} °C/dec, sig={stats['frac_significant']:.1%}")

    df = pd.DataFrame(rows)
    m = df[df["type"] == "meander"]["mean_trend_degC_dec"].mean()
    c = df[df["type"] == "control"]["mean_trend_degC_dec"].mean()
    print(f"\n  MEANDER mean: {m:+.4f}, CONTROL mean: {c:+.4f}, Diff: {m-c:+.4f}")

    df.to_csv(OUT_DIR / "sst_trend_stats.csv", index=False)
    return df


# =============================================================================
# 5. KS TEST (same as NB13)
# =============================================================================
def ks_test(slope_map):
    """KS test: meander vs control."""
    print("\nKS test...")
    meander_t, control_t = [], []

    for sk in SITES:
        fp = MASK_DIR / f"meander_envelope_mask_{sk}.nc"
        if not fp.exists(): continue
        ds = xr.open_dataset(fp)
        m = ds["meander_mask"].values.astype(bool)
        ds.close()
        t = slope_map[m]
        meander_t.append(t[np.isfinite(t)])

    for ck in CONTROLS:
        fp = MASK_DIR / f"control_mask_{ck}.nc"
        if not fp.exists(): continue
        ds = xr.open_dataset(fp)
        m = ds["control_mask"].values.astype(bool)
        ds.close()
        t = slope_map[m]
        control_t.append(t[np.isfinite(t)])

    m_all = np.concatenate(meander_t)
    c_all = np.concatenate(control_t)
    D, p = sp_stats.ks_2samp(m_all, c_all)

    print(f"  Meander: n={len(m_all)}, mean={np.mean(m_all):+.4f}")
    print(f"  Control: n={len(c_all)}, mean={np.mean(c_all):+.4f}")
    print(f"  D={D:.4f}, p={p:.2e}")

    pd.DataFrame([{"D": D, "p": p, "n_meander": len(m_all),
                    "n_control": len(c_all),
                    "mean_meander": np.mean(m_all),
                    "mean_control": np.mean(c_all)}]).to_csv(
        OUT_DIR / "ks_test_meander_vs_control.csv", index=False)


# =============================================================================
# 6. SST-EKE CORRELATION (same as NB13)
# =============================================================================
def sst_eke_corr(slope_map, lat, lon):
    """Spatial correlation between EKE and SST trends per site."""
    print("\nSST-EKE correlation...")
    rows = []

    for sk in SITES:
        mask_fp = MASK_DIR / f"meander_envelope_mask_{sk}.nc"
        eke_fp = EKE_DIR / f"spatial_eke_{sk}_trend.nc"
        if not mask_fp.exists() or not eke_fp.exists():
            continue

        ds_m = xr.open_dataset(mask_fp)
        mask = ds_m["meander_mask"].values.astype(bool)
        ds_m.close()

        ds_e = xr.open_dataset(eke_fp)
        eke_trend = ds_e["eke_trend"].values
        eke_lat = ds_e["latitude"].values
        eke_lon = ds_e["longitude"].values
        ds_e.close()

        eke_full = np.full_like(slope_map, np.nan)
        for i, la in enumerate(eke_lat):
            i_f = np.argmin(np.abs(lat - la))
            for j, lo in enumerate(eke_lon):
                j_f = np.argmin(np.abs(lon - lo))
                eke_full[i_f, j_f] = eke_trend[i, j]

        sst_in = slope_map[mask]
        eke_in = eke_full[mask]
        valid = np.isfinite(sst_in) & np.isfinite(eke_in)
        if valid.sum() < 10:
            continue

        r, p = sp_stats.pearsonr(eke_in[valid], sst_in[valid])
        rows.append({"site": sk, "n": int(valid.sum()), "R": r,
                      "R2": r**2, "p": p})
        print(f"  {sk}: R={r:+.3f}, R²={r**2:.3f}, p={p:.2e}")

    if rows:
        pd.DataFrame(rows).to_csv(OUT_DIR / "sst_eke_correlation.csv", index=False)


# =============================================================================
# 7. FIX 2: ALONG-ACC TRANSECT WITH CORRECTED LONGITUDE HANDLING
# =============================================================================
def compute_acc_transect_fixed(slope_map, lat, lon):
    """
    Along-ACC SST trend transect with corrected longitude convention.

    The circumpolar belt file uses -180 to +180 longitude.
    CP (150 to -150) and PAR (-150 to -80) are in the Western Hemisphere
    or cross the dateline. The fix handles this by matching longitudes
    using minimum angular distance rather than direct value comparison.
    """
    print("\nComputing along-ACC transect (FIXED longitude handling)...")

    env_fp = MASK_DIR / "circumpolar_meander_envelope.nc"
    if not env_fp.exists():
        print("  Envelope not found.")
        return

    ds_env = xr.open_dataset(env_fp)
    center_lat = ds_env["mean_center_lat"].values
    belt_lon = ds_env["longitude"].values
    ds_env.close()

    transect_trend = np.full(len(belt_lon), np.nan)

    for j, blon in enumerate(belt_lon):
        clat = center_lat[j]
        if np.isnan(clat):
            continue

        # Match longitude using minimum angular distance
        # This correctly handles the -180/+180 boundary
        lon_diff = np.abs(lon - blon)
        lon_diff = np.minimum(lon_diff, 360.0 - lon_diff)
        j_grid = np.argmin(lon_diff)

        i_grid = np.argmin(np.abs(lat - clat))

        if 0 <= i_grid < len(lat) and 0 <= j_grid < len(lon):
            transect_trend[j] = slope_map[i_grid, j_grid]

    # Save transect
    ds_t = xr.Dataset(
        {
            "sst_trend_along_acc": (["longitude"], transect_trend.astype(np.float32),
                                    {"units": "degC per decade"}),
            "acc_center_lat": (["longitude"], center_lat.astype(np.float32),
                              {"units": "degrees_north"}),
        },
        coords={"longitude": belt_lon},
        attrs={"title": "Along-ACC SST trend transect (corrected lon handling)"},
    )
    ds_t.to_netcdf(OUT_DIR / "sst_trend_along_acc_transect.nc")
    print(f"  Saved transect.")

    # Site-specific anomalies using the -180/+180 convention
    # CP: 150E to -150E (i.e., 150 to 210 in 0-360, but in -180/180: 150 to -150)
    # PAR: -150 to -80
    # SEIR: 130 to 152
    # SWIR: 15 to 45
    site_lon_ranges = {
        "CP":   {"lon_min": 150, "lon_max": -150, "wraps": True},
        "PAR":  {"lon_min": -150, "lon_max": -80, "wraps": False},
        "SEIR": {"lon_min": 130, "lon_max": 152, "wraps": False},
        "SWIR": {"lon_min": 15, "lon_max": 45, "wraps": False},
    }

    circ_mean = np.nanmean(transect_trend)
    print(f"  Circumpolar mean: {circ_mean:+.4f} °C/dec")
    print(f"  Valid transect points: {np.isfinite(transect_trend).sum()} / {len(transect_trend)}")

    for sk, info in site_lon_ranges.items():
        lo_min, lo_max, wraps = info["lon_min"], info["lon_max"], info["wraps"]
        if not wraps:
            mask = (belt_lon >= lo_min) & (belt_lon <= lo_max)
        else:
            mask = (belt_lon >= lo_min) | (belt_lon <= lo_max)

        vals = transect_trend[mask]
        valid_vals = vals[np.isfinite(vals)]

        if len(valid_vals) > 0:
            site_mean = np.nanmean(valid_vals)
            anomaly = site_mean - circ_mean
            print(f"  {sk}: mean={site_mean:+.4f} °C/dec, "
                  f"anomaly={anomaly:+.4f} °C/dec, "
                  f"n_valid={len(valid_vals)}")
        else:
            print(f"  {sk}: NO VALID DATA (mask selects {mask.sum()} points)")

    # Save site anomalies
    rows = []
    for sk, info in site_lon_ranges.items():
        lo_min, lo_max, wraps = info["lon_min"], info["lon_max"], info["wraps"]
        if not wraps:
            mask = (belt_lon >= lo_min) & (belt_lon <= lo_max)
        else:
            mask = (belt_lon >= lo_min) | (belt_lon <= lo_max)
        vals = transect_trend[mask]
        valid_vals = vals[np.isfinite(vals)]
        if len(valid_vals) > 0:
            rows.append({
                "site": sk,
                "mean_sst_trend": float(np.nanmean(valid_vals)),
                "anomaly_from_circumpolar": float(np.nanmean(valid_vals) - circ_mean),
                "n_valid": len(valid_vals),
                "circumpolar_mean": float(circ_mean),
            })

    if rows:
        pd.DataFrame(rows).to_csv(OUT_DIR / "along_acc_site_anomalies.csv", index=False)
        print(f"  Saved site anomalies CSV.")


# =============================================================================
# 8. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB14: Fix 1995 + PAR Transect + Recompute All Phase 1 Results")
    print("=" * 70)

    # Fix 1: Regrid 1995
    print("\n--- Step 1: Regrid 1995 ---")
    success = regrid_1995()
    if not success:
        print("  WARNING: 1995 regridding failed. Proceeding with 32 years.")

    # Verify all years present
    n_files = len(list(SST_REG_DIR.glob("sst_monthly_0125deg_*.nc")))
    print(f"  Regridded files available: {n_files}")

    # Load all SST
    print("\n--- Step 2: Load SST ---")
    sst_all = load_sst_timeseries()

    # Recompute trends with all years
    print("\n--- Step 3: Compute SST trends ---")
    slope_map, pvalue_map, sig_map, mean_sst_map, lat, lon = \
        compute_sst_trends(sst_all)
    del sst_all
    gc.collect()

    # Regional stats
    print("\n--- Step 4: Regional statistics ---")
    extract_regional_stats(slope_map, pvalue_map, sig_map, lat, lon)

    # KS test
    print("\n--- Step 5: KS test ---")
    ks_test(slope_map)

    # SST-EKE correlation
    print("\n--- Step 6: SST-EKE correlation ---")
    sst_eke_corr(slope_map, lat, lon)

    # Fixed transect
    print("\n--- Step 7: Along-ACC transect (FIXED) ---")
    compute_acc_transect_fixed(slope_map, lat, lon)

    print("\n" + "=" * 70)
    print("NB14 complete. All Phase 1 outputs updated.")
    print(f"Output: {OUT_DIR}")
    print("=" * 70)
