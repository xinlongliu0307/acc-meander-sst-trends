#!/usr/bin/env python3
"""
NB13_sst_trends.py
===================
Phase 1: Compute per-grid-point SST trends (1993–2025) and extract
meander-vs-control regional statistics for the GRL manuscript.

Inputs:
  - Regridded monthly SST files from NB12 (0.125°, in sst_regridded/)
  - Meander envelope masks (from setup_meander_sst_project.py)
  - Control region masks (from setup_meander_sst_project.py)
  - Spatial EKE trend maps (from NB03b_spatial_eke.py)

Outputs (saved to meander_sst_project/sst_trends/):
  - sst_trend_map.nc: per-grid-point SST trend, p-value, significance
  - sst_trend_stats.csv: regional statistics for manuscript Table/text
  - sst_eke_correlation.csv: EKE-SST spatial correlation per site

Run on NCI Gadi via PBS:
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB13_sst_trends.py

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
BASE_DIR    = Path("/g/data/gv90/xl1657/cmems_adt")
PROJ_DIR    = BASE_DIR / "meander_sst_project"
SST_DIR     = PROJ_DIR / "sst_regridded"
MASK_DIR    = PROJ_DIR / "masks"
EKE_DIR     = PROJ_DIR / "eke_spatial"
OUT_DIR     = PROJ_DIR / "sst_trends"
OUT_DIR.mkdir(parents=True, exist_ok=True)

YEARS = list(range(1993, 2026))

SITES = ["CP", "PAR", "SEIR", "SWIR"]
CONTROLS = ["CTRL_SE_PAC", "CTRL_S_IND", "CTRL_S_ATL", "CTRL_CP_SEIR"]


# =============================================================================
# 1. LOAD ALL REGRIDDED SST INTO A SINGLE TIME SERIES
# =============================================================================
def load_sst_timeseries():
    """Load all regridded monthly SST files and concatenate along time."""
    print("Loading regridded SST data...")
    datasets = []
    for year in YEARS:
        fp = SST_DIR / f"sst_monthly_0125deg_{year}.nc"
        if not fp.exists():
            print(f"  WARNING: {fp.name} not found, skipping.")
            continue
        ds = xr.open_dataset(fp)
        datasets.append(ds["sst"])
        print(f"  {year}: {ds.sizes.get('time', '?')} months")

    if len(datasets) == 0:
        raise FileNotFoundError(f"No regridded SST files found in {SST_DIR}")

    sst_all = xr.concat(datasets, dim="time")
    print(f"  Total: {sst_all.sizes['time']} months, "
          f"{sst_all.sizes['latitude']} lat × {sst_all.sizes['longitude']} lon")

    # Close individual datasets
    for ds in datasets:
        ds.close()

    return sst_all


# =============================================================================
# 2. COMPUTE PER-GRID-POINT SST TRENDS
# =============================================================================
def compute_sst_trends(sst_all):
    """
    Compute Sen's slope + Mann-Kendall trend at each grid point.
    Uses the same statistical framework as NB03/NB03b for consistency.
    """
    print("\nComputing per-grid-point SST trends...")
    t0 = _time.time()

    lat = sst_all.latitude.values
    lon = sst_all.longitude.values
    nlat = len(lat)
    nlon = len(lon)
    nt = sst_all.sizes["time"]

    # Load into memory as numpy array
    print(f"  Loading SST cube into memory ({nlat}×{nlon}×{nt})...")
    sst_data = sst_all.values  # (time, lat, lon)

    # Pre-compute calendar months for seasonal detrending
    time_idx = pd.DatetimeIndex(sst_all.time.values)
    cal_months = time_idx.month

    # Output arrays
    slope_map = np.full((nlat, nlon), np.nan, dtype=np.float32)
    pvalue_map = np.full((nlat, nlon), np.nan, dtype=np.float32)
    sig_map = np.zeros((nlat, nlon), dtype=np.int8)
    mean_sst_map = np.full((nlat, nlon), np.nan, dtype=np.float32)

    n_computed = 0
    n_skipped = 0

    print(f"  Processing {nlat} latitude rows...")

    for i in range(nlat):
        for j in range(nlon):
            ts = sst_data[:, i, j]
            valid = np.isfinite(ts)

            if valid.sum() < 60:  # require at least 5 years of data
                n_skipped += 1
                continue

            # Mean SST for reference
            mean_sst_map[i, j] = np.nanmean(ts)

            # Remove seasonal cycle
            ts_clean = ts[valid]
            cal_m = cal_months[valid]

            clim_12 = np.zeros(12)
            clim_n = np.zeros(12)
            for k in range(len(ts_clean)):
                cm = cal_m[k] - 1
                clim_12[cm] += ts_clean[k]
                clim_n[cm] += 1
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
                n_skipped += 1
                continue

        if (i + 1) % 20 == 0:
            elapsed = _time.time() - t0
            print(f"    Row {i+1}/{nlat} ({elapsed:.0f}s, "
                  f"{n_computed} computed, {n_skipped} skipped)")

    print(f"  Trends computed at {n_computed} grid points, {n_skipped} skipped")
    print(f"  Significant (p<0.05): {sig_map.sum()} grid points")
    print(f"  Total time: {(_time.time()-t0)/60:.1f} minutes")

    # Save trend map
    ds_trend = xr.Dataset(
        {
            "sst_trend": (["latitude", "longitude"], slope_map,
                         {"units": "degC per decade",
                          "long_name": "SST Sen slope per decade"}),
            "sst_pvalue": (["latitude", "longitude"], pvalue_map,
                          {"long_name": "Mann-Kendall p-value"}),
            "sst_significant": (["latitude", "longitude"], sig_map,
                               {"long_name": "Significant at p<0.05"}),
            "sst_mean": (["latitude", "longitude"], mean_sst_map,
                        {"units": "degC",
                         "long_name": "Time-mean SST (1993-2025)"}),
        },
        coords={
            "latitude": lat,
            "longitude": lon,
        },
        attrs={
            "title": "Southern Ocean SST trends (1993-2025)",
            "source": "CMEMS OSTIA Reprocessed, regridded to 0.125 deg",
            "method": "Sen slope + Mann-Kendall (Hamed-Rao when |ACF lag-1|>0.1)",
            "period": f"{YEARS[0]}-01 to {YEARS[-1]}-12",
        },
    )

    fp_out = OUT_DIR / "sst_trend_map.nc"
    ds_trend.to_netcdf(fp_out)
    print(f"  Saved: {fp_out}")

    return slope_map, pvalue_map, sig_map, mean_sst_map, lat, lon


# =============================================================================
# 3. EXTRACT REGIONAL STATISTICS
# =============================================================================
def extract_regional_stats(slope_map, pvalue_map, sig_map, lat, lon):
    """
    Extract area-weighted SST trend statistics for each meander site
    and each control region. Produces the core numbers for the manuscript.
    """
    print("\nExtracting regional SST trend statistics...")

    rows = []

    # Meander sites
    for site_key in SITES:
        fp = MASK_DIR / f"meander_envelope_mask_{site_key}.nc"
        if not fp.exists():
            print(f"  {site_key}: mask not found, skipping.")
            continue

        ds_mask = xr.open_dataset(fp)
        mask = ds_mask["meander_mask"].values.astype(bool)
        ds_mask.close()

        trends_in_mask = slope_map[mask]
        pvals_in_mask = pvalue_map[mask]
        sig_in_mask = sig_map[mask]

        valid = np.isfinite(trends_in_mask)
        if valid.sum() == 0:
            continue

        trends_valid = trends_in_mask[valid]

        # Compute area weights (cos latitude)
        lat_2d = np.broadcast_to(lat[:, np.newaxis], slope_map.shape)
        cos_weights = np.cos(np.deg2rad(lat_2d))
        weights_in_mask = cos_weights[mask][valid]
        weighted_mean = np.average(trends_valid, weights=weights_in_mask)

        stats = {
            "region": site_key,
            "type": "meander",
            "n_gridpoints": int(valid.sum()),
            "mean_trend_degC_dec": float(weighted_mean),
            "median_trend_degC_dec": float(np.median(trends_valid)),
            "std_trend_degC_dec": float(np.std(trends_valid)),
            "frac_significant": float(sig_in_mask[valid].sum() / valid.sum()),
            "frac_cooling": float((trends_valid < 0).sum() / len(trends_valid)),
            "frac_warming": float((trends_valid > 0).sum() / len(trends_valid)),
            "p10_trend": float(np.percentile(trends_valid, 10)),
            "p90_trend": float(np.percentile(trends_valid, 90)),
        }
        rows.append(stats)
        print(f"  {site_key}: mean={weighted_mean:+.4f} °C/dec, "
              f"sig={stats['frac_significant']:.1%}, "
              f"n={stats['n_gridpoints']}")

    # Control regions
    for ctrl_key in CONTROLS:
        fp = MASK_DIR / f"control_mask_{ctrl_key}.nc"
        if not fp.exists():
            print(f"  {ctrl_key}: mask not found, skipping.")
            continue

        ds_mask = xr.open_dataset(fp)
        mask = ds_mask["control_mask"].values.astype(bool)
        ds_mask.close()

        trends_in_mask = slope_map[mask]
        valid = np.isfinite(trends_in_mask)
        if valid.sum() == 0:
            continue

        trends_valid = trends_in_mask[valid]
        sig_in_mask = sig_map[mask][valid]

        lat_2d = np.broadcast_to(lat[:, np.newaxis], slope_map.shape)
        cos_weights = np.cos(np.deg2rad(lat_2d))
        weights_in_mask = cos_weights[mask][valid]
        weighted_mean = np.average(trends_valid, weights=weights_in_mask)

        stats = {
            "region": ctrl_key,
            "type": "control",
            "n_gridpoints": int(valid.sum()),
            "mean_trend_degC_dec": float(weighted_mean),
            "median_trend_degC_dec": float(np.median(trends_valid)),
            "std_trend_degC_dec": float(np.std(trends_valid)),
            "frac_significant": float(sig_in_mask.sum() / valid.sum()),
            "frac_cooling": float((trends_valid < 0).sum() / len(trends_valid)),
            "frac_warming": float((trends_valid > 0).sum() / len(trends_valid)),
            "p10_trend": float(np.percentile(trends_valid, 10)),
            "p90_trend": float(np.percentile(trends_valid, 90)),
        }
        rows.append(stats)
        print(f"  {ctrl_key}: mean={weighted_mean:+.4f} °C/dec, "
              f"sig={stats['frac_significant']:.1%}, "
              f"n={stats['n_gridpoints']}")

    df = pd.DataFrame(rows)

    # Compute meander vs control aggregate statistics
    meander_rows = df[df["type"] == "meander"]
    control_rows = df[df["type"] == "control"]

    if len(meander_rows) > 0 and len(control_rows) > 0:
        m_mean = meander_rows["mean_trend_degC_dec"].mean()
        c_mean = control_rows["mean_trend_degC_dec"].mean()
        diff = m_mean - c_mean
        print(f"\n  MEANDER aggregate mean: {m_mean:+.4f} °C/dec")
        print(f"  CONTROL aggregate mean: {c_mean:+.4f} °C/dec")
        print(f"  Difference (meander - control): {diff:+.4f} °C/dec")

    fp_out = OUT_DIR / "sst_trend_stats.csv"
    df.to_csv(fp_out, index=False)
    print(f"\n  Saved: {fp_out}")

    return df


# =============================================================================
# 4. MEANDER VS CONTROL STATISTICAL TEST
# =============================================================================
def meander_vs_control_test(slope_map, lat):
    """
    Kolmogorov-Smirnov test comparing meander and control SST trend
    distributions. Produces the p-value cited in the manuscript.
    """
    print("\nMeander vs Control: Kolmogorov-Smirnov test...")

    # Collect all meander grid-point trends
    meander_trends = []
    for site_key in SITES:
        fp = MASK_DIR / f"meander_envelope_mask_{site_key}.nc"
        if not fp.exists():
            continue
        ds_mask = xr.open_dataset(fp)
        mask = ds_mask["meander_mask"].values.astype(bool)
        ds_mask.close()
        t = slope_map[mask]
        meander_trends.append(t[np.isfinite(t)])

    # Collect all control grid-point trends
    control_trends = []
    for ctrl_key in CONTROLS:
        fp = MASK_DIR / f"control_mask_{ctrl_key}.nc"
        if not fp.exists():
            continue
        ds_mask = xr.open_dataset(fp)
        mask = ds_mask["control_mask"].values.astype(bool)
        ds_mask.close()
        t = slope_map[mask]
        control_trends.append(t[np.isfinite(t)])

    if len(meander_trends) == 0 or len(control_trends) == 0:
        print("  Insufficient data for KS test.")
        return

    m_all = np.concatenate(meander_trends)
    c_all = np.concatenate(control_trends)

    ks_stat, ks_p = sp_stats.ks_2samp(m_all, c_all)

    print(f"  Meander: n={len(m_all)}, mean={np.mean(m_all):+.4f} °C/dec")
    print(f"  Control: n={len(c_all)}, mean={np.mean(c_all):+.4f} °C/dec")
    print(f"  KS statistic: D = {ks_stat:.4f}")
    print(f"  KS p-value: p = {ks_p:.2e}")
    print(f"  Significant at p<0.05: {'YES' if ks_p < 0.05 else 'NO'}")

    # Save result
    ks_df = pd.DataFrame([{
        "test": "KS_2sample",
        "n_meander": len(m_all),
        "n_control": len(c_all),
        "mean_meander": np.mean(m_all),
        "mean_control": np.mean(c_all),
        "difference": np.mean(m_all) - np.mean(c_all),
        "ks_statistic": ks_stat,
        "ks_pvalue": ks_p,
        "significant_005": ks_p < 0.05,
    }])
    fp_out = OUT_DIR / "ks_test_meander_vs_control.csv"
    ks_df.to_csv(fp_out, index=False)
    print(f"  Saved: {fp_out}")


# =============================================================================
# 5. SST-EKE SPATIAL CORRELATION
# =============================================================================
def compute_sst_eke_correlation(slope_map, lat, lon):
    """
    For each meander site, compute the spatial correlation between
    per-grid-point EKE trends and SST trends within the meander envelope.
    """
    print("\nComputing SST-EKE spatial correlation per site...")

    rows = []
    for site_key in SITES:
        mask_fp = MASK_DIR / f"meander_envelope_mask_{site_key}.nc"
        eke_fp = EKE_DIR / f"spatial_eke_{site_key}_trend.nc"

        if not mask_fp.exists() or not eke_fp.exists():
            print(f"  {site_key}: missing files, skipping.")
            continue

        ds_mask = xr.open_dataset(mask_fp)
        mask_full = ds_mask["meander_mask"].values.astype(bool)
        ds_mask.close()

        ds_eke = xr.open_dataset(eke_fp)
        eke_trend = ds_eke["eke_trend"].values  # (nlat_site, nlon_site)
        eke_lat = ds_eke["latitude"].values
        eke_lon = ds_eke["longitude"].values
        ds_eke.close()

        # Map EKE trend onto the full grid
        eke_full = np.full_like(slope_map, np.nan)
        for i, la in enumerate(eke_lat):
            for j, lo in enumerate(eke_lon):
                i_full = np.argmin(np.abs(lat - la))
                j_full = np.argmin(np.abs(lon - lo))
                eke_full[i_full, j_full] = eke_trend[i, j]

        # Extract co-located values within meander envelope
        sst_in = slope_map[mask_full]
        eke_in = eke_full[mask_full]

        valid = np.isfinite(sst_in) & np.isfinite(eke_in)
        if valid.sum() < 10:
            print(f"  {site_key}: too few co-located points ({valid.sum()})")
            continue

        sst_v = sst_in[valid]
        eke_v = eke_in[valid]

        r, p = sp_stats.pearsonr(eke_v, sst_v)
        r2 = r**2

        rows.append({
            "site": site_key,
            "n_points": int(valid.sum()),
            "pearson_r": r,
            "r_squared": r2,
            "p_value": p,
            "significant": p < 0.05,
        })

        print(f"  {site_key}: R={r:+.3f}, R²={r2:.3f}, p={p:.2e}, n={valid.sum()}")

    if rows:
        df = pd.DataFrame(rows)
        fp_out = OUT_DIR / "sst_eke_correlation.csv"
        df.to_csv(fp_out, index=False)
        print(f"  Saved: {fp_out}")


# =============================================================================
# 6. ALONG-ACC SST TREND TRANSECT
# =============================================================================
def compute_acc_transect(slope_map, lat, lon):
    """
    Extract SST trend along the time-mean ACC axis (from the circumpolar
    belt file) at each longitude.
    """
    print("\nComputing along-ACC SST trend transect...")

    env_fp = MASK_DIR / "circumpolar_meander_envelope.nc"
    if not env_fp.exists():
        print("  Envelope file not found, skipping transect.")
        return

    ds_env = xr.open_dataset(env_fp)
    center_lat = ds_env["mean_center_lat"].values  # (nlon_belt,)
    belt_lon = ds_env["longitude"].values
    ds_env.close()

    # For each longitude in the belt, extract the SST trend at the center lat
    transect_trend = np.full(len(belt_lon), np.nan)
    transect_lon = belt_lon.copy()

    for j, blon in enumerate(belt_lon):
        clat = center_lat[j]
        if np.isnan(clat):
            continue

        # Find nearest grid point
        j_grid = np.argmin(np.abs(lon - blon))
        i_grid = np.argmin(np.abs(lat - clat))

        if 0 <= i_grid < len(lat) and 0 <= j_grid < len(lon):
            transect_trend[j] = slope_map[i_grid, j_grid]

    # Save transect
    ds_transect = xr.Dataset(
        {
            "sst_trend_along_acc": (["longitude"], transect_trend.astype(np.float32),
                                    {"units": "degC per decade",
                                     "long_name": "SST trend along ACC center"}),
            "acc_center_lat": (["longitude"], center_lat.astype(np.float32),
                              {"units": "degrees_north",
                               "long_name": "Mean ACC center latitude"}),
        },
        coords={"longitude": transect_lon},
        attrs={"title": "Along-ACC SST trend transect"},
    )
    fp_out = OUT_DIR / "sst_trend_along_acc_transect.nc"
    ds_transect.to_netcdf(fp_out)
    print(f"  Saved: {fp_out}")

    # Print summary statistics at meander sites
    site_lons = {
        "SWIR": (15, 45), "SEIR": (130, 152),
        "CP": (150, 210), "PAR": (210, 280),
    }
    circumpolar_mean = np.nanmean(transect_trend)
    print(f"  Circumpolar mean SST trend: {circumpolar_mean:+.4f} °C/dec")

    for sk, (lo0, lo1) in site_lons.items():
        if lo0 < lo1:
            mask = (transect_lon >= lo0) & (transect_lon <= lo1)
        else:
            mask = (transect_lon >= lo0) | (transect_lon <= lo1)
        site_mean = np.nanmean(transect_trend[mask])
        anomaly = site_mean - circumpolar_mean
        print(f"  {sk}: mean={site_mean:+.4f} °C/dec, "
              f"anomaly from circumpolar={anomaly:+.4f} °C/dec")


# =============================================================================
# 7. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB13: Phase 1 SST Trend Analysis")
    print(f"  SST source: {SST_DIR}")
    print(f"  Masks: {MASK_DIR}")
    print(f"  EKE trends: {EKE_DIR}")
    print(f"  Output: {OUT_DIR}")
    print("=" * 70)

    # Step 1: Load SST
    sst_all = load_sst_timeseries()

    # Step 2: Compute trends
    slope_map, pvalue_map, sig_map, mean_sst_map, lat, lon = \
        compute_sst_trends(sst_all)

    del sst_all
    gc.collect()

    # Step 3: Regional statistics
    extract_regional_stats(slope_map, pvalue_map, sig_map, lat, lon)

    # Step 4: KS test
    meander_vs_control_test(slope_map, lat)

    # Step 5: SST-EKE correlation
    compute_sst_eke_correlation(slope_map, lat, lon)

    # Step 6: Along-ACC transect
    compute_acc_transect(slope_map, lat, lon)

    print("\n" + "=" * 70)
    print("Phase 1 SST trend analysis complete.")
    print(f"All outputs in: {OUT_DIR}")
    print("=" * 70)
