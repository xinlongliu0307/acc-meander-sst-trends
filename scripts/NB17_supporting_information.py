#!/usr/bin/env python3
"""
NB17_supporting_information.py
===============================
Produce Supporting Information outputs for the GRL manuscript:

  1. Table S1: Control region definitions (CSV and LaTeX)
  2. Figure S1: Decadal decomposition of SST trends (three sub-periods)
  3. Table S2: Per-site and per-control SST trend statistics (full detail)

Run on NCI Gadi via PBS:
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB17_supporting_information.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import pymannkendall as mk
from statsmodels.tsa.stattools import acf
from scipy import stats as sp_stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
import gc
import time as _time

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR   = Path("/g/data/gv90/xl1657/cmems_adt")
PROJ_DIR   = BASE_DIR / "meander_sst_project"
SST_DIR    = PROJ_DIR / "sst_regridded"
MASK_DIR   = PROJ_DIR / "masks"
TREND_DIR  = PROJ_DIR / "sst_trends"
FIG_DIR    = PROJ_DIR / "figures"
SI_DIR     = PROJ_DIR / "supporting_information"
SI_DIR.mkdir(parents=True, exist_ok=True)

SITES = {
    "CP":   {"name": "Campbell Plateau",       "lon_range": "150°E–150°W",
             "lat_range": "57°S–46°S", "wraps": True},
    "PAR":  {"name": "Pacific-Antarctic Ridge", "lon_range": "150°W–80°W",
             "lat_range": "60°S–48°S", "wraps": False},
    "SEIR": {"name": "Southeast Indian Ridge",  "lon_range": "130°E–152°E",
             "lat_range": "56°S–44°S", "wraps": False},
    "SWIR": {"name": "Southwest Indian Ridge",  "lon_range": "15°E–45°E",
             "lat_range": "58°S–44°S", "wraps": False},
}

CONTROLS = {
    "CTRL_SE_PAC":  {"name": "Southeast Pacific",
                     "lon_range": "100°W–80°W", "lat_range": "60°S–48°S",
                     "rationale": "Flat seafloor between PAR and Drake Passage; "
                                  "no major topographic features"},
    "CTRL_S_IND":   {"name": "Central South Indian",
                     "lon_range": "90°E–110°E", "lat_range": "56°S–44°S",
                     "rationale": "Quiescent ACC section between Kerguelen "
                                  "and SEIR meanders"},
    "CTRL_S_ATL":   {"name": "South Atlantic",
                     "lon_range": "10°W–10°E", "lat_range": "58°S–44°S",
                     "rationale": "Flat abyssal plain south of the Mid-Atlantic "
                                  "Ridge; minimal topographic interaction"},
    "CTRL_CP_SEIR": {"name": "CP–SEIR gap",
                     "lon_range": "180°–165°W", "lat_range": "57°S–46°S",
                     "rationale": "Transition zone between Campbell Plateau "
                                  "and Southeast Indian Ridge meanders"},
}

# Colour-blind-safe palette
CB = {"CP": "#0072B2", "PAR": "#E69F00", "SEIR": "#CC79A7", "SWIR": "#009E73"}

# Sub-periods for decadal decomposition
SUB_PERIODS = [
    ("1993–2005", 1993, 2005),
    ("2006–2015", 2006, 2015),
    ("2016–2025", 2016, 2025),
]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "axes.labelsize": 10,
    "figure.dpi": 150, "mathtext.default": "regular",
})


# =============================================================================
# TABLE S1: CONTROL REGION DEFINITIONS
# =============================================================================
def produce_table_s1():
    """Generate Table S1 in CSV and LaTeX format."""
    print("\nProducing Table S1: Control region definitions...")

    rows = []
    for ck, ctrl in CONTROLS.items():
        # Load mask to get grid point count
        fp = MASK_DIR / f"control_mask_{ck}.nc"
        n_pts = 0
        if fp.exists():
            ds = xr.open_dataset(fp)
            mask = ds["control_mask"].values.astype(bool)
            n_pts = mask.sum()
            ds.close()

        rows.append({
            "Region ID": ck,
            "Name": ctrl["name"],
            "Longitude range": ctrl["lon_range"],
            "Latitude range": ctrl["lat_range"],
            "Grid points": int(n_pts),
            "Selection rationale": ctrl["rationale"],
        })

    df = pd.DataFrame(rows)

    # Save CSV
    csv_fp = SI_DIR / "table_s1_control_regions.csv"
    df.to_csv(csv_fp, index=False)
    print(f"  Saved CSV: {csv_fp}")

    # Generate LaTeX table
    latex_fp = SI_DIR / "table_s1_control_regions.tex"
    with open(latex_fp, "w") as f:
        f.write("\\begin{table}[ht]\n")
        f.write("\\centering\n")
        f.write("\\caption{Definition of quiescent ACC control regions used for "
                "statistical comparison with standing meander sites. Each region "
                "was selected at longitudes not associated with major topographic "
                "features or documented standing meanders, spanning a latitude "
                "band comparable to the ACC core.}\n")
        f.write("\\label{tab:control_regions}\n")
        f.write("\\small\n")
        f.write("\\begin{tabular}{llllp{5.5cm}}\n")
        f.write("\\hline\n")
        f.write("Name & Longitude & Latitude & Grid points & "
                "Selection rationale \\\\\n")
        f.write("\\hline\n")
        for _, row in df.iterrows():
            f.write(f"{row['Name']} & {row['Longitude range']} & "
                    f"{row['Latitude range']} & {row['Grid points']:,} & "
                    f"{row['Selection rationale']} \\\\\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\end{table}\n")

    print(f"  Saved LaTeX: {latex_fp}")
    print(f"  Table S1:")
    print(df.to_string(index=False))

    return df


# =============================================================================
# TABLE S2: DETAILED PER-REGION SST TREND STATISTICS
# =============================================================================
def produce_table_s2():
    """Generate detailed statistics table for all regions."""
    print("\nProducing Table S2: Detailed regional SST trend statistics...")

    # Load the full trend stats CSV produced by NB14
    stats_fp = TREND_DIR / "sst_trend_stats.csv"
    if not stats_fp.exists():
        print(f"  ERROR: {stats_fp} not found. Run NB14 first.")
        return None

    df = pd.read_csv(stats_fp)
    csv_fp = SI_DIR / "table_s2_detailed_stats.csv"
    df.to_csv(csv_fp, index=False)
    print(f"  Saved: {csv_fp}")
    print(f"  Table S2:")
    print(df.to_string(index=False))

    return df


# =============================================================================
# FIGURE S1: DECADAL DECOMPOSITION
# =============================================================================
def load_sst_subperiod(year_start, year_end):
    """Load regridded monthly SST for a sub-period."""
    datasets = []
    for year in range(year_start, year_end + 1):
        fp = SST_DIR / f"sst_monthly_0125deg_{year}.nc"
        if not fp.exists():
            continue
        ds = xr.open_dataset(fp)
        datasets.append(ds["sst"])

    if len(datasets) == 0:
        return None
    sst_all = xr.concat(datasets, dim="time")
    for ds in datasets:
        ds.close()
    return sst_all


def compute_trends_subperiod(sst_all):
    """Compute per-grid-point SST trends for a sub-period."""
    lat = sst_all.latitude.values
    lon = sst_all.longitude.values
    nlat, nlon = len(lat), len(lon)
    sst_data = sst_all.values
    time_idx = pd.DatetimeIndex(sst_all.time.values)
    cal_months = time_idx.month

    slope_map = np.full((nlat, nlon), np.nan, dtype=np.float32)

    for i in range(nlat):
        for j in range(nlon):
            ts = sst_data[:, i, j]
            valid = np.isfinite(ts)
            if valid.sum() < 24:  # require 2 years
                continue

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
                acf_vals = acf(ts_anom, nlags=3, fft=True)
                lag1 = acf_vals[1]
            except Exception:
                lag1 = 0.0

            try:
                if abs(lag1) > 0.1:
                    result = mk.hamed_rao_modification_test(ts_anom)
                else:
                    result = mk.original_test(ts_anom)
                slope_map[i, j] = result.slope * 12 * 10  # per decade
            except Exception:
                continue

    return slope_map, lat, lon


def extract_meander_control_means(slope_map, lat, lon):
    """Extract area-weighted meander and control mean trends."""
    lat_2d = np.broadcast_to(lat[:, np.newaxis], slope_map.shape)
    cos_w = np.cos(np.deg2rad(lat_2d))

    meander_means = []
    for sk in ["CP", "PAR", "SEIR", "SWIR"]:
        fp = MASK_DIR / f"meander_envelope_mask_{sk}.nc"
        if not fp.exists():
            continue
        ds = xr.open_dataset(fp)
        mask = ds["meander_mask"].values.astype(bool)
        ds.close()
        t = slope_map[mask]
        w = cos_w[mask]
        valid = np.isfinite(t)
        if valid.sum() > 0:
            meander_means.append(np.average(t[valid], weights=w[valid]))

    control_means = []
    for ck in ["CTRL_SE_PAC", "CTRL_S_IND", "CTRL_S_ATL", "CTRL_CP_SEIR"]:
        fp = MASK_DIR / f"control_mask_{ck}.nc"
        if not fp.exists():
            continue
        ds = xr.open_dataset(fp)
        mask = ds["control_mask"].values.astype(bool)
        ds.close()
        t = slope_map[mask]
        w = cos_w[mask]
        valid = np.isfinite(t)
        if valid.sum() > 0:
            control_means.append(np.average(t[valid], weights=w[valid]))

    m_agg = np.mean(meander_means) if meander_means else np.nan
    c_agg = np.mean(control_means) if control_means else np.nan

    return m_agg, c_agg


def produce_figure_s1():
    """
    Figure S1: Decadal decomposition of the meander-vs-control SST trend
    contrast across three sub-periods (1993-2005, 2006-2015, 2016-2025).
    """
    print("\nProducing Figure S1: Decadal decomposition...")

    results = []
    for label, y0, y1 in SUB_PERIODS:
        print(f"\n  Processing {label} ({y0}–{y1})...")
        t0 = _time.time()

        sst = load_sst_subperiod(y0, y1)
        if sst is None:
            print(f"    No data for {label}")
            continue

        n_months = sst.sizes["time"]
        print(f"    Loaded: {n_months} months")

        slope_map, lat, lon = compute_trends_subperiod(sst)
        del sst
        gc.collect()

        m_mean, c_mean = extract_meander_control_means(slope_map, lat, lon)
        diff = m_mean - c_mean

        results.append({
            "period": label,
            "meander_mean": m_mean,
            "control_mean": c_mean,
            "difference": diff,
        })

        elapsed = (_time.time() - t0) / 60
        print(f"    Meander: {m_mean:+.4f}, Control: {c_mean:+.4f}, "
              f"Diff: {diff:+.4f} °C/dec ({elapsed:.1f} min)")

    if len(results) == 0:
        print("  ERROR: No sub-period results computed")
        return

    df = pd.DataFrame(results)
    df.to_csv(SI_DIR / "decadal_decomposition.csv", index=False)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    periods = [r["period"] for r in results]
    x = np.arange(len(periods))
    w = 0.35

    # Panel (a): Meander vs control trends per sub-period
    bars_m = ax1.bar(x - w/2, [r["meander_mean"] for r in results],
                     w, color="#0072B2", alpha=0.7, label="Meander sites")
    bars_c = ax1.bar(x + w/2, [r["control_mean"] for r in results],
                     w, color="0.6", alpha=0.7, label="Control regions")

    ax1.axhline(0, color="black", lw=0.4, alpha=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(periods, fontsize=9)
    ax1.set_ylabel("Mean SST trend (°C per decade)", fontsize=10)
    ax1.legend(fontsize=8, framealpha=0.9)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.text(0.02, 0.95, "(a)", transform=ax1.transAxes,
            fontsize=13, fontweight="bold", va="top")

    # Panel (b): Meander-control difference per sub-period
    diffs = [r["difference"] for r in results]
    colors = ["#E69F00" if d < 0 else "#0072B2" for d in diffs]
    ax2.bar(x, diffs, 0.5, color=colors, alpha=0.7, edgecolor="black", lw=0.5)
    ax2.axhline(0, color="black", lw=0.4, alpha=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(periods, fontsize=9)
    ax2.set_ylabel("Meander – Control SST trend\ndifference (°C per decade)",
                   fontsize=10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(0.02, 0.95, "(b)", transform=ax2.transAxes,
            fontsize=13, fontweight="bold", va="top")

    # Annotate each bar with value
    for i, d in enumerate(diffs):
        ax2.text(i, d + (0.005 if d >= 0 else -0.005),
                f"{d:+.03f}", ha="center",
                va="bottom" if d >= 0 else "top", fontsize=8)

    fig.tight_layout(w_pad=2.0)

    fp_out = SI_DIR / "fig_s1_decadal_decomposition.pdf"
    fig.savefig(fp_out, dpi=300, bbox_inches="tight")
    fig.savefig(fp_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {fp_out}")


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB17: Supporting Information for GRL Meander SST Trends Paper")
    print(f"  Output: {SI_DIR}")
    print("=" * 70)

    produce_table_s1()
    produce_table_s2()
    produce_figure_s1()

    print("\n" + "=" * 70)
    print("All Supporting Information outputs generated.")
    print(f"Output directory: {SI_DIR}")
    print("=" * 70)
