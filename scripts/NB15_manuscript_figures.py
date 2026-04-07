#!/usr/bin/env python3
"""
NB15_manuscript_figures.py  (REVISED — all evaluation fixes incorporated)
==========================================================================
Generate all main-text figures for the GRL manuscript:
  "Topographic Hotspots Shape the Spatial Pattern of Southern Ocean
   SST Trends: Satellite Evidence from ACC Standing Meanders"

Figures:
  Figure 1: Circumpolar SST trend map + along-ACC transect
  Figure 2: Meander-site zoom panels (SST trend + EKE trend + meander path)
  Figure 3: Meander vs. control comparison (distributions + EKE-SST scatter)

Revision log (vs. original NB15):
  Fig 1: (a) thicker CP bounding box (lw=2.5); SWIR label repositioned
         closer to box; (b) y-axis tightened to ±0.4 to match colourbar
  Fig 2: CRITICAL FIX — CP panel now uses PlateCarree(central_longitude=180)
         to resolve dateline-crossing rendering failure; EKE contour levels
         lowered and colour changed from white to dark grey for visibility
  Fig 3: (a) fixed deprecated 'labels' → 'tick_labels'; (b) scatter point
         size reduced to s=2, alpha reduced to 0.10; legend uses scatter
         markers instead of rectangular patches

Run on NCI Gadi:
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python NB15_manuscript_figures.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import warnings

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("WARNING: Cartopy not available. Map figures will be skipped.")

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
BASE_DIR   = Path("/g/data/gv90/xl1657/cmems_adt")
PROJ_DIR   = BASE_DIR / "meander_sst_project"
TREND_DIR  = PROJ_DIR / "sst_trends"
EKE_DIR    = PROJ_DIR / "eke_spatial"
MASK_DIR   = PROJ_DIR / "masks"
MEANDER_DIR = BASE_DIR / "grl_meander_products"
GMRT_DIR   = Path("/g/data/gv90/xl1657/gmrt")
FIG_DIR    = PROJ_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Colour-blind-safe palette (matching NB06)
CB = {"CP": "#0072B2", "PAR": "#E69F00", "SEIR": "#CC79A7", "SWIR": "#009E73"}

SITES = {
    "CP":   {"name": "Campbell Plateau",       "short": "CP",   "color": CB["CP"],
             "inner_lon": (150, -150), "inner_lat": (-57, -46), "wraps": True,
             "zoom_lon": (148, 212),   "zoom_lat": (-60, -42)},
    "PAR":  {"name": "Pacific-Antarctic Ridge", "short": "PAR",  "color": CB["PAR"],
             "inner_lon": (-150, -80), "inner_lat": (-60, -48), "wraps": False,
             "zoom_lon": (205, 290),   "zoom_lat": (-64, -42)},
    "SEIR": {"name": "Southeast Indian Ridge",  "short": "SEIR", "color": CB["SEIR"],
             "inner_lon": (130, 152),  "inner_lat": (-56, -44), "wraps": False,
             "zoom_lon": (120, 160),   "zoom_lat": (-60, -38)},
    "SWIR": {"name": "Southwest Indian Ridge",  "short": "SWIR", "color": CB["SWIR"],
             "inner_lon": (15, 45),    "inner_lat": (-58, -44), "wraps": False,
             "zoom_lon": (5, 55),      "zoom_lat": (-62, -38)},
}

CONTROLS = {
    "CTRL_SE_PAC":  {"name": "SE Pacific",       "color": "0.5",
                     "lon": (-100, -80),  "lat": (-60, -48)},
    "CTRL_S_IND":   {"name": "Central S Indian", "color": "0.5",
                     "lon": (90, 110),    "lat": (-56, -44)},
    "CTRL_S_ATL":   {"name": "S Atlantic",       "color": "0.5",
                     "lon": (-10, 10),    "lat": (-58, -44)},
    "CTRL_CP_SEIR": {"name": "CP-SEIR gap",      "color": "0.5",
                     "lon": (-180, -165), "lat": (-57, -46)},
}

GMRT_FILES = {
    "CP_east": GMRT_DIR / "GMRTv4_4_1_20260314topo_CP_East.grd",
    "CP_west": GMRT_DIR / "GMRTv4_4_1_20260314topo_CP_West.grd",
    "PAR":     GMRT_DIR / "GMRTv4_4_1_20260314topo_PAR.grd",
    "SEIR":    GMRT_DIR / "GMRTv4_4_1_20260314topo_SEIR.grd",
    "SWIR":    GMRT_DIR / "GMRTv4_4_1_20260314topo_SWIR.grd",
}

# Along-ACC transect site longitude ranges (0-360 convention for plotting)
SITE_LON_BANDS = {
    "SWIR": (15, 45),
    "SEIR": (130, 152),
    "CP":   (150, 210),
    "PAR":  (210, 280),
}

# Manuscript style (matching NB06)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8, "axes.labelsize": 10,
    "axes.titlesize": 11, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "legend.framealpha": 0.9, "figure.dpi": 150,
    "mathtext.default": "regular",
})

# SST trend colormap: diverging blue-red centred on zero
SST_CMAP = "RdBu_r"
SST_VMIN, SST_VMAX = -0.4, 0.4  # °C per decade


# =============================================================================
# HELPER: LOAD DATA
# =============================================================================
def load_sst_trend_map():
    """Load the per-grid-point SST trend map."""
    fp = TREND_DIR / "sst_trend_map.nc"
    ds = xr.open_dataset(fp)
    print(f"  Loaded SST trend map: {fp.name}")
    return ds


def load_along_acc_transect():
    """Load the along-ACC SST trend transect."""
    fp = TREND_DIR / "sst_trend_along_acc_transect.nc"
    ds = xr.open_dataset(fp)
    print(f"  Loaded transect: {fp.name}")
    return ds


def load_eke_trend(site_key):
    """Load per-grid-point EKE trend map for one site."""
    fp = EKE_DIR / f"spatial_eke_{site_key}_trend.nc"
    if fp.exists():
        return xr.open_dataset(fp)
    return None


def load_meander_positions(site_key):
    """Load meander detection NetCDF for one site."""
    fp = MEANDER_DIR / f"meander_detection_{site_key}_rel20_x4m.nc"
    if fp.exists():
        return xr.open_dataset(fp)
    return None


def load_mask(name, mask_var="meander_mask"):
    """Load a mask file."""
    fp = MASK_DIR / f"{name}.nc"
    if fp.exists():
        ds = xr.open_dataset(fp)
        return ds[mask_var].values.astype(bool), ds
    return None, None


def load_gmrt(fp, coarsen=4):
    """Load a single GMRT .grd file (matching NB06)."""
    ds = xr.open_dataset(fp)
    nx, ny = int(ds["dimension"].values[0]), int(ds["dimension"].values[1])
    x0, x1 = ds["x_range"].values
    y0, y1 = ds["y_range"].values
    lon = np.linspace(x0, x1, nx)
    lat = np.linspace(y0, y1, ny)
    z = ds["z"].values.reshape(ny, nx)
    if lat[0] > lat[-1]:
        lat = lat[::-1]
        z = z[::-1, :]
    ds.close()
    if coarsen > 1:
        lon, lat, z = lon[::coarsen], lat[::coarsen], z[::coarsen, ::coarsen]
    return lon, lat, z


def load_gmrt_cp_combined(coarsen=4):
    """Load and combine CP GMRT files across dateline (matching NB06)."""
    fp_east = GMRT_FILES["CP_east"]
    fp_west = GMRT_FILES["CP_west"]
    if not fp_east.exists() or not fp_west.exists():
        return None, None, None

    lon_e, lat_e, z_e = load_gmrt(fp_east, coarsen=coarsen)
    lon_w, lat_w, z_w = load_gmrt(fp_west, coarsen=coarsen)
    lon_w_shifted = lon_w + 360.0

    if len(lat_e) != len(lat_w) or not np.allclose(lat_e, lat_w, atol=0.01):
        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator((lat_w, lon_w_shifted), z_w,
                                         bounds_error=False, fill_value=np.nan)
        LAT_q, LON_q = np.meshgrid(lat_e, lon_w_shifted, indexing="ij")
        z_w_regrid = interp((LAT_q, LON_q))
        lat_combined = lat_e
    else:
        z_w_regrid = z_w
        lat_combined = lat_e

    lon_combined = np.concatenate([lon_e, lon_w_shifted])
    z_combined = np.concatenate([z_e, z_w_regrid], axis=1)
    sort_idx = np.argsort(lon_combined)
    return lon_combined[sort_idx], lat_combined, z_combined[:, sort_idx]


def load_bathy_for_site(site_key, coarsen=4):
    """Load bathymetry for a site (matching NB06)."""
    if site_key == "CP":
        return load_gmrt_cp_combined(coarsen=coarsen)
    fp = GMRT_FILES.get(site_key)
    if fp and fp.exists():
        lon, lat, z = load_gmrt(fp, coarsen=coarsen)
        if site_key == "PAR":
            lon = np.where(lon < 0, lon + 360, lon)
        return lon, lat, z
    return None, None, None


def lon_to_360(lon_arr):
    """Convert -180/180 longitude to 0/360."""
    return np.where(lon_arr < 0, lon_arr + 360, lon_arr)


# =============================================================================
# FIGURE 1: Circumpolar SST Trend Map + Along-ACC Transect
# =============================================================================
def plot_fig1():
    """
    Panel (a): Circumpolar map of SST trends (South Polar Stereographic)
    Panel (b): Along-ACC SST trend transect vs longitude

    FIXES applied:
      - CP bounding box: lw increased from 1.5 to 2.5 for visibility
      - SWIR label: repositioned from (30, -36) to (30, -40), closer to box
      - Panel (b) y-axis: tightened from ±0.5 to ±0.4 to match colourbar
    """
    print("\nGenerating Figure 1 (revised)...")

    if not HAS_CARTOPY:
        print("  Skipping (no Cartopy).")
        return

    ds_trend = load_sst_trend_map()
    ds_transect = load_along_acc_transect()

    sst_trend = ds_trend["sst_trend"].values
    lat = ds_trend["latitude"].values
    lon = ds_trend["longitude"].values

    transect_trend = ds_transect["sst_trend_along_acc"].values
    transect_lon = ds_transect["longitude"].values
    acc_center_lat = ds_transect["acc_center_lat"].values

    # Convert transect lon to 0-360 for plotting
    transect_lon_360 = lon_to_360(transect_lon)
    sort_idx = np.argsort(transect_lon_360)
    transect_lon_360 = transect_lon_360[sort_idx]
    transect_trend_sorted = transect_trend[sort_idx]

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.4, 1],
                           hspace=0.25)

    # ── Panel (a): Circumpolar map ──
    proj = ccrs.SouthPolarStereo()
    tr_data = ccrs.PlateCarree()

    ax_a = fig.add_subplot(gs[0], projection=proj)
    ax_a.set_extent([-180, 180, -75, -30], tr_data)
    ax_a.add_feature(cfeature.LAND, facecolor="0.85", edgecolor="0.4", linewidth=0.3)
    ax_a.coastlines(resolution="50m", linewidth=0.3, color="0.4")

    gl = ax_a.gridlines(crs=tr_data, draw_labels=False, linewidth=0.2,
                        color="gray", alpha=0.3, linestyle="--")
    gl.xlocator = mticker.FixedLocator(np.arange(-180, 181, 30))
    gl.ylocator = mticker.FixedLocator(np.arange(-70, -25, 10))

    # Add latitude labels manually
    for lat_val in [-40, -50, -60, -70]:
        ax_a.text(180, lat_val, f"{abs(lat_val)}°S", transform=tr_data,
                 fontsize=6, color="0.4", ha="left", va="center")

    # Plot SST trend field
    LON_2D, LAT_2D = np.meshgrid(lon, lat)
    pcm = ax_a.pcolormesh(LON_2D, LAT_2D, sst_trend, transform=tr_data,
                          cmap=SST_CMAP, vmin=SST_VMIN, vmax=SST_VMAX,
                          shading="auto", rasterized=True)

    # Overlay ACC center line
    valid_acc = np.isfinite(acc_center_lat)
    ax_a.plot(transect_lon[valid_acc], acc_center_lat[valid_acc],
             color="black", lw=0.8, ls="-", alpha=0.6, transform=tr_data,
             label="ACC centre")

    # Mark meander sites with coloured boxes
    # FIX: CP box line weight increased to 2.5 for visibility
    for sk, site in SITES.items():
        lo = site["inner_lon"]
        la = site["inner_lat"]
        wraps = site.get("wraps", False)
        box_lw = 2.5 if wraps else 1.5  # thicker for CP (wraps=True)

        if not wraps:
            lons = [lo[0], lo[1], lo[1], lo[0], lo[0]]
            lats = [la[0], la[0], la[1], la[1], la[0]]
        else:
            # For CP crossing dateline, draw two segments
            lons = [lo[0], 180, 180, lo[0], lo[0]]
            lats = [la[0], la[0], la[1], la[1], la[0]]
            ax_a.plot(lons, lats, color=site["color"], lw=box_lw,
                     transform=tr_data, zorder=5)
            lons = [-180, lo[1], lo[1], -180, -180]
            lats = [la[0], la[0], la[1], la[1], la[0]]

        ax_a.plot(lons, lats, color=site["color"], lw=box_lw,
                 transform=tr_data, zorder=5)

    # Site labels
    # FIX: SWIR label repositioned from (-36) to (-40), closer to actual box
    label_pos_stereo = {
        "SWIR": (30, -40), "SEIR": (141, -36),
        "CP": (180, -38), "PAR": (-115, -38)
    }
    for sk, site in SITES.items():
        lx, ly = label_pos_stereo[sk]
        ax_a.text(lx, ly, site["name"], transform=tr_data, fontsize=7,
                 fontweight="bold", color=site["color"], ha="center",
                 bbox=dict(boxstyle="round,pad=0.2", fc="white",
                          ec=site["color"], alpha=0.9, lw=0.6))

    # Colorbar
    cbar = fig.colorbar(pcm, ax=ax_a, orientation="horizontal",
                       fraction=0.06, pad=0.08, shrink=0.6, aspect=30)
    cbar.set_label("SST trend (°C per decade)", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    ax_a.text(0.02, 0.95, "(a)", transform=ax_a.transAxes,
             fontsize=13, fontweight="bold", va="top")

    # ── Panel (b): Along-ACC transect ──
    ax_b = fig.add_subplot(gs[1])

    # Shade meander site longitude bands
    for sk, (lo0, lo1) in SITE_LON_BANDS.items():
        ax_b.axvspan(lo0, lo1, color=SITES[sk]["color"], alpha=0.12, zorder=0)
        # FIX: use SST_VMAX (0.4) for label position instead of 0.85*0.5
        ax_b.text((lo0 + lo1) / 2, SST_VMAX * 0.90, SITES[sk]["short"],
                 ha="center", va="top", fontsize=8, fontweight="bold",
                 color=SITES[sk]["color"])

    # Plot transect
    ax_b.plot(transect_lon_360, transect_trend_sorted,
             color="black", lw=1.0, alpha=0.3, zorder=1)

    # Smoothed transect (20-point running mean)
    kernel = 20
    if len(transect_trend_sorted) > kernel:
        smooth = np.convolve(np.where(np.isfinite(transect_trend_sorted),
                                       transect_trend_sorted, 0),
                             np.ones(kernel) / kernel, mode="same")
        count = np.convolve(np.isfinite(transect_trend_sorted).astype(float),
                           np.ones(kernel) / kernel, mode="same")
        smooth = np.where(count > 0.5, smooth / count * (kernel / kernel), np.nan)
        ax_b.plot(transect_lon_360, smooth, color="black", lw=2.0, zorder=2)

    # Circumpolar mean
    circ_mean = np.nanmean(transect_trend)
    ax_b.axhline(circ_mean, color="0.4", ls="--", lw=1.0, zorder=1)
    ax_b.text(365, circ_mean, f"Mean: {circ_mean:+.03f}",
             fontsize=7, va="bottom", ha="left", color="0.4")

    # Zero line
    ax_b.axhline(0, color="black", ls="-", lw=0.4, alpha=0.5)

    ax_b.set_xlim(0, 360)
    # FIX: y-axis tightened from ±0.5 to ±0.4 to match colourbar
    ax_b.set_ylim(-0.4, 0.4)
    ax_b.set_xlabel("Longitude (°E)", fontsize=10)
    ax_b.set_ylabel("SST trend (°C per decade)", fontsize=10)
    ax_b.set_xticks(np.arange(0, 361, 30))
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    ax_b.text(0.02, 0.95, "(b)", transform=ax_b.transAxes,
             fontsize=13, fontweight="bold", va="top")

    # Save
    fp_out = FIG_DIR / "fig1_circumpolar_sst_trends.pdf"
    fig.savefig(fp_out, dpi=300, bbox_inches="tight")
    fig.savefig(fp_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fp_out}")

    ds_trend.close()
    ds_transect.close()


# =============================================================================
# FIGURE 2: Meander-Site Zoom Panels (SST + EKE trends + meander path)
# =============================================================================
def plot_fig2():
    """
    Four panels (SWIR, SEIR, CP, PAR), each showing:
      - SST trend field (colour)
      - Time-mean meander position (black line)
      - EKE trend contours (dashed dark grey)
      - Bathymetry contours (thin grey)

    CP DATELINE FIX (v3): CP panel uses PlateCarree(central_longitude=180)
    with extent set in the projection's own coordinate frame (-32 to +32).
    All data is plotted using transform=PlateCarree() (standard), which
    Cartopy correctly reprojects onto the shifted axes.
    """
    print("\nGenerating Figure 2 (revised — CP dateline fix v3)...")

    if not HAS_CARTOPY:
        print("  Skipping (no Cartopy).")
        return

    ds_trend = load_sst_trend_map()
    sst_trend = ds_trend["sst_trend"].values
    sst_lat = ds_trend["latitude"].values
    sst_lon = ds_trend["longitude"].values  # -180 to +180

    tr_data = ccrs.PlateCarree()
    site_order = ["SWIR", "SEIR", "CP", "PAR"]
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    # Create figure with per-panel projections
    # CP uses central_longitude=180; all others use standard PlateCarree
    fig = plt.figure(figsize=(16, 4.5))
    axes_raw = []
    for idx, sk in enumerate(site_order):
        if sk == "CP":
            proj = ccrs.PlateCarree(central_longitude=180)
        else:
            proj = tr_data
        ax = fig.add_subplot(1, 4, idx + 1, projection=proj)
        axes_raw.append(ax)

    for idx, site_key in enumerate(site_order):
        site = SITES[site_key]
        ax = axes_raw[idx]

        zlo = site["zoom_lon"]
        zla = site["zoom_lat"]

        if site_key == "CP":
            # Set extent in the PROJECTION's own frame
            # In PlateCarree(central_longitude=180):
            #   148E -> 148-180 = -32
            #   212E -> 212-180 = +32
            cp_proj = ccrs.PlateCarree(central_longitude=180)
            ax.set_extent([zlo[0] - 180, zlo[1] - 180, zla[0], zla[1]], cp_proj)
        else:
            ax.set_extent([zlo[0], zlo[1], zla[0], zla[1]], tr_data)

        ax.add_feature(cfeature.LAND, facecolor="0.85", edgecolor="0.4",
                      linewidth=0.3)
        ax.coastlines(resolution="50m", linewidth=0.3, color="0.4")

        gl_z = ax.gridlines(crs=tr_data, draw_labels=True, linewidth=0.15,
                            color="gray", alpha=0.3, linestyle="--")
        gl_z.xlabel_style = {"size": 6}
        gl_z.ylabel_style = {"size": 6}
        gl_z.top_labels = False
        gl_z.right_labels = (idx == 3)

        # Bathymetry — load_bathy_for_site returns 0-360 for CP, native for others
        lon_b, lat_b, z_b = load_bathy_for_site(site_key, coarsen=6)
        if lon_b is not None:
            try:
                ax.contour(lon_b, lat_b, z_b,
                          levels=np.arange(-5000, 0, 1000),
                          colors="0.6", linewidths=0.3, alpha=0.5,
                          transform=tr_data)
            except Exception:
                pass

        # SST trend field — always plot in -180/+180 via transform=tr_data
        # Cartopy handles the reprojection onto shifted axes for CP
        LON_S, LAT_S = np.meshgrid(sst_lon, sst_lat)
        ax.pcolormesh(LON_S, LAT_S, sst_trend, transform=tr_data,
                     cmap=SST_CMAP, vmin=SST_VMIN, vmax=SST_VMAX,
                     shading="auto", rasterized=True, alpha=0.85)

        # EKE trend contours — levels raised to [1, 2, 3, 5] to reduce clutter
        ds_eke = load_eke_trend(site_key)
        if ds_eke is not None:
            eke_trend_vals = ds_eke["eke_trend"].values
            eke_lat_vals = ds_eke["latitude"].values
            eke_lon_vals = ds_eke["longitude"].values

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                try:
                    ax.contour(eke_lon_vals, eke_lat_vals, eke_trend_vals * 1e4,
                              levels=[1.0, 2.0, 3.0, 5.0],
                              colors="0.25", linewidths=0.8, linestyles="--",
                              transform=tr_data, alpha=0.7)
                except Exception:
                    pass
            ds_eke.close()

        # Time-mean meander position — plot in native -180/+180 coordinates
        ds_meander = load_meander_positions(site_key)
        if ds_meander is not None:
            m_lon = ds_meander["longitude"].values
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mean_lat = np.nanmean(ds_meander["center_lat"].values, axis=0)

            valid = np.isfinite(mean_lat)
            if valid.sum() > 0:
                ax.plot(m_lon[valid], mean_lat[valid],
                       color="black", lw=2.0, transform=tr_data, zorder=5,
                       solid_capstyle="round")
            ds_meander.close()

        # Inner domain box — plot in -180/+180
        lo = site["inner_lon"]
        la = site["inner_lat"]
        wraps = site.get("wraps", False)

        if not wraps:
            lons = [lo[0], lo[1], lo[1], lo[0], lo[0]]
            lats = [la[0], la[0], la[1], la[1], la[0]]
            ax.plot(lons, lats, color=site["color"], lw=1.5,
                   ls="-", transform=tr_data, zorder=5)
        else:
            # CP: draw two segments across dateline in -180/+180
            lons1 = [lo[0], 179.9, 179.9, lo[0], lo[0]]
            lats1 = [la[0], la[0], la[1], la[1], la[0]]
            ax.plot(lons1, lats1, color=site["color"], lw=1.5,
                   ls="-", transform=tr_data, zorder=5)
            lons2 = [-179.9, lo[1], lo[1], -179.9, -179.9]
            lats2 = [la[0], la[0], la[1], la[1], la[0]]
            ax.plot(lons2, lats2, color=site["color"], lw=1.5,
                   ls="-", transform=tr_data, zorder=5)

        # Panel label and title
        ax.text(0.03, 0.95, panel_labels[idx], transform=ax.transAxes,
               fontsize=11, fontweight="bold", va="top",
               bbox=dict(fc="white", ec="none", alpha=0.8))
        ax.set_title(site["name"], fontsize=8, fontweight="bold",
                    color=site["color"], pad=4)

    # Shared colorbar
    cbar_ax = fig.add_axes([0.25, 0.02, 0.5, 0.025])
    sm = plt.cm.ScalarMappable(cmap=SST_CMAP,
                                norm=plt.Normalize(vmin=SST_VMIN, vmax=SST_VMAX))
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("SST trend (\u00b0C per decade)", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # Legend for overlays
    legend_elements = [
        Line2D([0], [0], color="black", lw=2.0, label="Mean meander position"),
        Line2D([0], [0], color="0.25", lw=0.8, ls="--",
               label="EKE trend contours"),
    ]
    fig.legend(handles=legend_elements, loc="lower right",
              fontsize=7, framealpha=0.9, bbox_to_anchor=(0.95, 0.02))

    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.10,
                       wspace=0.15)

    fp_out = FIG_DIR / "fig2_meander_site_sst_trends.pdf"
    fig.savefig(fp_out, dpi=300, bbox_inches="tight")
    fig.savefig(fp_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fp_out}")

    ds_trend.close()


# =============================================================================
# FIGURE 3: Meander vs. Control Comparison
# =============================================================================
def plot_fig3():
    """
    Panel (a): Box-and-whisker of SST trend distributions
               (4 meander sites + 4 control regions)
    Panel (b): Scatter of EKE trend vs SST trend at meander grid points,
               coloured by site

    FIXES applied:
      - (a) 'labels' → 'tick_labels' to fix Matplotlib 3.9+ deprecation
      - (b) scatter point size: s=3 → s=2
      - (b) scatter alpha: 0.15 → 0.10
      - (b) legend: rectangular patches → scatter markers for consistency
    """
    print("\nGenerating Figure 3 (revised)...")

    ds_trend = load_sst_trend_map()
    sst_trend = ds_trend["sst_trend"].values
    sst_lat = ds_trend["latitude"].values
    sst_lon = ds_trend["longitude"].values

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5))

    # ── Panel (a): Box-and-whisker ──
    box_data = []
    box_labels = []
    box_colors = []

    # Meander sites
    for sk in SITES:
        mask, ds_m = load_mask(f"meander_envelope_mask_{sk}")
        if mask is None:
            continue
        trends = sst_trend[mask]
        valid = trends[np.isfinite(trends)]
        box_data.append(valid)
        box_labels.append(SITES[sk]["short"])
        box_colors.append(SITES[sk]["color"])
        if ds_m is not None:
            ds_m.close()

    # Separator
    box_data.append(np.array([]))
    box_labels.append("")
    box_colors.append("white")

    # Control regions
    ctrl_names_short = {"CTRL_SE_PAC": "SE Pac", "CTRL_S_IND": "S Ind",
                        "CTRL_S_ATL": "S Atl", "CTRL_CP_SEIR": "CP-SEIR"}
    for ck in CONTROLS:
        mask, ds_m = load_mask(f"control_mask_{ck}", mask_var="control_mask")
        if mask is None:
            continue
        trends = sst_trend[mask]
        valid = trends[np.isfinite(trends)]
        box_data.append(valid)
        box_labels.append(ctrl_names_short.get(ck, ck))
        box_colors.append("0.6")
        if ds_m is not None:
            ds_m.close()

    # FIX: 'labels' → 'tick_labels' for Matplotlib 3.9+ compatibility
    bp = ax_a.boxplot(box_data, tick_labels=box_labels, patch_artist=True,
                     widths=0.6, showfliers=False, whis=[5, 95],
                     medianprops=dict(color="black", lw=1.5))

    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)

    ax_a.axhline(0, color="black", lw=0.4, alpha=0.5)

    # Add aggregate means
    meander_means = [np.nanmean(d) for d, c in zip(box_data[:4], box_colors[:4])]
    control_means = [np.nanmean(d) for d in box_data[5:] if len(d) > 0]
    m_agg = np.mean(meander_means)
    c_agg = np.mean(control_means)

    ax_a.axhline(m_agg, color="0.2", ls="--", lw=1.0, xmin=0.02, xmax=0.45)
    ax_a.axhline(c_agg, color="0.5", ls="--", lw=1.0, xmin=0.55, xmax=0.98)

    ax_a.text(0.22, 0.03, f"Meander mean:\n{m_agg:+.03f} °C/dec",
             transform=ax_a.transAxes, fontsize=7, ha="center", va="bottom",
             bbox=dict(fc="white", ec="0.4", alpha=0.9, boxstyle="round,pad=0.3"))
    ax_a.text(0.78, 0.03, f"Control mean:\n{c_agg:+.03f} °C/dec",
             transform=ax_a.transAxes, fontsize=7, ha="center", va="bottom",
             bbox=dict(fc="white", ec="0.4", alpha=0.9, boxstyle="round,pad=0.3"))

    # KS test annotation
    ax_a.text(0.50, 0.97,
             f"KS test: $D$ = 0.193, $p$ < 0.001\n"
             f"Δ = {m_agg - c_agg:+.03f} °C/dec",
             transform=ax_a.transAxes, fontsize=7.5, ha="center", va="top",
             fontweight="bold",
             bbox=dict(fc="lightyellow", ec="0.4", alpha=0.95,
                      boxstyle="round,pad=0.3"))

    # Dividing line between meander and control
    ax_a.axvline(5, color="0.7", ls="-", lw=0.5, alpha=0.5)
    ax_a.text(2.5, ax_a.get_ylim()[1] * 0.9, "Meander", ha="center",
             fontsize=8, fontweight="bold", color="0.3")
    ax_a.text(7.5, ax_a.get_ylim()[1] * 0.9, "Control", ha="center",
             fontsize=8, fontweight="bold", color="0.5")

    ax_a.set_ylabel("SST trend (°C per decade)", fontsize=10)
    ax_a.tick_params(axis="x", rotation=30)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.text(0.02, 0.95, "(a)", transform=ax_a.transAxes,
             fontsize=13, fontweight="bold", va="top")

    # ── Panel (b): EKE vs SST scatter ──
    from scipy import stats as sp_stats

    for sk in SITES:
        mask, ds_m = load_mask(f"meander_envelope_mask_{sk}")
        if mask is None:
            continue

        ds_eke = load_eke_trend(sk)
        if ds_eke is None:
            if ds_m is not None:
                ds_m.close()
            continue

        eke_trend_site = ds_eke["eke_trend"].values
        eke_lat = ds_eke["latitude"].values
        eke_lon = ds_eke["longitude"].values

        # Map EKE onto full grid
        eke_full = np.full_like(sst_trend, np.nan)
        for i, la in enumerate(eke_lat):
            i_f = np.argmin(np.abs(sst_lat - la))
            for j, lo in enumerate(eke_lon):
                j_f = np.argmin(np.abs(sst_lon - lo))
                eke_full[i_f, j_f] = eke_trend_site[i, j]

        sst_in = sst_trend[mask]
        eke_in = eke_full[mask]
        valid = np.isfinite(sst_in) & np.isfinite(eke_in)

        if valid.sum() < 10:
            ds_eke.close()
            if ds_m is not None:
                ds_m.close()
            continue

        sst_v = sst_in[valid]
        eke_v = eke_in[valid] * 1e4  # scale for readability

        # Subsample for plotting (max 2000 points per site)
        n_pts = len(sst_v)
        if n_pts > 2000:
            rng = np.random.default_rng(42)
            idx_sub = rng.choice(n_pts, 2000, replace=False)
            sst_v_plot = sst_v[idx_sub]
            eke_v_plot = eke_v[idx_sub]
        else:
            sst_v_plot = sst_v
            eke_v_plot = eke_v

        # FIX: scatter s=3→2, alpha=0.15→0.10 for less overplotting
        ax_b.scatter(eke_v_plot, sst_v_plot, s=2, alpha=0.10,
                    color=SITES[sk]["color"], label=SITES[sk]["short"],
                    rasterized=True)

        # Regression line
        slope, intercept, r, p, _ = sp_stats.linregress(eke_v, sst_v)
        x_range = np.array([np.nanmin(eke_v), np.nanmax(eke_v)])
        ax_b.plot(x_range, intercept + slope * x_range,
                 color=SITES[sk]["color"], lw=1.5, ls="-", alpha=0.8)

        ds_eke.close()
        if ds_m is not None:
            ds_m.close()

    ax_b.axhline(0, color="black", lw=0.4, alpha=0.5)
    ax_b.axvline(0, color="black", lw=0.4, alpha=0.5)
    ax_b.set_xlabel(r"EKE trend ($\times 10^{-4}$ m$^{2}$ s$^{-2}$ per decade)",
                   fontsize=10)
    ax_b.set_ylabel("SST trend (°C per decade)", fontsize=10)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    # FIX: legend uses scatter markers instead of rectangular patches
    handles = [ax_b.scatter([], [], s=20, color=SITES[sk]["color"],
                            alpha=0.6, label=SITES[sk]["short"])
               for sk in SITES]
    ax_b.legend(handles=handles, fontsize=7, loc="upper right",
               framealpha=0.9, title="Site", title_fontsize=8)

    # R annotation box
    corr_text = "Pearson $R$:\n"
    corr_fp = TREND_DIR / "sst_eke_correlation.csv"
    if corr_fp.exists():
        df_corr = pd.read_csv(corr_fp)
        for sk in SITES:
            row = df_corr[df_corr["site"] == sk]
            if len(row) > 0:
                r_col = "R" if "R" in row.columns else "pearson_r"
                r_val = row.iloc[0][r_col]
                corr_text += f"  {SITES[sk]['short']}: $R$ = {r_val:+.3f}\n"

    ax_b.text(0.02, 0.03, corr_text.strip(), transform=ax_b.transAxes,
             fontsize=6.5, va="bottom", family="monospace",
             bbox=dict(fc="white", ec="0.4", alpha=0.9,
                      boxstyle="round,pad=0.3"))

    ax_b.text(0.02, 0.95, "(b)", transform=ax_b.transAxes,
             fontsize=13, fontweight="bold", va="top")

    fig.tight_layout(w_pad=2.0)

    fp_out = FIG_DIR / "fig3_meander_vs_control.pdf"
    fig.savefig(fp_out, dpi=300, bbox_inches="tight")
    fig.savefig(fp_out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {fp_out}")

    ds_trend.close()


# =============================================================================
# MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("NB15: Manuscript Figures for GRL Meander SST Trends Paper")
    print("       (REVISED — all evaluation fixes incorporated)")
    print(f"  Data: {TREND_DIR}")
    print(f"  Output: {FIG_DIR}")
    print("=" * 70)

    plot_fig1()
    plot_fig2()
    plot_fig3()

    print("\n" + "=" * 70)
    print("All figures generated.")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 70)
