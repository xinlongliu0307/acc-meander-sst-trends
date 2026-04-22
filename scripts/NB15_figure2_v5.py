#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NB15_figure2_v5.py
===================

Generates Figure 2 for the JGR:Oceans manuscript. This version (v5)
adds targeted EKE contour clipping refinements relative to v4:

(1) Campbell Plateau (panel c): EKE contours are now clipped strictly
    to the longitude range 150 degrees East to 150 degrees West,
    matching the envelope longitude boundary. Any contours that
    appeared outside this range in the previous version due to
    meshgrid floating-point edge effects are now excluded.

(2) Pacific-Antarctic Ridge (panel d): EKE contours are now clipped
    to an eastern boundary of 82 degrees West rather than the
    envelope boundary of 80 degrees West. This 2-degree buffer
    excludes the noise contours that appeared near the South
    American coast at approximately 80 to 65 degrees West, which
    represent boundary artefacts rather than physically meaningful
    EKE trends.

All other aspects of the figure remain identical to v4, including
the envelope fill rendering, longitude label placement on all four
panels, and the reduced row spacing.

Author: Xinlong Liu
Target: JGR:Oceans
"""

import warnings
from pathlib import Path

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LongitudeFormatter, LatitudeFormatter

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE = Path("/g/data/gv90/xl1657/cmems_adt/meander_sst_project")
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SST_TREND_MAP = BASE / "sst_trends" / "sst_trend_map.nc"
EKE_DIR = BASE / "eke_spatial"

# SITES dictionary now includes a separate "eke_clip_lon" key that
# specifies the longitude range for EKE contour clipping, which can
# differ from the envelope_lon when tighter clipping is required.
SITES = {
    "SWIR": {
        "name": "Southwest Indian Ridge",
        "envelope_lon": (15.0, 45.0),
        "envelope_lat": (-58.0, -44.0),
        "eke_clip_lon": (15.0, 45.0),
        "eke_clip_lat": (-58.0, -44.0),
        "zoom_extent": [0.0, 60.0, -64.0, -38.0],
        "colour": "#2ca02c",
        "projection_lon": 0.0,
        "eke_trend_file": EKE_DIR / "spatial_eke_SWIR_trend.nc",
        "panel_label": "(a)",
        "xticks": [0, 15, 30, 45, 60],
    },
    "SEIR": {
        "name": "Southeast Indian Ridge",
        "envelope_lon": (130.0, 152.0),
        "envelope_lat": (-56.0, -44.0),
        "eke_clip_lon": (130.0, 152.0),
        "eke_clip_lat": (-56.0, -44.0),
        "zoom_extent": [115.0, 165.0, -62.0, -38.0],
        "colour": "#e377c2",
        "projection_lon": 0.0,
        "eke_trend_file": EKE_DIR / "spatial_eke_SEIR_trend.nc",
        "panel_label": "(b)",
        "xticks": [120, 135, 150, 165],
    },
    "CP": {
        "name": "Campbell Plateau",
        "envelope_lon": (150.0, -150.0),
        "envelope_lat": (-57.0, -46.0),
        # v5 change: strict clip to 150 E – 150 W for EKE contours
        "eke_clip_lon": (150.0, -150.0),
        "eke_clip_lat": (-57.0, -46.0),
        "zoom_extent": [-35.0, 35.0, -64.0, -40.0],
        "colour": "#1f77b4",
        "projection_lon": 180.0,
        "eke_trend_file": EKE_DIR / "spatial_eke_CP_trend.nc",
        "panel_label": "(c)",
        "xticks": [-30, -15, 0, 15, 30],
    },
    "PAR": {
        "name": "Pacific-Antarctic Ridge",
        "envelope_lon": (-150.0, -80.0),
        "envelope_lat": (-60.0, -48.0),
        # v5 change: tighter eastern boundary at 82 W instead of 80 W
        # to exclude coastal South American noise contours
        "eke_clip_lon": (-150.0, -82.0),
        "eke_clip_lat": (-60.0, -48.0),
        "zoom_extent": [-165.0, -65.0, -67.0, -43.0],
        "colour": "#ff7f0e",
        "projection_lon": 0.0,
        "eke_trend_file": EKE_DIR / "spatial_eke_PAR_trend.nc",
        "panel_label": "(d)",
        "xticks": [-160, -140, -120, -100, -80],
    },
}

PANEL_ORDER = ["SWIR", "SEIR", "CP", "PAR"]

SST_TREND_VMIN, SST_TREND_VMAX = -0.4, 0.4
SST_TREND_CMAP = "RdBu_r"

EKE_CONTOUR_LEVELS = [1e-4, 2e-4, 3e-4, 5e-4]

ENVELOPE_FILL_ALPHA = 0.18

FS_PANEL_TITLE = 15
FS_TICK = 13
FS_CBAR_LABEL = 14
FS_CBAR_TICK = 13
FS_LEGEND = 13


# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

def draw_envelope_with_fill(ax, lon_range, lat_range, colour,
                            linewidth=3.2, zorder=6, fill_alpha=0.18):
    """Draw the meander envelope as a filled region with coloured border."""
    lon_min, lon_max = lon_range
    lat_min, lat_max = lat_range

    if lon_max < lon_min:
        # Dateline crossing
        for (a, b) in [(lon_min, 180.0), (-180.0, lon_max)]:
            rect_fill = Rectangle(
                (a, lat_min), b - a, lat_max - lat_min,
                fill=True, facecolor=colour,
                edgecolor="none", alpha=fill_alpha,
                zorder=zorder - 1,
                transform=ccrs.PlateCarree(),
            )
            ax.add_patch(rect_fill)
            rect_border = Rectangle(
                (a, lat_min), b - a, lat_max - lat_min,
                fill=False, edgecolor=colour,
                linewidth=linewidth, linestyle="-",
                zorder=zorder,
                transform=ccrs.PlateCarree(),
            )
            ax.add_patch(rect_border)
    else:
        rect_fill = Rectangle(
            (lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
            fill=True, facecolor=colour,
            edgecolor="none", alpha=fill_alpha,
            zorder=zorder - 1,
            transform=ccrs.PlateCarree(),
        )
        ax.add_patch(rect_fill)
        rect_border = Rectangle(
            (lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
            fill=False, edgecolor=colour,
            linewidth=linewidth, linestyle="-",
            zorder=zorder,
            transform=ccrs.PlateCarree(),
        )
        ax.add_patch(rect_border)


def clip_eke_to_region(eke_lon, eke_lat, eke_trend, clip_lon, clip_lat):
    """
    Return a copy of eke_trend with values outside the specified
    longitude and latitude range set to NaN, so that contours are
    drawn only inside the clip region.
    """
    lon2d, lat2d = np.meshgrid(eke_lon, eke_lat)
    lat_min, lat_max = clip_lat
    lon_min, lon_max = clip_lon

    lat_mask = (lat2d >= lat_min) & (lat2d <= lat_max)

    # Handle longitude wrapping: convert any EKE longitudes in the
    # 0-360 range to -180 to 180 if necessary
    lon2d_norm = ((lon2d + 180.0) % 360.0) - 180.0

    if lon_max < lon_min:
        # Dateline crossing: inside means lon >= lon_min OR lon <= lon_max
        lon_mask = (lon2d_norm >= lon_min) | (lon2d_norm <= lon_max)
    else:
        lon_mask = (lon2d_norm >= lon_min) & (lon2d_norm <= lon_max)

    inside_region = lat_mask & lon_mask
    clipped = np.where(inside_region, eke_trend, np.nan)
    return clipped


def safe_open_dataset(path):
    if path is None or not Path(path).exists():
        print(f"  [warning] file not found: {path}")
        return None
    try:
        return xr.open_dataset(path)
    except Exception as e:
        print(f"  [warning] failed to open {path}: {e}")
        return None


def detect_eke_variable(eke_ds):
    candidates = [
        "eke_trend", "trend", "slope", "eke_slope",
        "sen_slope", "eke_sen_slope",
    ]
    for v in candidates:
        if v in eke_ds.data_vars:
            return v
    if len(eke_ds.data_vars) > 0:
        return list(eke_ds.data_vars)[0]
    return None


def add_longitude_latitude_labels(ax, xticks, projection):
    """Add explicit longitude and latitude labels."""
    ax.set_xticks(xticks, crs=projection)

    extent = ax.get_extent(crs=projection)
    lat_min_ext, lat_max_ext = extent[2], extent[3]
    lat_start = int(np.ceil(lat_min_ext / 5.0) * 5)
    lat_end = int(np.floor(lat_max_ext / 5.0) * 5)
    lat_ticks = list(range(lat_start, lat_end + 1, 5))
    ax.set_yticks(lat_ticks, crs=ccrs.PlateCarree())

    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    ax.tick_params(axis="both", labelsize=FS_TICK)


# ==============================================================================
# FIGURE 2
# ==============================================================================

def make_figure2():
    print("=" * 70)
    print("Generating Figure 2 (v5): Meander-Site SST Trend Patterns")
    print("=" * 70)

    ds_sst = xr.open_dataset(SST_TREND_MAP)
    trend_map = ds_sst["sst_trend"].values
    lat_map = ds_sst["latitude"].values
    lon_map = ds_sst["longitude"].values

    print(f"  SST trend map shape: {trend_map.shape}")

    fig = plt.figure(figsize=(16, 12))
    fig.subplots_adjust(
        left=0.06, right=0.97,
        top=0.95, bottom=0.16,
        wspace=0.12, hspace=0.08,
    )

    mesh_handle = None

    for idx, key in enumerate(PANEL_ORDER):
        site = SITES[key]
        row = idx // 2
        col = idx % 2
        subplot_idx = row * 2 + col + 1

        projection = ccrs.PlateCarree(
            central_longitude=site["projection_lon"]
        )

        ax = fig.add_subplot(2, 2, subplot_idx, projection=projection)
        ax.set_extent(site["zoom_extent"], crs=projection)

        # SST trend field
        mesh = ax.pcolormesh(
            lon_map, lat_map, trend_map,
            cmap=SST_TREND_CMAP,
            vmin=SST_TREND_VMIN, vmax=SST_TREND_VMAX,
            transform=ccrs.PlateCarree(),
            shading="auto", rasterized=True,
            zorder=1,
        )
        mesh_handle = mesh

        ax.add_feature(
            cfeature.LAND, facecolor="#f0f0f0",
            edgecolor="black", linewidth=0.5, zorder=3,
        )
        ax.coastlines(resolution="50m", linewidth=0.6, zorder=4)

        # Envelope fill and border
        draw_envelope_with_fill(
            ax,
            site["envelope_lon"], site["envelope_lat"],
            colour=site["colour"],
            linewidth=3.2, zorder=6,
            fill_alpha=ENVELOPE_FILL_ALPHA,
        )

        # EKE trend contours clipped to the EKE clip region
        eke_ds = safe_open_dataset(site["eke_trend_file"])
        if eke_ds is not None:
            eke_var = detect_eke_variable(eke_ds)
            if eke_var is not None:
                print(f"  {key}: using EKE variable '{eke_var}', "
                      f"clip lon {site['eke_clip_lon']}, "
                      f"clip lat {site['eke_clip_lat']}")
                eke_trend = eke_ds[eke_var].values

                lon_name = "longitude" if "longitude" in eke_ds.coords else "lon"
                lat_name = "latitude" if "latitude" in eke_ds.coords else "lat"
                eke_lon = eke_ds[lon_name].values
                eke_lat = eke_ds[lat_name].values

                # v5: use eke_clip_lon/lat instead of envelope_lon/lat
                clipped = clip_eke_to_region(
                    eke_lon, eke_lat, eke_trend,
                    site["eke_clip_lon"], site["eke_clip_lat"],
                )

                try:
                    ax.contour(
                        eke_lon, eke_lat, clipped,
                        levels=EKE_CONTOUR_LEVELS,
                        colors="#222222",
                        linewidths=1.1,
                        linestyles="--",
                        transform=ccrs.PlateCarree(),
                        zorder=5,
                    )
                except Exception as e:
                    print(f"  [warning] EKE contour failed for {key}: {e}")
            eke_ds.close()

        # Longitude and latitude labels
        try:
            add_longitude_latitude_labels(ax, site["xticks"], projection)
        except Exception as e:
            print(f"  [warning] label setup failed for {key}: {e}")
            gl = ax.gridlines(
                draw_labels=True, linewidth=0.4,
                color="grey", alpha=0.5, linestyle=":",
            )
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": FS_TICK}
            gl.ylabel_style = {"size": FS_TICK}

        # Interior gridlines without labels
        ax.gridlines(
            draw_labels=False, linewidth=0.3,
            color="grey", alpha=0.35, linestyle=":",
            zorder=2,
        )

        title_text = f"{site['panel_label']} {site['name']}"
        ax.set_title(
            title_text, fontsize=FS_PANEL_TITLE,
            fontweight="bold", loc="left", pad=6,
        )

    # Shared colourbar
    cbar_ax = fig.add_axes([0.25, 0.09, 0.50, 0.020])
    cbar = fig.colorbar(
        mesh_handle, cax=cbar_ax, orientation="horizontal",
        extend="both",
    )
    cbar.set_label(
        "Sea surface temperature trend (°C per decade)",
        fontsize=FS_CBAR_LABEL, labelpad=5,
    )
    cbar.ax.tick_params(labelsize=FS_CBAR_TICK)

    # Shared legend
    legend_handles = []
    for key in PANEL_ORDER:
        site = SITES[key]
        legend_handles.append(
            mpatches.Patch(
                facecolor=site["colour"], edgecolor=site["colour"],
                linewidth=3.2, alpha=ENVELOPE_FILL_ALPHA + 0.15,
                label=f"{site['name']} meander envelope",
            )
        )
    legend_handles.append(
        plt.Line2D(
            [0], [0], color="#222222", linewidth=1.2, linestyle="--",
            label="Eddy kinetic energy trend contours (inside envelope)",
        )
    )

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.50, 0.005),
        ncol=3, fontsize=FS_LEGEND,
        framealpha=0.95, frameon=True,
        columnspacing=1.6, handletextpad=0.8,
    )

    out_pdf = FIG_DIR / "fig2_meander_site_sst_trends.pdf"
    out_png = FIG_DIR / "fig2_meander_site_sst_trends.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    ds_sst.close()

    print(f"  -> {out_pdf}")
    print(f"  -> {out_png}")
    print("Figure 2 v5 done.\n")


if __name__ == "__main__":
    make_figure2()
    print("=" * 70)
    print("NB15 Figure 2 v5 generation complete.")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 70)
