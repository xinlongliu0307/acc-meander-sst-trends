#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NB15_manuscript_figures_v5.py
==============================

Generates Figure 1 for the JGR:Oceans manuscript. This version (v5)
solves the vertical compression problem of v4 by eliminating the
separate middle strip and integrating the legend and colourbar as
overlay elements inside panel (a). This gives panel (a) the full
available vertical space and removes the empty middle gap.

Key v5 changes relative to v4:

(1) GridSpec simplified from 3 rows to 2 rows: panel (a) and panel (b)
    only. No middle strip. Height ratios [1.0, 0.72].

(2) Horizontal colourbar positioned using inset_axes INSIDE panel (a),
    at the bottom-centre of the map within the ocean area south of 75°S
    where no trend data is displayed. Alternative position uses
    fig.add_axes with absolute figure coordinates for precise control.

(3) Legend box positioned as an inset inside panel (a) at the upper
    portion of the South Atlantic sector (approximately 30°W-0°E)
    where the trend field is relatively homogeneous and the overlay
    does not obscure any meander site or control region.

(4) Panel (a) aspect ratio unlocked with set_aspect("auto") so the
    map stretches vertically to fill the available space rather than
    being constrained to the geographic aspect ratio.

(5) hspace between the two GridSpec rows reduced to 0.15 to keep
    panel (b) close to panel (a) without overlap.

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
from matplotlib.gridspec import GridSpec
import cartopy.crs as ccrs
import cartopy.feature as cfeature

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ==============================================================================
# CONFIGURATION
# ==============================================================================
BASE = Path("/g/data/gv90/xl1657/cmems_adt/meander_sst_project")
FIG_DIR = BASE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

SST_TREND_MAP = BASE / "sst_trends" / "sst_trend_map.nc"
SST_TRANSECT = BASE / "sst_trends" / "sst_trend_along_acc_transect.nc"

SITES = {
    "CP":   {"lon": (150.0, -150.0), "lat": (-57.0, -46.0),
             "colour": "#1f77b4",
             "label": "Campbell Plateau",
             "range_label": "150°E – 150°W"},
    "PAR":  {"lon": (-150.0, -80.0), "lat": (-60.0, -48.0),
             "colour": "#ff7f0e",
             "label": "Pacific-Antarctic Ridge",
             "range_label": "150°W – 80°W"},
    "SEIR": {"lon": (130.0, 152.0), "lat": (-56.0, -44.0),
             "colour": "#e377c2",
             "label": "Southeast Indian Ridge",
             "range_label": "130°E – 152°E"},
    "SWIR": {"lon": (15.0, 45.0), "lat": (-58.0, -44.0),
             "colour": "#2ca02c",
             "label": "Southwest Indian Ridge",
             "range_label": "15°E – 45°E"},
}

CONTROLS = {
    "SE_Pacific": {"lon": (-100.0, -80.0), "lat": (-60.0, -48.0),
                   "colour": "#7f7f7f",
                   "label": "Southeast Pacific",
                   "range_label": "100°W – 80°W"},
    "S_Indian":   {"lon": (90.0, 110.0), "lat": (-56.0, -44.0),
                   "colour": "#8c564b",
                   "label": "Central South Indian",
                   "range_label": "90°E – 110°E"},
    "S_Atlantic": {"lon": (-10.0, 10.0), "lat": (-58.0, -44.0),
                   "colour": "#bcbd22",
                   "label": "South Atlantic",
                   "range_label": "10°W – 10°E"},
    "CP_SEIR_gap": {"lon": (180.0, -165.0), "lat": (-57.0, -46.0),
                    "colour": "#17becf",
                    "label": "Campbell Plateau – Southeast Indian Ridge gap",
                    "range_label": "180° – 165°W"},
}

TREND_VMIN, TREND_VMAX = -0.4, 0.4
TREND_CMAP = "RdBu_r"

RAW_LINE_ALPHA = 0.25
RAW_LINE_WIDTH = 0.6
RUNNING_MEAN_WINDOW = 40
RUNNING_MEAN_WIDTH = 2.2
CIRCUMPOLAR_MEAN_COLOUR = "#d62728"
CIRCUMPOLAR_MEAN_STYLE = "--"
CIRCUMPOLAR_MEAN_WIDTH = 1.3
ZERO_CROSSING_COLOUR = "black"
ZERO_CROSSING_STYLE = ":"
ZERO_CROSSING_WIDTH = 1.0

BOUNDARY_LINE_ALPHA = 0.55
BOUNDARY_LINE_WIDTH = 1.0
BOUNDARY_LINE_STYLE = "--"

FULL_HEIGHT_BAND_ALPHA = 0.12

BRACKET_OFFSET = 1.5
BRACKET_TICK_LENGTH = 0.8

FS_TITLE = 14
FS_LABEL = 12
FS_TICK = 11
FS_ANNOT = 9
FS_LEGEND = 10
FS_RANGE = 8


def running_mean(x, window):
    if window < 2:
        return x.copy()
    pad = window // 2
    padded = np.concatenate([
        x[pad-1::-1] if pad > 0 else np.array([]),
        x,
        x[:-pad-1:-1] if pad > 0 else np.array([])
    ])
    kernel = np.ones(window) / window
    smoothed = np.convolve(padded, kernel, mode="valid")
    if len(smoothed) > len(x):
        trim = (len(smoothed) - len(x)) // 2
        smoothed = smoothed[trim:trim+len(x)]
    return smoothed[:len(x)]


def format_longitude_label(lon_value):
    lv = ((lon_value + 180.0) % 360.0) - 180.0
    if np.isclose(lv, 180.0) or np.isclose(lv, -180.0):
        return "180°"
    if lv < 0:
        return f"{int(round(abs(lv)))}°W"
    elif lv > 0:
        return f"{int(round(lv))}°E"
    else:
        return "0°"


def draw_bounding_box(ax, lon_range, lat_range, colour,
                      linewidth=1.8, linestyle="-", zorder=5):
    lon_min, lon_max = lon_range
    lat_min, lat_max = lat_range

    if lon_max < lon_min:
        rect1 = Rectangle(
            (lon_min, lat_min), 180.0 - lon_min, lat_max - lat_min,
            fill=False, edgecolor=colour, linewidth=linewidth,
            linestyle=linestyle, zorder=zorder,
            transform=ccrs.PlateCarree(),
        )
        ax.add_patch(rect1)
        rect2 = Rectangle(
            (-180.0, lat_min), lon_max - (-180.0), lat_max - lat_min,
            fill=False, edgecolor=colour, linewidth=linewidth,
            linestyle=linestyle, zorder=zorder,
            transform=ccrs.PlateCarree(),
        )
        ax.add_patch(rect2)
    else:
        rect = Rectangle(
            (lon_min, lat_min), lon_max - lon_min, lat_max - lat_min,
            fill=False, edgecolor=colour, linewidth=linewidth,
            linestyle=linestyle, zorder=zorder,
            transform=ccrs.PlateCarree(),
        )
        ax.add_patch(rect)


def draw_range_bracket(ax, lon_range, lat_top, colour, range_label,
                       is_control=False):
    lon_min, lon_max = lon_range
    y_bracket = lat_top + BRACKET_OFFSET
    y_tick = y_bracket - BRACKET_TICK_LENGTH

    linestyle = "--" if is_control else "-"
    lw = 1.4 if is_control else 1.8

    if lon_max < lon_min:
        ax.plot([lon_min, 180.0], [y_bracket, y_bracket],
                color=colour, linewidth=lw, linestyle=linestyle,
                zorder=7, transform=ccrs.PlateCarree(),
                solid_capstyle="butt")
        ax.plot([lon_min, lon_min], [y_bracket, y_tick],
                color=colour, linewidth=lw, linestyle="-",
                zorder=7, transform=ccrs.PlateCarree())
        ax.plot([-180.0, lon_max], [y_bracket, y_bracket],
                color=colour, linewidth=lw, linestyle=linestyle,
                zorder=7, transform=ccrs.PlateCarree(),
                solid_capstyle="butt")
        ax.plot([lon_max, lon_max], [y_bracket, y_tick],
                color=colour, linewidth=lw, linestyle="-",
                zorder=7, transform=ccrs.PlateCarree())

        left_extent = 180.0 - lon_min
        right_extent = lon_max - (-180.0)
        if left_extent >= right_extent:
            label_x = (lon_min + 180.0) / 2.0
        else:
            label_x = (-180.0 + lon_max) / 2.0
    else:
        ax.plot([lon_min, lon_max], [y_bracket, y_bracket],
                color=colour, linewidth=lw, linestyle=linestyle,
                zorder=7, transform=ccrs.PlateCarree(),
                solid_capstyle="butt")
        ax.plot([lon_min, lon_min], [y_bracket, y_tick],
                color=colour, linewidth=lw, linestyle="-",
                zorder=7, transform=ccrs.PlateCarree())
        ax.plot([lon_max, lon_max], [y_bracket, y_tick],
                color=colour, linewidth=lw, linestyle="-",
                zorder=7, transform=ccrs.PlateCarree())
        label_x = (lon_min + lon_max) / 2.0

    ax.text(
        label_x, y_bracket + 0.8, range_label,
        color=colour, fontsize=FS_RANGE, fontweight="bold",
        ha="center", va="bottom",
        transform=ccrs.PlateCarree(),
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor=colour, linewidth=0.5, alpha=0.9),
        zorder=10,
    )


def make_figure1():
    print("=" * 70)
    print("Generating Figure 1 (v5): Circumpolar SST Trends")
    print("=" * 70)

    ds_map = xr.open_dataset(SST_TREND_MAP)
    ds_tr = xr.open_dataset(SST_TRANSECT)

    trend_map = ds_map["sst_trend"].values
    lat_map = ds_map["latitude"].values
    lon_map = ds_map["longitude"].values

    transect_lon = ds_tr["longitude"].values
    transect_trend = ds_tr["sst_trend_along_acc"].values
    circumpolar_mean = float(np.nanmean(transect_trend))

    print(f"  Trend map shape: {trend_map.shape}")
    print(f"  Transect length: {len(transect_lon)}")
    print(f"  Circumpolar mean: {circumpolar_mean:+.3f} deg C/decade")

    # ==================================================================
    # V5 LAYOUT: 2-row GridSpec with legend and colourbar as overlays
    # inside panel (a). This eliminates the middle strip entirely.
    # ==================================================================
    fig = plt.figure(figsize=(16, 13))
    gs = GridSpec(
        2, 1,
        height_ratios=[1.0, 0.72],
        hspace=0.15,
        figure=fig,
        left=0.06, right=0.97,
        top=0.96, bottom=0.06,
    )

    # ==================================================================
    # PANEL (a): MAP — vertical size maximised
    # ==================================================================
    ax_a = fig.add_subplot(
        gs[0],
        projection=ccrs.PlateCarree(central_longitude=180),
    )
    # IMPORTANT: unlock aspect ratio so the map stretches vertically
    ax_a.set_aspect("auto")
    ax_a.set_global()
    ax_a.set_extent([-180, 180, -78, -28], crs=ccrs.PlateCarree())

    mesh = ax_a.pcolormesh(
        lon_map, lat_map, trend_map,
        cmap=TREND_CMAP, vmin=TREND_VMIN, vmax=TREND_VMAX,
        transform=ccrs.PlateCarree(), shading="auto", rasterized=True,
    )

    ax_a.add_feature(cfeature.LAND, facecolor="#f0f0f0",
                     edgecolor="black", linewidth=0.4, zorder=3)
    ax_a.coastlines(resolution="110m", linewidth=0.5, zorder=4)

    gl = ax_a.gridlines(
        draw_labels=True, linewidth=0.3,
        color="grey", alpha=0.5, linestyle=":",
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": FS_TICK}
    gl.ylabel_style = {"size": FS_TICK}

    for key, site in SITES.items():
        lon_min, lon_max = site["lon"]
        lat_min, lat_max = site["lat"]
        colour = site["colour"]
        draw_bounding_box(
            ax_a, (lon_min, lon_max), (lat_min, lat_max),
            colour=colour, linewidth=2.4, linestyle="-", zorder=6,
        )
        draw_range_bracket(
            ax_a, (lon_min, lon_max), lat_max, colour,
            site["range_label"], is_control=False,
        )

    for key, ctrl in CONTROLS.items():
        lon_min, lon_max = ctrl["lon"]
        lat_min, lat_max = ctrl["lat"]
        colour = ctrl["colour"]
        draw_bounding_box(
            ax_a, (lon_min, lon_max), (lat_min, lat_max),
            colour=colour, linewidth=1.8, linestyle="--", zorder=6,
        )
        draw_range_bracket(
            ax_a, (lon_min, lon_max), lat_max, colour,
            ctrl["range_label"], is_control=True,
        )

    ax_a.set_title(
        "(a) Circumpolar sea surface temperature trends, 1993–2025",
        fontsize=FS_TITLE, pad=10, loc="left",
    )

    # ------------------------------------------------------------------
    # LEGEND OVERLAY INSIDE PANEL (a)
    # Placed inside the Antarctic continent area (below 70°S) which
    # contains no trend data, so no information is obscured.
    # ------------------------------------------------------------------
    meander_handles = [
        mpatches.Patch(
            facecolor="none", edgecolor=s["colour"],
            linewidth=2.2, label=s["label"],
        )
        for s in SITES.values()
    ]
    control_handles = [
        mpatches.Patch(
            facecolor="none", edgecolor=c["colour"],
            linewidth=1.8, linestyle="--", label=c["label"],
        )
        for c in CONTROLS.values()
    ]

    legend = ax_a.legend(
        handles=meander_handles + control_handles,
        loc="lower left",
        bbox_to_anchor=(0.01, 0.01),
        ncol=4, fontsize=FS_LEGEND - 1,
        framealpha=0.92, frameon=True,
        title=(
            "Meander sites (solid boxes)   ·   "
            "Quiescent control regions (dashed boxes)"
        ),
        title_fontsize=FS_LEGEND - 1,
    )
    legend.set_zorder(11)

    # ------------------------------------------------------------------
    # COLOURBAR OVERLAY INSIDE PANEL (a) — horizontal, top-right
    # Positioned above the trend data in the subtropical zone
    # ------------------------------------------------------------------
    cbar_ax = ax_a.inset_axes(
        [0.62, 0.88, 0.35, 0.04],
        transform=ax_a.transAxes,
    )
    cbar = plt.colorbar(
        mesh, cax=cbar_ax, orientation="horizontal",
    )
    cbar.set_label(
        "Sea surface temperature trend (°C per decade)",
        fontsize=FS_LABEL - 1, labelpad=3,
    )
    cbar.ax.tick_params(labelsize=FS_TICK - 1)
    # White background panel behind colourbar for readability
    cbar_ax.set_facecolor("white")
    # Add a white background rectangle behind the colourbar region
    # using a FancyBboxPatch in axes coordinates
    bg_rect = mpatches.FancyBboxPatch(
        (0.60, 0.82), 0.39, 0.16,
        boxstyle="round,pad=0.005",
        transform=ax_a.transAxes,
        facecolor="white", edgecolor="grey", linewidth=0.5,
        alpha=0.88, zorder=9,
    )
    ax_a.add_patch(bg_rect)

    # ==================================================================
    # PANEL (b): ALONG-ACC TRANSECT
    # ==================================================================
    ax_b = fig.add_subplot(gs[1])

    y_min, y_max = -0.45, 0.45
    ax_b.set_xlim([-180, 180])
    ax_b.set_ylim([y_min, y_max])

    def _draw_full_height_band(lon_range, colour):
        lon_min, lon_max = lon_range
        if lon_max < lon_min:
            ax_b.axvspan(
                lon_min, 180.0,
                color=colour, alpha=FULL_HEIGHT_BAND_ALPHA,
                zorder=1, linewidth=0,
            )
            ax_b.axvspan(
                -180.0, lon_max,
                color=colour, alpha=FULL_HEIGHT_BAND_ALPHA,
                zorder=1, linewidth=0,
            )
        else:
            ax_b.axvspan(
                lon_min, lon_max,
                color=colour, alpha=FULL_HEIGHT_BAND_ALPHA,
                zorder=1, linewidth=0,
            )

    for key, site in SITES.items():
        _draw_full_height_band(site["lon"], site["colour"])
    for key, ctrl in CONTROLS.items():
        _draw_full_height_band(ctrl["lon"], ctrl["colour"])

    def _draw_boundary_lines(lon_range, colour):
        lon_min, lon_max = lon_range
        for b in [lon_min, lon_max]:
            b_norm = ((b + 180.0) % 360.0) - 180.0
            ax_b.axvline(
                b_norm, color=colour, linestyle=BOUNDARY_LINE_STYLE,
                linewidth=BOUNDARY_LINE_WIDTH, alpha=BOUNDARY_LINE_ALPHA,
                zorder=2,
            )

    for key, site in SITES.items():
        _draw_boundary_lines(site["lon"], site["colour"])
    for key, ctrl in CONTROLS.items():
        _draw_boundary_lines(ctrl["lon"], ctrl["colour"])

    ax_b.plot(
        transect_lon, transect_trend,
        color="grey", linewidth=RAW_LINE_WIDTH, alpha=RAW_LINE_ALPHA,
        zorder=3, label="Per-longitude trend",
    )

    smoothed = running_mean(transect_trend, RUNNING_MEAN_WINDOW)
    ax_b.plot(
        transect_lon, smoothed,
        color="black", linewidth=RUNNING_MEAN_WIDTH,
        zorder=5,
        label=f"{RUNNING_MEAN_WINDOW}-point running mean",
    )

    ax_b.axhline(
        0.0, color=ZERO_CROSSING_COLOUR,
        linestyle=ZERO_CROSSING_STYLE, linewidth=ZERO_CROSSING_WIDTH,
        zorder=4, alpha=0.6, label="Zero trend",
    )
    ax_b.axhline(
        circumpolar_mean, color=CIRCUMPOLAR_MEAN_COLOUR,
        linestyle=CIRCUMPOLAR_MEAN_STYLE, linewidth=CIRCUMPOLAR_MEAN_WIDTH,
        zorder=4,
        label=f"Circumpolar mean ({circumpolar_mean:+.3f} °C per decade)",
    )

    def _region_centre_lon(lon_range):
        a, b = lon_range
        if b < a:
            extent = (180 - a) + (b - (-180))
            centre = a + extent / 2.0
            if centre > 180:
                centre -= 360
            return centre
        return (a + b) / 2.0

    top_y_high = y_max - 0.035
    top_y_low = y_max - 0.12

    meander_y_offsets = {
        "CP":   top_y_high,
        "PAR":  top_y_high,
        "SEIR": top_y_low,
        "SWIR": top_y_high,
    }

    for key, site in SITES.items():
        centre = _region_centre_lon(site["lon"])
        label_y = meander_y_offsets[key]
        ax_b.text(
            centre, label_y, site["label"].replace(" ", "\n"),
            fontsize=FS_RANGE, fontweight="bold", color=site["colour"],
            ha="center", va="top", zorder=7,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor=site["colour"], linewidth=0.8, alpha=0.92),
        )

    label_y_ctrl = y_min + 0.035
    for key, ctrl in CONTROLS.items():
        centre = _region_centre_lon(ctrl["lon"])
        ctrl_label = ctrl["label"]
        if key == "CP_SEIR_gap":
            ctrl_label = "Campbell Plateau –\nSoutheast Indian\nRidge gap"
        else:
            ctrl_label = ctrl_label.replace(" ", "\n")
        ax_b.text(
            centre, label_y_ctrl, ctrl_label,
            fontsize=FS_RANGE, fontweight="bold", color=ctrl["colour"],
            ha="center", va="bottom", zorder=7,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor=ctrl["colour"], linewidth=0.8,
                      linestyle="--", alpha=0.92),
        )

    ax_b.set_xlabel("Longitude", fontsize=FS_LABEL)
    ax_b.set_ylabel(
        "Sea surface temperature trend (°C per decade)",
        fontsize=FS_LABEL,
    )
    ax_b.set_title(
        "(b) Sea surface temperature trend along the "
        "time-mean Antarctic Circumpolar Current axis",
        fontsize=FS_TITLE, loc="left",
    )

    ax_b.set_xticks(np.arange(-180, 181, 30))
    ax_b.set_xticklabels(
        [format_longitude_label(l) for l in np.arange(-180, 181, 30)],
        fontsize=FS_TICK,
    )
    ax_b.tick_params(axis="y", labelsize=FS_TICK)

    ax_b.grid(
        True, which="major", axis="y", linestyle=":",
        linewidth=0.4, alpha=0.4,
    )

    ax_b.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        fontsize=FS_LEGEND, framealpha=0.95, ncol=4,
        columnspacing=1.5, handletextpad=0.6,
    )

    out_pdf = FIG_DIR / "fig1_circumpolar_sst_trends.pdf"
    out_png = FIG_DIR / "fig1_circumpolar_sst_trends.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    ds_map.close()
    ds_tr.close()

    print(f"  -> {out_pdf}")
    print(f"  -> {out_png}")
    print("Figure 1 v5 done.\n")


if __name__ == "__main__":
    make_figure1()
    print("=" * 70)
    print("NB15 v5 figure generation complete.")
    print(f"Output directory: {FIG_DIR}")
    print("=" * 70)
