#!/usr/bin/env python3
"""
setup_meander_sst_project.py
=============================
Task C: Create directory structure and meander envelope masks for the
new GRL manuscript on meander-modulated SST trends.

Generates:
  1. Project directory tree under meander_sst_project/
  2. Meander envelope masks from circumpolar_belt_rel20_x3m.nc
  3. Control (quiescent) region masks for null comparison
  4. Site-specific bounding-box masks on the 0.125° ADT grid

Run on NCI Gadi:
  cd /g/data/gv90/xl1657/cmems_adt/notebooks/
  source /g/data/gv90/xl1657/venvs/cmems_py311/bin/activate
  python setup_meander_sst_project.py

Author: Xinlong Liu, IMAS, University of Tasmania
"""

from pathlib import Path
import numpy as np
import xarray as xr
import json

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
BASE_DIR = Path("/g/data/gv90/xl1657/cmems_adt")
BELT_FP  = BASE_DIR / "circumpolar_belt_rel20_x3m.nc"
PROJ_DIR = BASE_DIR / "meander_sst_project"

# ADT grid parameters (0.125°, 90S–30S)
LAT_GRID = np.arange(-89.9375, -29.9, 0.125)   # 480 points
LON_GRID = np.arange(-179.9375, 180.0, 0.125)   # 2880 points

# Site bounding boxes — must match NB02/NB03 exactly
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

# Control (quiescent ACC) regions — away from major topographic features
# Each spans roughly the same latitude band as the ACC (~58S to ~45S)
CONTROLS = {
    "CTRL_SE_PAC": {"name": "SE Pacific (quiescent)",
                    "lon": (-100, -80), "lat": (-60, -48), "wraps": False},
    "CTRL_S_IND":  {"name": "Central S Indian (quiescent)",
                    "lon": (90, 110),   "lat": (-56, -44), "wraps": False},
    "CTRL_S_ATL":  {"name": "S Atlantic (quiescent)",
                    "lon": (10, -10),   "lat": (-58, -44), "wraps": True},
    "CTRL_CP_SEIR":{"name": "Between CP and SEIR (quiescent)",
                    "lon": (-180, -165),"lat": (-57, -46), "wraps": False},
}


# =============================================================================
# 1. CREATE DIRECTORY TREE
# =============================================================================
def create_directories():
    """Create the project directory structure."""
    dirs = [
        PROJ_DIR,
        PROJ_DIR / "masks",
        PROJ_DIR / "sst_trends",
        PROJ_DIR / "eke_spatial",
        PROJ_DIR / "joint_analysis",
        PROJ_DIR / "figures",
        PROJ_DIR / "si_figures",
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        print(f"  Created: {d}")


# =============================================================================
# 2. GENERATE MEANDER ENVELOPE MASKS FROM CIRCUMPOLAR BELT
# =============================================================================
def generate_meander_envelope_masks():
    """
    From the circumpolar belt file, extract the time-mean meander envelope
    (northernmost and southernmost monthly positions across 1993–2025)
    and create a binary mask on the ADT grid for each site.
    """
    print("\nGenerating meander envelope masks...")

    if not BELT_FP.exists():
        print(f"  ERROR: Circumpolar belt file not found: {BELT_FP}")
        print("  Falling back to site bounding-box masks only.")
        return None

    ds = xr.open_dataset(BELT_FP)
    print(f"  Loaded: {BELT_FP.name}")
    print(f"  Dimensions: {dict(ds.dims)}")
    print(f"  Variables: {list(ds.data_vars)}")

    belt_lon = ds["longitude"].values
    center_lat = ds["meander_center_lat"].values    # (month, lon)
    north_lat  = ds["meander_north_lat"].values     # (month, lon)
    south_lat  = ds["meander_south_lat"].values     # (month, lon)
    n_months = center_lat.shape[0]

    # Compute the envelope: across ALL months, find the most extreme
    # northward and southward positions at each longitude
    envelope_north = np.nanmax(north_lat, axis=0)   # (lon,)
    envelope_south = np.nanmin(south_lat, axis=0)   # (lon,)
    envelope_center = np.nanmean(center_lat, axis=0) # (lon,)

    # Also compute the time-mean meander center for along-ACC transects
    mean_center = np.nanmean(center_lat, axis=0)    # (lon,)

    # Build the full-domain envelope dataset
    ds_envelope = xr.Dataset(
        {
            "envelope_north_lat": (["longitude"], envelope_north),
            "envelope_south_lat": (["longitude"], envelope_south),
            "envelope_center_lat": (["longitude"], envelope_center),
            "mean_center_lat": (["longitude"], mean_center),
        },
        coords={"longitude": belt_lon},
        attrs={
            "title": "Circumpolar meander envelope from monthly positions",
            "source": str(BELT_FP),
            "description": "North/south extremes across all months define the envelope",
        },
    )

    fp_env = PROJ_DIR / "masks" / "circumpolar_meander_envelope.nc"
    ds_envelope.to_netcdf(fp_env)
    print(f"  Saved circumpolar envelope: {fp_env}")

    # Now create per-site binary masks on the ADT grid
    for site_key, site in SITES.items():
        mask = _make_site_envelope_mask(
            belt_lon, envelope_north, envelope_south,
            site["inner_lon"], site["inner_lat"], site["wraps"],
        )
        print(f"  {site_key}: {np.sum(mask)} grid points in meander envelope")

        ds_mask = xr.Dataset(
            {
                "meander_mask": (["latitude", "longitude"], mask.astype(np.int8)),
            },
            coords={
                "latitude": LAT_GRID,
                "longitude": LON_GRID,
            },
            attrs={
                "site": site_key,
                "site_name": site["name"],
                "description": "1 = inside meander envelope, 0 = outside",
                "inner_lon": str(site["inner_lon"]),
                "inner_lat": str(site["inner_lat"]),
            },
        )
        fp = PROJ_DIR / "masks" / f"meander_envelope_mask_{site_key}.nc"
        ds_mask.to_netcdf(fp)
        print(f"    Saved: {fp}")

    ds.close()
    return ds_envelope


def _make_site_envelope_mask(belt_lon, env_north, env_south,
                             inner_lon, inner_lat, wraps):
    """
    Create a 2D binary mask (nlat, nlon) on the ADT grid where
    the grid point falls within the meander envelope AND the site
    bounding box.
    """
    nlat = len(LAT_GRID)
    nlon = len(LON_GRID)
    mask = np.zeros((nlat, nlon), dtype=bool)

    # Longitude selection on ADT grid
    if not wraps:
        lon_sel = (LON_GRID >= inner_lon[0]) & (LON_GRID <= inner_lon[1])
    else:
        lon_sel = (LON_GRID >= inner_lon[0]) | (LON_GRID <= inner_lon[1])

    # Latitude bounding box
    lat_sel = (LAT_GRID >= inner_lat[0]) & (LAT_GRID <= inner_lat[1])

    for j in range(nlon):
        if not lon_sel[j]:
            continue

        # Find nearest longitude in the belt file
        this_lon = LON_GRID[j]
        lon_diffs = np.abs(belt_lon - this_lon)
        # Handle wrapping for CP
        if wraps:
            lon_diffs = np.minimum(lon_diffs, 360 - lon_diffs)
        nearest_idx = np.argmin(lon_diffs)

        n_lat_val = env_north[nearest_idx]
        s_lat_val = env_south[nearest_idx]

        if np.isnan(n_lat_val) or np.isnan(s_lat_val):
            continue

        for i in range(nlat):
            if not lat_sel[i]:
                continue
            if s_lat_val <= LAT_GRID[i] <= n_lat_val:
                mask[i, j] = True

    return mask


# =============================================================================
# 3. GENERATE CONTROL REGION MASKS
# =============================================================================
def generate_control_masks():
    """Create binary masks for quiescent ACC control regions."""
    print("\nGenerating control region masks...")

    for ctrl_key, ctrl in CONTROLS.items():
        mask = np.zeros((len(LAT_GRID), len(LON_GRID)), dtype=np.int8)
        lat_sel = (LAT_GRID >= ctrl["lat"][0]) & (LAT_GRID <= ctrl["lat"][1])

        if not ctrl.get("wraps", False):
            lon_sel = (LON_GRID >= ctrl["lon"][0]) & (LON_GRID <= ctrl["lon"][1])
        else:
            lon_sel = (LON_GRID >= ctrl["lon"][0]) | (LON_GRID <= ctrl["lon"][1])

        mask[np.ix_(lat_sel, lon_sel)] = 1

        ds_ctrl = xr.Dataset(
            {"control_mask": (["latitude", "longitude"], mask)},
            coords={"latitude": LAT_GRID, "longitude": LON_GRID},
            attrs={
                "region": ctrl_key,
                "region_name": ctrl["name"],
                "lon_range": str(ctrl["lon"]),
                "lat_range": str(ctrl["lat"]),
            },
        )
        fp = PROJ_DIR / "masks" / f"control_mask_{ctrl_key}.nc"
        ds_ctrl.to_netcdf(fp)
        n_pts = int(mask.sum())
        print(f"  {ctrl_key}: {n_pts} grid points -> {fp}")


# =============================================================================
# 4. SAVE CONFIGURATION METADATA
# =============================================================================
def save_config():
    """Save project configuration as JSON for reproducibility."""
    config = {
        "project": "GRL Meander SST Trends",
        "base_dir": str(BASE_DIR),
        "project_dir": str(PROJ_DIR),
        "belt_file": str(BELT_FP),
        "adt_grid": {"lat_range": [-89.9375, -30.0625], "lon_range": [-179.9375, 179.9375],
                     "spacing_deg": 0.125, "nlat": 480, "nlon": 2880},
        "sites": {k: {"name": v["name"], "inner_lon": v["inner_lon"],
                       "inner_lat": v["inner_lat"], "wraps": v["wraps"]}
                  for k, v in SITES.items()},
        "controls": {k: {"name": v["name"], "lon": v["lon"],
                         "lat": v["lat"], "wraps": v.get("wraps", False)}
                     for k, v in CONTROLS.items()},
    }
    fp = PROJ_DIR / "project_config.json"
    with open(fp, "w") as f:
        json.dump(config, f, indent=2)
    print(f"\nSaved config: {fp}")


# =============================================================================
# 5. MAIN
# =============================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("Meander SST Project Setup")
    print("=" * 70)

    create_directories()
    generate_meander_envelope_masks()
    generate_control_masks()
    save_config()

    print("\n" + "=" * 70)
    print("Setup complete. Directory structure:")
    print(f"  {PROJ_DIR}/")
    print(f"    masks/           <- meander envelope + control masks")
    print(f"    sst_trends/      <- SST trend fields (Phase 1)")
    print(f"    eke_spatial/     <- spatial EKE fields from NB03b")
    print(f"    joint_analysis/  <- meander vs SST joint outputs (Phase 2)")
    print(f"    figures/         <- manuscript figures (Phase 4)")
    print(f"    si_figures/      <- Supporting Information figures")
    print("=" * 70)
