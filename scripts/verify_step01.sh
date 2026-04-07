#!/bin/bash
# ==============================================================================
# Step 0.1 — Verify Existing Derived Datasets on NCI Gadi
# Project: GRL Manuscript — Meander SST Trends
# Author: Xinlong Liu
# Environment: NCI Gadi, project gv90, user xl1657
#
# Usage:
#   chmod +x verify_step01.sh
#   ./verify_step01.sh          # Run Stage A (discovery) + Stage B (integrity)
#   ./verify_step01.sh --scan   # Run Stage A only (discovery scan)
#
# This script does NOT modify any files. It is read-only and diagnostic.
# ==============================================================================

set -euo pipefail

# --- Configuration -----------------------------------------------------------
BASE_DIR="/g/data/gv90/xl1657/cmems_adt"
VENV_DIR="/g/data/gv90/xl1657/venvs/cmems_py311"
LOG_FILE="verify_step01_$(date +%Y%m%d_%H%M%S).log"
SITES=("CP" "PAR" "SEIR" "SWIR")

# Colour codes for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No colour

# --- Helper functions --------------------------------------------------------
log() {
    echo -e "$1" | tee -a "$LOG_FILE"
}

check_exists() {
    local path="$1"
    local label="$2"
    if [ -e "$path" ]; then
        log "  ${GREEN}[FOUND]${NC}  $label"
        log "          -> $path"
        if [ -f "$path" ]; then
            local size
            size=$(du -h "$path" | cut -f1)
            local mtime
            mtime=$(stat -c '%y' "$path" 2>/dev/null | cut -d'.' -f1)
            log "          -> Size: $size | Last modified: $mtime"
        fi
        return 0
    else
        log "  ${RED}[MISSING]${NC} $label"
        log "          -> Expected: $path"
        return 1
    fi
}

check_dir_exists() {
    local path="$1"
    local label="$2"
    if [ -d "$path" ]; then
        local nfiles
        nfiles=$(find "$path" -maxdepth 1 -type f | wc -l)
        local ndirs
        ndirs=$(find "$path" -maxdepth 1 -type d | wc -l)
        ndirs=$((ndirs - 1))  # exclude self
        log "  ${GREEN}[FOUND]${NC}  $label"
        log "          -> $path ($nfiles files, $ndirs subdirs)"
        return 0
    else
        log "  ${RED}[MISSING]${NC} $label"
        log "          -> Expected: $path"
        return 1
    fi
}

section_header() {
    log ""
    log "=============================================================================="
    log "${CYAN}$1${NC}"
    log "=============================================================================="
}

# ==============================================================================
# STAGE A: Discovery Scan
# Purpose: Map the full directory tree so we know what exists and where.
# ==============================================================================
stage_a() {
    section_header "STAGE A: Directory Discovery Scan"
    log "Timestamp: $(date)"
    log "Base directory: $BASE_DIR"
    log ""

    # A.1 — Check base directory exists
    log "--- A.1: Base Directory ---"
    if ! check_dir_exists "$BASE_DIR" "Base data directory"; then
        log "${RED}CRITICAL: Base directory does not exist. Aborting.${NC}"
        exit 1
    fi

    # A.2 — Full directory tree (2 levels deep)
    log ""
    log "--- A.2: Directory Tree (3 levels deep) ---"
    log ""
    tree -L 3 -d "$BASE_DIR" 2>/dev/null >> "$LOG_FILE" || \
        find "$BASE_DIR" -maxdepth 3 -type d | sort >> "$LOG_FILE"
    log "  (see $LOG_FILE for full tree)"

    # A.3 — List all file types and counts
    log ""
    log "--- A.3: File Type Summary ---"
    log ""
    log "  Extension | Count | Total Size"
    log "  ----------|-------|----------"
    # Find all file extensions and count them
    find "$BASE_DIR" -type f | sed 's/.*\.//' | sort | uniq -c | sort -rn | \
    while read -r count ext; do
        total_size=$(find "$BASE_DIR" -type f -name "*.$ext" -exec du -ch {} + 2>/dev/null | tail -1 | cut -f1)
        printf "  %-10s | %5d | %s\n" ".$ext" "$count" "$total_size" | tee -a "$LOG_FILE"
    done

    # A.4 — Search for site-specific files
    log ""
    log "--- A.4: Site-Specific File Discovery ---"
    for site in "${SITES[@]}"; do
        log ""
        log "  ${YELLOW}Site: $site${NC}"
        # Case-insensitive search for files containing site name
        local found
        found=$(find "$BASE_DIR" -type f \( -iname "*${site}*" -o -iname "*$(echo "$site" | tr '[:upper:]' '[:lower:]')*" \) 2>/dev/null | head -20)
        if [ -n "$found" ]; then
            echo "$found" | while read -r f; do
                local size
                size=$(du -h "$f" | cut -f1)
                log "    $size  $f"
            done
            local total
            total=$(find "$BASE_DIR" -type f \( -iname "*${site}*" -o -iname "*$(echo "$site" | tr '[:upper:]' '[:lower:]')*" \) 2>/dev/null | wc -l)
            if [ "$total" -gt 20 ]; then
                log "    ... and $((total - 20)) more files"
            fi
        else
            log "    ${RED}No files found containing '$site' in name${NC}"
        fi
    done

    # A.5 — Search for key pipeline scripts and outputs
    log ""
    log "--- A.5: Key Pipeline Components ---"
    log ""

    local search_terms=("NB01" "NB02" "NB03" "NB04" "NB05" "NB06" "NB07" "NB08" "NB09" "NB10"
                        "patch_width" "combine_cp" "corrected"
                        "meander_pos" "position" "eke" "width" "speed"
                        "ERA5" "era5" "wind" "trend" "mask" "section"
                        "threshold" "sensitivity")

    for term in "${search_terms[@]}"; do
        local matches
        matches=$(find "$BASE_DIR" -type f -iname "*${term}*" 2>/dev/null | head -5)
        if [ -n "$matches" ]; then
            log "  ${GREEN}[$term]${NC}"
            echo "$matches" | while read -r f; do
                log "    $(du -h "$f" | cut -f1)  $f"
            done
            local total
            total=$(find "$BASE_DIR" -type f -iname "*${term}*" 2>/dev/null | wc -l)
            if [ "$total" -gt 5 ]; then
                log "    ... ($total total matches)"
            fi
        fi
    done

    # A.6 — Search for NetCDF files and inspect headers of the first few
    log ""
    log "--- A.6: NetCDF File Inventory ---"
    local nc_files
    nc_files=$(find "$BASE_DIR" -type f \( -name "*.nc" -o -name "*.nc4" -o -name "*.netcdf" \) 2>/dev/null | sort)
    local nc_count
    nc_count=$(echo "$nc_files" | grep -c '.' 2>/dev/null || echo 0)
    log "  Total NetCDF files found: $nc_count"
    log ""

    if [ "$nc_count" -gt 0 ]; then
        log "  First 10 NetCDF files (with ncdump -h summary):"
        echo "$nc_files" | head -10 | while read -r ncf; do
            log ""
            log "  ${YELLOW}$ncf${NC}"
            log "  Size: $(du -h "$ncf" | cut -f1)"
            # Try ncdump -h for header info
            if command -v ncdump &> /dev/null; then
                ncdump -h "$ncf" 2>/dev/null | head -30 | while read -r line; do
                    log "    $line"
                done
                log "    ..."
            else
                log "    (ncdump not available; load module: module load netcdf)"
            fi
        done
    fi

    # A.7 — Search for MATLAB .mat files
    log ""
    log "--- A.7: MATLAB File Inventory ---"
    local mat_files
    mat_files=$(find "$BASE_DIR" -type f -name "*.mat" 2>/dev/null | sort)
    local mat_count
    mat_count=$(echo "$mat_files" | grep -c '.' 2>/dev/null || echo 0)
    log "  Total .mat files found: $mat_count"
    if [ "$mat_count" -gt 0 ]; then
        echo "$mat_files" | head -20 | while read -r mf; do
            log "    $(du -h "$mf" | cut -f1)  $mf"
        done
    fi

    # A.8 — Check Python virtual environment
    log ""
    log "--- A.8: Python Virtual Environment ---"
    check_dir_exists "$VENV_DIR" "Python 3.11 venv (cmems_py311)"
    if [ -f "$VENV_DIR/bin/python" ]; then
        local pyver
        pyver=$("$VENV_DIR/bin/python" --version 2>&1)
        log "  Python version: $pyver"
        log ""
        log "  Key packages:"
        for pkg in xarray numpy scipy matplotlib cartopy pandas netCDF4 xesmf cdo; do
            local ver
            ver=$("$VENV_DIR/bin/python" -c "import $pkg; print($pkg.__version__)" 2>/dev/null) || ver="NOT INSTALLED"
            log "    $pkg: $ver"
        done
    fi

    # A.9 — Disk usage summary
    log ""
    log "--- A.9: Disk Usage Summary ---"
    log ""
    du -h --max-depth=2 "$BASE_DIR" 2>/dev/null | sort -rh | head -20 | while read -r line; do
        log "  $line"
    done

    log ""
    log "--- A.10: Broader Search for ERA5 / Wind Data ---"
    log "  Searching outside cmems_adt/ for ERA5 or wind-related data..."
    local era5_broader
    era5_broader=$(find "/g/data/gv90/xl1657" -maxdepth 3 -type d \( -iname "*era5*" -o -iname "*wind*" -o -iname "*reanalysis*" \) 2>/dev/null)
    if [ -n "$era5_broader" ]; then
        echo "$era5_broader" | while read -r d; do
            log "  ${GREEN}[FOUND]${NC} $d"
        done
    else
        log "  ${YELLOW}No ERA5/wind directories found under /g/data/gv90/xl1657/${NC}"
    fi
}


# ==============================================================================
# STAGE B: Targeted Integrity Checks
# Purpose: Validate specific files that the working plan depends on.
#
# >>> IMPORTANT: Update the file paths below after Stage A reveals the
# >>> actual naming conventions and directory structure. Paths marked with
# >>> [UPDATE] are placeholders that need your input.
# ==============================================================================
stage_b() {
    section_header "STAGE B: Targeted Integrity Checks"
    log "Timestamp: $(date)"
    log ""
    log "${YELLOW}NOTE: Paths marked [PLACEHOLDER] require updating after Stage A.${NC}"
    log "${YELLOW}Edit this script and replace placeholders with actual paths.${NC}"
    log ""

    local pass=0
    local fail=0
    local warn=0

    # --- B.1: Raw CMEMS DUACS L4 ADT data ---
    log "--- B.1: Raw CMEMS DUACS L4 ADT Data ---"
    log "  Checking for ADT data files (0.125 deg, daily, 1993-2025)..."

    # [UPDATE] Adjust the glob pattern to match your actual ADT file structure
    # Common patterns: one file per year, one file per variable, or one big file
    local adt_dir="$BASE_DIR"  # [UPDATE] if ADT is in a subdirectory
    local adt_count
    adt_count=$(find "$adt_dir" -maxdepth 2 -type f \( -name "*adt*" -o -name "*ADT*" -o -name "*sealevel*" -o -name "*ssh*" -o -name "*SSH*" \) 2>/dev/null | wc -l)
    if [ "$adt_count" -gt 0 ]; then
        log "  ${GREEN}[OK]${NC} Found $adt_count ADT-related files"
        ((pass++))
    else
        log "  ${YELLOW}[WARN]${NC} No ADT files found with standard naming — check Stage A output"
        ((warn++))
    fi

    # --- B.2: Meander position fields (4 sites) ---
    log ""
    log "--- B.2: Monthly Meander Position Fields ---"
    for site in "${SITES[@]}"; do
        # [UPDATE] Replace with actual file paths after Stage A
        # Example patterns to try:
        local found=0
        for pattern in \
            "$BASE_DIR/*${site}*position*" \
            "$BASE_DIR/*${site}*pos*" \
            "$BASE_DIR/*meander*${site}*" \
            "$BASE_DIR/${site}/*position*" \
            "$BASE_DIR/${site}/*meander*" \
            "$BASE_DIR/positions/*${site}*" \
            "$BASE_DIR/*$(echo "$site" | tr '[:upper:]' '[:lower:]')*position*" \
            "$BASE_DIR/*$(echo "$site" | tr '[:upper:]' '[:lower:]')*pos*"; do
            local matches
            matches=$(ls $pattern 2>/dev/null | head -1)
            if [ -n "$matches" ]; then
                log "  ${GREEN}[FOUND]${NC} $site meander positions: $matches"
                found=1
                ((pass++))
                break
            fi
        done
        if [ "$found" -eq 0 ]; then
            log "  ${RED}[MISSING]${NC} $site meander position file — update path after Stage A"
            ((fail++))
        fi
    done

    # --- B.3: Per-grid-point EKE anomaly fields (corrected pipeline) ---
    log ""
    log "--- B.3: Corrected EKE Anomaly Fields ---"
    # Search for NB03 corrected outputs and general EKE files
    for pattern in \
        "$BASE_DIR/*NB03*corrected*" \
        "$BASE_DIR/*eke*corrected*" \
        "$BASE_DIR/*EKE*corrected*" \
        "$BASE_DIR/*corrected*eke*" \
        "$BASE_DIR/*NB03*"; do
        local matches
        matches=$(ls $pattern 2>/dev/null | head -3)
        if [ -n "$matches" ]; then
            log "  ${GREEN}[FOUND]${NC} EKE corrected output:"
            echo "$matches" | while read -r f; do
                log "          $f ($(du -h "$f" | cut -f1))"
            done
            ((pass++))
            break
        fi
    done
    # Also check for the log file that confirms the corrected pipeline
    if [ -f "$BASE_DIR/NB03_corrected.log" ]; then
        log "  ${GREEN}[FOUND]${NC} NB03_corrected.log (pipeline validation log)"
    else
        local logmatch
        logmatch=$(find "$BASE_DIR" -name "*NB03*corrected*log*" -o -name "*corrected*.log" 2>/dev/null | head -1)
        if [ -n "$logmatch" ]; then
            log "  ${GREEN}[FOUND]${NC} Corrected pipeline log: $logmatch"
        else
            log "  ${YELLOW}[WARN]${NC} NB03_corrected.log not found at expected location"
            ((warn++))
        fi
    fi

    # --- B.4: Half-peak-height width time series ---
    log ""
    log "--- B.4: Half-Peak-Height Width Time Series (NB02) ---"
    for pattern in \
        "$BASE_DIR/*NB02*patch_width*" \
        "$BASE_DIR/*patch_width*" \
        "$BASE_DIR/*half_peak*width*" \
        "$BASE_DIR/*width*" \
        "$BASE_DIR/*NB02*"; do
        local matches
        matches=$(ls $pattern 2>/dev/null | head -5)
        if [ -n "$matches" ]; then
            log "  ${GREEN}[FOUND]${NC} Width outputs:"
            echo "$matches" | while read -r f; do
                log "          $f ($(du -h "$f" | cut -f1))"
            done
            ((pass++))
            break
        fi
    done

    # Check for patch_width scripts
    for script in \
        "$BASE_DIR/patch_width_all_thresholds.py" \
        "$BASE_DIR/NB02_patch_width.py" \
        "$BASE_DIR/*patch_width*.py"; do
        local smatch
        smatch=$(ls $script 2>/dev/null | head -1)
        if [ -n "$smatch" ]; then
            log "  ${GREEN}[FOUND]${NC} Width script: $smatch"
            break
        fi
    done

    # --- B.5: Section boundary masks ---
    log ""
    log "--- B.5: Section Boundary / Meander Envelope Masks ---"
    for pattern in \
        "$BASE_DIR/*mask*" \
        "$BASE_DIR/*section*" \
        "$BASE_DIR/*boundary*" \
        "$BASE_DIR/*envelope*" \
        "$BASE_DIR/*upstream*" \
        "$BASE_DIR/*downstream*"; do
        local matches
        matches=$(ls $pattern 2>/dev/null | head -5)
        if [ -n "$matches" ]; then
            log "  ${GREEN}[FOUND]${NC} Mask/section files:"
            echo "$matches" | while read -r f; do
                log "          $f ($(du -h "$f" | cut -f1))"
            done
            ((pass++))
            break
        fi
    done

    # --- B.6: ERA5 wind speed trend fields ---
    log ""
    log "--- B.6: ERA5 Wind Speed Trend Fields ---"
    for searchdir in "$BASE_DIR" "/g/data/gv90/xl1657"; do
        for pattern in \
            "$searchdir/*ERA5*wind*trend*" \
            "$searchdir/*era5*wind*trend*" \
            "$searchdir/*wind*speed*trend*" \
            "$searchdir/*ERA5*" \
            "$searchdir/*era5*"; do
            local matches
            matches=$(ls $pattern 2>/dev/null | head -3)
            if [ -n "$matches" ]; then
                log "  ${GREEN}[FOUND]${NC} ERA5/wind files:"
                echo "$matches" | while read -r f; do
                    log "          $f ($(du -h "$f" | cut -f1))"
                done
                ((pass++))
                break 2
            fi
        done
    done

    # --- B.7: Threshold sensitivity outputs (SI Table S1 / Figure S1) ---
    log ""
    log "--- B.7: Threshold Sensitivity Outputs ---"
    for pattern in \
        "$BASE_DIR/*threshold*" \
        "$BASE_DIR/*sensitivity*" \
        "$BASE_DIR/*15*40*" \
        "$BASE_DIR/*SI*"; do
        local matches
        matches=$(ls $pattern 2>/dev/null | head -5)
        if [ -n "$matches" ]; then
            log "  ${GREEN}[FOUND]${NC} Threshold/sensitivity files:"
            echo "$matches" | while read -r f; do
                log "          $f ($(du -h "$f" | cut -f1))"
            done
            ((pass++))
            break
        fi
    done

    # --- B.8: GitHub repository sync check ---
    log ""
    log "--- B.8: GitHub Repository Status ---"
    local repo_dir
    repo_dir=$(find "/g/data/gv90/xl1657" -maxdepth 3 -type d -name "acc-meander-circumpolar-trends" 2>/dev/null | head -1)
    if [ -n "$repo_dir" ]; then
        log "  ${GREEN}[FOUND]${NC} Repository: $repo_dir"
        if [ -d "$repo_dir/.git" ]; then
            local branch
            branch=$(cd "$repo_dir" && git branch --show-current 2>/dev/null || echo "unknown")
            local last_commit
            last_commit=$(cd "$repo_dir" && git log -1 --format="%H %s" 2>/dev/null || echo "unknown")
            local status
            status=$(cd "$repo_dir" && git status --porcelain 2>/dev/null | wc -l)
            log "  Branch: $branch"
            log "  Last commit: $last_commit"
            log "  Uncommitted changes: $status files"
        fi
    else
        log "  ${YELLOW}[WARN]${NC} Repository acc-meander-circumpolar-trends not found on Gadi"
        ((warn++))
    fi

    # --- Summary ---
    log ""
    section_header "VERIFICATION SUMMARY"
    log "  ${GREEN}Passed:${NC}  $pass checks"
    log "  ${RED}Failed:${NC}  $fail checks"
    log "  ${YELLOW}Warnings:${NC} $warn checks"
    log ""
    log "  Full log saved to: $LOG_FILE"
    log ""
    if [ "$fail" -gt 0 ]; then
        log "  ${RED}ACTION REQUIRED: Review failed checks above.${NC}"
        log "  After running Stage A, update the placeholder paths in this script"
        log "  and re-run to confirm all datasets are accessible."
    else
        log "  ${GREEN}All critical checks passed. Ready to proceed to Step 0.2.${NC}"
    fi
}


# ==============================================================================
# STAGE C: NetCDF Integrity Spot-Check (runs inside Python venv)
# Purpose: Open a sample of NetCDF files and verify they load without errors,
# have expected dimensions, and cover the expected time range.
# ==============================================================================
stage_c() {
    section_header "STAGE C: NetCDF Integrity Spot-Check (Python)"
    log "Activating virtual environment: $VENV_DIR"

    if [ ! -f "$VENV_DIR/bin/python" ]; then
        log "${RED}Python venv not found. Skipping Stage C.${NC}"
        return
    fi

    "$VENV_DIR/bin/python" << 'PYEOF'
import os
import sys
import glob

base = "/g/data/gv90/xl1657/cmems_adt"

# Find all NetCDF files
nc_files = sorted(glob.glob(os.path.join(base, "**", "*.nc"), recursive=True))
nc_files += sorted(glob.glob(os.path.join(base, "**", "*.nc4"), recursive=True))

print(f"\n  Total NetCDF files found: {len(nc_files)}")

if len(nc_files) == 0:
    print("  No NetCDF files to inspect.")
    sys.exit(0)

try:
    import xarray as xr
    import numpy as np
except ImportError as e:
    print(f"  ERROR: Cannot import required package: {e}")
    print("  Install with: pip install xarray netCDF4 numpy")
    sys.exit(1)

# Inspect up to 10 representative files
sample = nc_files[:10] if len(nc_files) > 10 else nc_files

for f in sample:
    print(f"\n  --- {os.path.basename(f)} ---")
    print(f"      Path: {f}")
    print(f"      Size: {os.path.getsize(f) / 1e6:.1f} MB")
    try:
        ds = xr.open_dataset(f)
        print(f"      Variables: {list(ds.data_vars)}")
        print(f"      Dimensions: {dict(ds.dims)}")
        coords = list(ds.coords)
        print(f"      Coordinates: {coords}")

        # Check for time dimension
        time_names = [c for c in coords if c.lower() in ('time', 't', 'date')]
        if time_names:
            t = ds[time_names[0]]
            print(f"      Time range: {str(t.values[0])[:10]} to {str(t.values[-1])[:10]}")
            print(f"      Time steps: {len(t)}")

        # Check for lat/lon
        lat_names = [c for c in coords if c.lower() in ('latitude', 'lat', 'y')]
        lon_names = [c for c in coords if c.lower() in ('longitude', 'lon', 'x')]
        if lat_names:
            lat = ds[lat_names[0]]
            print(f"      Lat range: {float(lat.min()):.2f} to {float(lat.max()):.2f}")
        if lon_names:
            lon = ds[lon_names[0]]
            print(f"      Lon range: {float(lon.min()):.2f} to {float(lon.max()):.2f}")

        # Check for NaN fraction in first variable
        first_var = list(ds.data_vars)[0]
        sample_data = ds[first_var].isel(
            **{d: 0 for d in ds[first_var].dims if d.lower() not in ('latitude', 'lat', 'longitude', 'lon', 'y', 'x')}
        )
        nan_frac = float(np.isnan(sample_data.values).mean())
        print(f"      NaN fraction (first time step, '{first_var}'): {nan_frac:.3f}")

        ds.close()
        print(f"      Status: OK")

    except Exception as e:
        print(f"      ERROR: {e}")

print("\n  Stage C complete.\n")
PYEOF
}


# ==============================================================================
# Main execution
# ==============================================================================
main() {
    log "======================================================================"
    log "Step 0.1 — Verify Existing Derived Datasets on NCI Gadi"
    log "Date: $(date)"
    log "User: $(whoami)"
    log "Host: $(hostname)"
    log "======================================================================"

    # Always run Stage A
    stage_a

    # Run Stage B and C unless --scan flag is passed
    if [ "${1:-}" != "--scan" ]; then
        stage_b
        stage_c
    else
        log ""
        log "${YELLOW}Scan-only mode: Stages B and C skipped.${NC}"
        log "Run without --scan to perform full verification."
    fi

    log ""
    log "======================================================================"
    log "Verification complete. Full log: $LOG_FILE"
    log "======================================================================"
}

main "$@"
