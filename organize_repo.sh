#!/bin/bash
# =============================================================================
# organize_repo.sh
# Reorganizes the Zenodo archive files into a clean GitHub repository structure
#
# Usage:
#   1. Unzip 19398337.zip into a temporary directory
#   2. Run this script from the repository root
#   3. It will move files from the temp directory into the proper structure
#
# Run: bash organize_repo.sh /path/to/unzipped/files
# =============================================================================

set -e

SRC="${1:-.}"
REPO="$(pwd)"

echo "======================================================================"
echo "Organizing repository: ${REPO}"
echo "Source files: ${SRC}"
echo "======================================================================"

# Create directory structure
echo "Creating directory structure..."
mkdir -p scripts
mkdir -p pbs
mkdir -p data/masks
mkdir -p data/eke_trends
mkdir -p data/sst_trends
mkdir -p data/statistics
mkdir -p manuscript
mkdir -p logs

# ── Scripts ──
echo "Moving analysis scripts..."
for f in setup_meander_sst_project.py download_cmems_sst.py \
         NB03b_spatial_eke.py NB12_regrid_sst.py NB13_sst_trends.py \
         NB14_fix_1995_and_par_transect.py NB15_manuscript_figures.py \
         NB16_fix_1995_complete.py NB17_supporting_information.py \
         verify_step01.sh; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" scripts/ && echo "  scripts/${f}"
done

# Remove duplicate NB15 (keep only the final version as NB15_manuscript_figures.py)
if [ -f "${SRC}/NB15_manuscript_figures_final.py" ]; then
  cp "${SRC}/NB15_manuscript_figures_final.py" scripts/NB15_manuscript_figures.py
  echo "  scripts/NB15_manuscript_figures.py (from final version)"
fi

# ── PBS scripts ──
echo "Moving PBS scripts..."
for f in download_sst.pbs run_spatial_eke.pbs run_regrid.pbs \
         run_sst_trends.pbs run_fix.pbs; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" pbs/ && echo "  pbs/${f}"
done

# ── Data: Masks ──
echo "Moving mask files..."
for f in meander_envelope_mask_CP.nc meander_envelope_mask_PAR.nc \
         meander_envelope_mask_SEIR.nc meander_envelope_mask_SWIR.nc \
         control_mask_CTRL_SE_PAC.nc control_mask_CTRL_S_IND.nc \
         control_mask_CTRL_S_ATL.nc control_mask_CTRL_CP_SEIR.nc \
         circumpolar_meander_envelope.nc; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" data/masks/ && echo "  data/masks/${f}"
done

# ── Data: EKE trends ──
echo "Moving EKE trend files..."
for f in spatial_eke_CP_trend.nc spatial_eke_PAR_trend.nc \
         spatial_eke_SEIR_trend.nc spatial_eke_SWIR_trend.nc; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" data/eke_trends/ && echo "  data/eke_trends/${f}"
done

# ── Data: SST trends ──
echo "Moving SST trend files..."
for f in sst_trend_map.nc sst_trend_along_acc_transect.nc; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" data/sst_trends/ && echo "  data/sst_trends/${f}"
done

# ── Data: Statistics (CSV) ──
echo "Moving statistics files..."
for f in sst_trend_stats.csv ks_test_meander_vs_control.csv \
         sst_eke_correlation.csv along_acc_site_anomalies.csv \
         decadal_decomposition.csv table_s1_control_regions.csv \
         table_s2_detailed_stats.csv; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" data/statistics/ && echo "  data/statistics/${f}"
done

# ── Manuscript ──
echo "Moving manuscript files..."
[ -f "${SRC}/manuscript.tex" ] && cp "${SRC}/manuscript.tex" manuscript/ && echo "  manuscript/manuscript.tex"
[ -f "${SRC}/table_s1_control_regions.tex" ] && cp "${SRC}/table_s1_control_regions.tex" manuscript/ && echo "  manuscript/table_s1_control_regions.tex"

# ── Logs ──
echo "Moving log files..."
for f in spatial_eke.out spatial_eke.err regrid_sst.out regrid_sst.err \
         sst_trends.out sst_trends.err fix_1995_par.out fix_1995_par.err; do
  [ -f "${SRC}/${f}" ] && cp "${SRC}/${f}" logs/ && echo "  logs/${f}"
done

# ── Config ──
echo "Moving config files..."
[ -f "${SRC}/project_config.json" ] && cp "${SRC}/project_config.json" ./ && echo "  project_config.json"

# ── Note about large monthly EKE files ──
echo ""
echo "======================================================================"
echo "NOTE: The following large files are NOT included in the repository"
echo "because they are intermediate products regenerable from the scripts:"
echo "  spatial_eke_CP_monthly.nc  (53 MB)"
echo "  spatial_eke_PAR_monthly.nc (67 MB)"
echo "  spatial_eke_SEIR_monthly.nc (21 MB)"
echo "  spatial_eke_SWIR_monthly.nc (34 MB)"
echo ""
echo "To regenerate, run: python scripts/NB03b_spatial_eke.py"
echo "======================================================================"

# Summary
echo ""
echo "Repository structure:"
find . -not -path './.git/*' -not -name '.git' | sort | head -60
echo ""
echo "Total size:"
du -sh .
echo ""
echo "Done. Ready for git init and first commit."
