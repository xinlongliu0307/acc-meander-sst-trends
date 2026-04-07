# ACC Standing Meander SST Trends

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19398336.svg)](https://doi.org/10.5281/zenodo.19398336)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Python analysis scripts and derived data products for the manuscript:

> **Liu, X.** (2026). Topographic Hotspots Shape the Spatial Pattern of Southern Ocean SST Trends: Satellite Evidence from ACC Standing Meanders. *Geophysical Research Letters*.

This study co-locates 33 years (1993–2025) of satellite sea surface temperature (SST) trends with structural trends of Antarctic Circumpolar Current (ACC) standing meanders at four major topographic hotspots. We demonstrate that SST trends at meander sites differ significantly from quiescent ACC sections, with meander sites warming 0.047 °C per decade less than control regions (*D* = 0.181, *p* < 0.001).

## Study Sites

| Site | Abbreviation | Longitude Range | Dynamical Regime |
|------|-------------|-----------------|------------------|
| Campbell Plateau | CP | 150°E–150°W | Plateau-controlled |
| Pacific-Antarctic Ridge | PAR | 150°W–80°W | Ridge-controlled |
| Southeast Indian Ridge | SEIR | 130°E–152°E | Ridge-controlled |
| Southwest Indian Ridge | SWIR | 15°E–45°E | Ridge-controlled |

## Repository Structure

```
acc-meander-sst-trends/
├── README.md
├── LICENSE
├── .gitignore
├── environment.yml
├── project_config.json
│
├── scripts/                    # Analysis pipeline (numbered by execution order)
│   ├── setup_meander_sst_project.py   # Project directory setup and mask generation
│   ├── download_cmems_sst.py          # CMEMS OSTIA SST data download
│   ├── NB03b_spatial_eke.py           # Spatial EKE computation per site
│   ├── NB12_regrid_sst.py             # SST regridding (0.05° → 0.125°)
│   ├── NB13_sst_trends.py             # Per-grid-point SST trend computation
│   ├── NB14_fix_1995_and_par_transect.py  # 1995 fix + along-ACC transect
│   ├── NB15_manuscript_figures.py     # Main-text figure generation
│   ├── NB16_fix_1995_complete.py      # 1995 SST re-download and reprocessing
│   ├── NB17_supporting_information.py # Supporting Information outputs
│   └── verify_step01.sh               # Infrastructure verification
│
├── pbs/                        # NCI Gadi PBS job submission scripts
│   ├── download_sst.pbs
│   ├── run_spatial_eke.pbs
│   ├── run_regrid.pbs
│   ├── run_sst_trends.pbs
│   └── run_fix.pbs
│
├── data/                       # Derived data products (masks, trends, statistics)
│   ├── masks/
│   │   ├── meander_envelope_mask_CP.nc
│   │   ├── meander_envelope_mask_PAR.nc
│   │   ├── meander_envelope_mask_SEIR.nc
│   │   ├── meander_envelope_mask_SWIR.nc
│   │   ├── control_mask_CTRL_SE_PAC.nc
│   │   ├── control_mask_CTRL_S_IND.nc
│   │   ├── control_mask_CTRL_S_ATL.nc
│   │   ├── control_mask_CTRL_CP_SEIR.nc
│   │   └── circumpolar_meander_envelope.nc
│   │
│   ├── eke_trends/
│   │   ├── spatial_eke_CP_trend.nc
│   │   ├── spatial_eke_PAR_trend.nc
│   │   ├── spatial_eke_SEIR_trend.nc
│   │   └── spatial_eke_SWIR_trend.nc
│   │
│   ├── sst_trends/
│   │   ├── sst_trend_map.nc
│   │   └── sst_trend_along_acc_transect.nc
│   │
│   └── statistics/
│       ├── sst_trend_stats.csv
│       ├── ks_test_meander_vs_control.csv
│       ├── sst_eke_correlation.csv
│       ├── along_acc_site_anomalies.csv
│       ├── decadal_decomposition.csv
│       ├── table_s1_control_regions.csv
│       └── table_s2_detailed_stats.csv
│
├── manuscript/
│   └── table_s1_control_regions.tex
│
└── logs/                       # PBS job output logs (for reproducibility)
    ├── spatial_eke.out
    ├── spatial_eke.err
    ├── regrid_sst.out
    ├── regrid_sst.err
    ├── sst_trends.out
    ├── sst_trends.err
    ├── fix_1995_par.out
    └── fix_1995_par.err
```

## Analysis Pipeline

The analysis proceeds in the following order:

| Step | Script | Description | PBS Script |
|------|--------|-------------|------------|
| 0 | `verify_step01.sh` | Verify NCI Gadi infrastructure | — |
| 1 | `setup_meander_sst_project.py` | Create directory structure, meander envelope masks, control region masks | — |
| 2 | `download_cmems_sst.py` | Download CMEMS OSTIA SST (0.05°, daily, 1993–2025) | `download_sst.pbs` |
| 3 | `NB03b_spatial_eke.py` | Compute per-grid-point monthly EKE and Sen's slope trends | `run_spatial_eke.pbs` |
| 4 | `NB12_regrid_sst.py` | Regrid SST from 0.05° to 0.125° monthly | `run_regrid.pbs` |
| 5 | `NB13_sst_trends.py` | Compute per-grid-point SST trends (Sen's slope + Mann-Kendall) | `run_sst_trends.pbs` |
| 6 | `NB14_fix_1995_and_par_transect.py` | Recompute with all 33 years + corrected along-ACC transect | `run_fix.pbs` |
| 7 | `NB15_manuscript_figures.py` | Generate Figures 1–3 for the manuscript | — |
| 8 | `NB17_supporting_information.py` | Generate Supporting Information tables and Figure S1 | — |

## Data Sources

| Dataset | Resolution | Period | Source |
|---------|-----------|--------|--------|
| CMEMS DUACS L4 ADT | 0.125°, daily | 1993–2025 | [DOI: 10.48670/moi-00148](https://doi.org/10.48670/moi-00148) |
| CMEMS OSTIA SST | 0.05°, daily | 1993–2025 | [DOI: 10.48670/moi-00168](https://doi.org/10.48670/moi-00168) |
| GMRT Bathymetry | Variable | — | [GMRT](https://www.gmrt.org) |

## Requirements

- Python 3.11+
- Key packages: `xarray`, `numpy`, `scipy`, `matplotlib`, `cartopy`, `pymannkendall`, `statsmodels`, `pandas`
- NCI Gadi account with access to project `gv90` (for PBS job submission)
- `copernicusmarine` CLI tool (for SST data download)

## Installation

```bash
git clone https://github.com/xinlongliu0307/acc-meander-sst-trends.git
cd acc-meander-sst-trends
conda env create -f environment.yml
conda activate meander-sst
```

## Key Results

- **Meander–control SST trend difference:** −0.047 °C/decade (*D* = 0.181, *p* < 0.001)
- **PAR meander envelope:** −0.008 °C/decade (net cooling; 13.4% significant)
- **CP meander envelope:** +0.217 °C/decade (strongest warming; 84.4% significant)
- **EKE–SST sign reversal:** CP (*R* = +0.104) vs. PAR (*R* = −0.082)

## Citation

If you use these scripts or data products, please cite:

```bibtex
@article{liu2026meander_sst,
  title={Topographic Hotspots Shape the Spatial Pattern of Southern Ocean SST Trends: Satellite Evidence from ACC Standing Meanders},
  author={Liu, Xinlong},
  journal={Geophysical Research Letters},
  year={2026},
  publisher={Wiley Online Library}
}
```

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

## Acknowledgments

Computational resources were provided by the Australian National Computational Infrastructure (NCI). This study used E.U. Copernicus Marine Service Information. XL acknowledges support from the University of Tasmania and the Australian Antarctic Program Partnership.
