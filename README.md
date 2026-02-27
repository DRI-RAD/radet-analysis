# RADET - Analysis

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![R 4.4+](https://img.shields.io/badge/R-4.4+-276DC3.svg?logo=r&logoColor=white)](https://cran.r-project.org/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green)](LICENSE)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)](https://earthengine.google.com/)
[![EarthArXiv Preprint](https://img.shields.io/badge/EarthArXiv-10.31223%2FX51B4P-blue)](https://doi.org/10.31223/X51B4P)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18225224.svg)](https://doi.org/10.5281/zenodo.18225224)

This repository contains analysis scripts for the **RADET model (Radiation Advection Diffusivity-independent Evapotranspiration)**, a Google Earth Engine-based approach for estimating actual evapotranspiration (ET) ([Kim et al., 2026](https://doi.org/10.31223/X51B4P)). For full model documentation, design details, input requirements, and supported collections, see the [radet-beta](https://github.com/DRI-RAD/radet-beta) repository.

### Analysis

The `analysis/` folder contains the following:

- [runtime_comparison.ipynb](analysis/runtime_comparison.ipynb) — Compare runtimes of OpenET models
- [eecu_analysis.py](analysis/eecu_analysis.py) — Analyze Earth Engine Compute Unit (EECU) usage across OpenET models

## Project Structure

```
radet-analysis/
├── analysis/
│   ├── README.md
│   ├── analysis_R/        # R scripts for data merging and figure generation
│   │   ├── README.md
│   │   ├── 1_merge_data.r
│   │   ├── 2_daily_figures.R
│   │   ├── 3_merge_data_monthly.r
│   │   ├── 4_monthly_figures.R
│   │   └── 5_KGE_improve_map.R
│   └── runtime_analysis/  # Runtime and EECU analysis
│       ├── README.md
│       ├── runtime_comparison.ipynb
│       ├── eecu_analysis.py
│       ├── eecu_data/     # Raw EECU input data
│       └── eecu_output/   # Generated analysis results and plots
├── .gitignore
├── LICENSE
└── README.md
```

## Python Dependencies

- [radet-beta](https://github.com/DRI-RAD/radet-beta) # main RADET model dependency
- [openet-sims](https://pypi.org/project/openet-sims/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-ssebop](https://pypi.org/project/openet-ssebop/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-ptjpl](https://pypi.org/project/openet-ptjpl/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-geesebal](https://pypi.org/project/openet-geesebal/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-disalexi](https://pypi.org/project/openet-disalexi/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [pandas](https://pypi.org/project/pandas/) # For analysis scripts
- [seaborn](https://seaborn.pydata.org/) # For analysis scripts

## Installation

### 1. Download and Install Anaconda/Miniconda

Either [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is required for managing Python packages (Python >= 3.10 recommended).

**Windows users:** After installation, open Anaconda Prompt and run `conda init powershell` to add conda to PowerShell.

**Linux/Mac users:** Ensure conda is added to your PATH (typically automatic). Restart your shell if needed.

Update conda: `conda update conda`

### 2. Create the Conda Environment

Create and activate a new conda environment:

```bash
conda create -y -n radet-analysis python=3.12
conda activate radet-analysis
pip install git+https://github.com/DRI-RAD/radet-beta.git
pip install pandas notebook seaborn openet-sims openet-ssebop openet-ptjpl openet-geesebal openet-disalexi
```

### Google Earth Engine Authentication

This project uses the Google Earth Engine (GEE) Python API for geospatial data extraction.

1. Install [Google Cloud CLI](https://cloud.google.com/sdk/docs/install-sdk)
2. Create a GCloud project (e.g., `gee-radet`) with GEE API enabled at https://console.cloud.google.com/
3. Configure the project:
   ```bash
   gcloud config set project gee-radet
   gcloud auth application-default set-quota-project gee-radet  # if prompted
   earthengine authenticate
   ```

See the [Earth Engine Python installation guide](https://developers.google.com/earth-engine/guides/python_install) for details.

## R Dependencies
Install the following packages before running the R scripts in [analysis_R](analysis/analysis_R/):

```r
install.packages(c(
  "tidyverse", "stringr", "lubridate", "Metrics", "hydroGOF",
  "ggpubr", "cowplot", "ggpmisc", "readr", "sf", "ggplot2",
  "rnaturalearth", "rnaturalearthdata", "ggspatial", "tigris", "scales"
))
```
## References

Kim, Y., Huntington, J. L., Comini de Andrade, B., Johnson, M. S., Volk, J. M., Majumdar, S., Morton, C., & ReVelle, P. (2026). Thermodynamically constrained surface energy balance using medium-resolution remote sensing for efficient evapotranspiration mapping. *EarthArXiv (preprint)*. [https://doi.org/10.31223/X51B4P](https://doi.org/10.31223/X51B4P)
