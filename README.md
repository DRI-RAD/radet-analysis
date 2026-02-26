# RADET - Analysis

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-Apache%202.0-green)](LICENSE)
[![GEE](https://img.shields.io/badge/Google%20Earth%20Engine-4285F4?logo=google-earth&logoColor=white)](https://earthengine.google.com/)
[![EarthArXiv Preprint](https://img.shields.io/badge/EarthArXiv-10.31223%2FX51B4P-blue)](https://doi.org/10.31223/X51B4P)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18225226.svg)](https://doi.org/10.5281/zenodo.18225226)

This repository contains analysis scripts for the **RADET model (Radiation Advection Diffusivity-independent Evapotranspiration)**, a Google Earth Engine-based approach for estimating actual evapotranspiration (ET) ([Kim et al., 2026](https://doi.org/10.31223/X51B4P)). For full model documentation, design details, input requirements, and supported collections, see the [radet-beta](https://github.com/DRI-RAD/radet-beta) repository.

### Analysis

The `analysis/` folder contains the following:

- [runtime_comparison.ipynb](analysis/runtime_comparison.ipynb) — Compare runtimes of OpenET models
- [eecu_analysis.py](analysis/eecu_analysis.py) — Analyze Earth Engine Compute Unit (EECU) usage across OpenET models

## Note on RADET Code

The `radet/` directory in this repository is cloned from [https://github.com/DRI-RAD/radet-beta](https://github.com/DRI-RAD/radet-beta). This duplicated code will be removed once a standalone Python package is published.

## Project Structure

```
radet-analysis/
├── radet/                 # Cloned from https://github.com/DRI-RAD/radet-beta
│   ├── __init__.py        #   (will be removed once a Python package is published)
│   ├── collection.py
│   ├── image.py
│   ├── interpolate.py
│   ├── landsat.py
│   ├── model.py
│   └── utils.py
├── analysis/
│   ├── README.md
│   ├── runtime_comparison.ipynb
│   ├── eecu_analysis.py
│   ├── eecu_data/         # Raw EECU input data
│   └── eecu_output/       # Generated analysis results and plots
├── .gitignore
├── LICENSE
└── README.md
```

## Dependencies

- [earthengine-api](https://github.com/google/earthengine-api) # main RADET model dependency
- [openet-core](https://pypi.org/project/openet-core/) # main RADET model dependency
- [openet-sims](https://pypi.org/project/openet-sims/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-ssebop](https://pypi.org/project/openet-ssebop/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-ptjpl](https://pypi.org/project/openet-ptjpl/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-geesebal](https://pypi.org/project/openet-geesebal/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [openet-disalexi](https://pypi.org/project/openet-disalexi/) # Only for runtime comparisons ([runtime_comparison.ipynb](analysis/runtime_comparison.ipynb))
- [pandas](https://pypi.org/project/pandas/) # For analysis scripts
- [seaborn](https://seaborn.pydata.org/) # For analysis scripts

## Installation

```
pip install earthengine-api openet-core pandas seaborn openet-sims openet-ptjpl openet-ssebop openet-disalexi openet-geesebal
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

## References

Kim, Y., Huntington, J. L., Comini de Andrade, B., Johnson, M. S., Volk, J. M., Majumdar, S., Morton, C., & ReVelle, P. (2026). Thermodynamically constrained closed-form surface energy balance using medium-resolution remote sensing for efficient evapotranspiration mapping. *EarthArXiv (preprint)*. [https://doi.org/10.31223/X51B4P](https://doi.org/10.31223/X51B4P)
