# RADET - Analysis

This folder contains R scripts for analyzing the performance of the [RADET model](https://github.com/DRI-RAD/radet-beta).

## Scripts

- [1_merge_data.r](1_merge_data.r) — This script integrates daily flux observations (https://zenodo.org/record/7636781), OpenET data (https://zenodo.org/record/7636781), and RADET data (https://doi.org/10.5281/zenodo.18225226). Before running the script, users must ensure that all required datasets are downloaded and properly stored.
- [2_daily_figures.R](2_daily_figures.R) — This script generates the daily-scale figures used in [Kim et al. (2026)](https://doi.org/10.31223/X51B4P). Script 1 must be executed before running this script.
- [3_merge_data_monthly.r](3_merge_data_monthly.r) — This script integrates monthly flux observations (https://zenodo.org/record/7636781), OpenET data (https://zenodo.org/record/7636781), and RADET data (https://doi.org/10.5281/zenodo.18225226). Before running the script, users must ensure that all required datasets are downloaded and properly stored.
- [4_monthly_figures.R](4_monthly_figures.R) — This script generates the monthly-scale figures used in [Kim et al. (2026)](https://doi.org/10.31223/X51B4P). Script 3 must be executed before running this script.
- [5_KGE_improve_map.R](5_KGE_improve_map.R) — This script generates the Map figures used in [Kim et al. (2026)](https://doi.org/10.31223/X51B4P). Scripts 2 and 4 must be executed before running this script.

## Required R Packages

Install the following packages before running the scripts:

```r
install.packages(c(
  "tidyverse", "stringr", "lubridate", "Metrics", "hydroGOF",
  "ggpubr", "cowplot", "ggpmisc", "readr", "sf", "ggplot2",
  "rnaturalearth", "rnaturalearthdata", "ggspatial", "tigris", "scales"
))
```

## References

Kim, Y., Huntington, J. L., Comini de Andrade, B., Johnson, M. S., Volk, J. M., Majumdar, S., Morton, C., & ReVelle, P. (2026). Thermodynamically constrained surface energy balance using medium-resolution remote sensing for efficient evapotranspiration mapping. *EarthArXiv (preprint)*. [https://doi.org/10.31223/X51B4P](https://doi.org/10.31223/X51B4P)
