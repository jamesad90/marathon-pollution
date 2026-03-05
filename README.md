# Air Pollution and Marathon Performance

Code repository for the research study examining the association between ambient air pollution (NO2 and PM2.5) and marathon finish times across the six World Marathon Majors (2010–2024).

This repository is provided in accordance with the BMJ's data sharing and code availability policy.

## Study Overview

- **Design:** Retrospective observational study
- **Exposure:** Race-day NO2 and PM2.5 concentrations (Copernicus CAMS Global Reanalysis EAC4)
- **Outcome:** Individual marathon finish time
- **Marathons:** Boston, London, Berlin, Chicago, New York City, Tokyo
- **Period:** 2010–2024
- **Statistical approach:** Linear mixed models with random intercepts per marathon-year event, adjusted for heat index (natural spline), wind speed, and secular trends; linear quantile mixed models at multiple performance percentiles; subgroup analyses by sex, age, and performance level

## Repository Structure

```
marathon-pollution/
├── api_scraper.py                     # Runner results collection (runzy.com API)
├── marathon_pollution_data.py         # Pollution data download (Copernicus CAMS/CDS API)
├── get_marathon_weather_openmeteo.py  # Weather data collection (Open-Meteo API)
├── merge_comprehensive_data.py        # Merges runner, pollution, and weather data
├── mixed_models_analysis.R            # Primary LMM and LQMM analyses
├── age_ability_mixed_models.R         # Subgroup analyses by age and performance level
├── sensitivity_analyses.R             # Sensitivity and robustness checks
├── two_pollutant_models.R             # Two-pollutant models and pollutant correlations
├── random_intercept_tests.R           # Likelihood ratio tests and marathon-specific models
├── pollution_nc_files/                # Raw CAMS netCDF files
└── marathon gpx files/                # GPX route files for each marathon
```

## Data Pipeline

Scripts should be run in the following order:

1. **`api_scraper.py`** — Scrapes individual runner finish times (name, time, sex, age, place) for all marathons and years from the runzy.com/marathonguide.com API.

2. **`marathon_pollution_data.py`** — Downloads race-day PM2.5, PM10, and NO2 concentrations from Copernicus CAMS EAC4 at four time points per day across a bounding box around each marathon city. Requires a Copernicus CDS API key (see below).

3. **`get_marathon_weather_openmeteo.py`** — Retrieves hourly temperature, humidity, wind speed, wind direction, precipitation, and surface pressure from the Open-Meteo historical archive, averaged over race hours (06:00–16:00).

4. **`merge_comprehensive_data.py`** — Joins runner results with pollution and weather data by marathon and year. Converts finish times to minutes and produces the final analysis-ready dataset.

## Statistical Analyses

All R scripts read from the merged dataset and write results to a `results/` directory.

| Script | Description |
|--------|-------------|
| `mixed_models_analysis.R` | Primary sex-stratified LMMs and LQMMs (1st–90th percentile) for NO2 and PM2.5 |
| `age_ability_mixed_models.R` | Stratified analyses by age group (18–29 to 60+) and performance level (elite to back-of-pack) |
| `sensitivity_analyses.R` | Alternative heat index specifications, exclusion of extreme pollution events, alternative exposure windows |
| `two_pollutant_models.R` | Simultaneous two-pollutant models; NO2–PM2.5 correlations |
| `random_intercept_tests.R` | Likelihood ratio tests for random effects; ICC; marathon-specific models |

## Requirements

### Python

- Python 3
- `pandas`, `numpy`, `requests`, `xarray`, `netCDF4`, `cdsapi`

### R

- `data.table`, `lme4`, `lmerTest`, `lqmm`, `splines`, `ggplot2`, `dplyr`, `parallel`

### Copernicus CDS API

To reproduce the pollution data download, a Copernicus Climate Data Store account and API key are required. Configure credentials in `~/.cdsapirc` as described at https://cds.climate.copernicus.eu/api-how-to.

## Data Availability

Raw runner results were obtained from publicly available marathon results databases. Pollution data are from the Copernicus Atmosphere Monitoring Service (CAMS) Global Reanalysis (EAC4). Weather data are from the Open-Meteo historical weather API.

## License

This code is made available for the purposes of transparency and reproducibility of the associated research study.
