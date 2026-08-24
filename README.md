# Air Pollution and Marathon Performance

Code repository for the research study by Donaldson et al. examining the association between ambient air pollution (NO2 and PM2.5) and marathon finish times across the six World Marathon Majors (2010–2024).

This repository is provided in accordance with the BMJ's data sharing and code availability policy.

## Study Overview

- **Design:** Retrospective observational study
- **Exposure:** Race-day NO2 and PM2.5 concentrations (Copernicus CAMS Global Reanalysis EAC4)
- **Outcome:** Individual marathon finish time
- **Marathons:** Boston, London, Berlin, Chicago, New York City, Tokyo
- **Period:** 2010–2024
- **Statistical approach:** Linear mixed models with random intercepts per marathon-year event, adjusted for heat index (natural spline), wind speed, and secular trends; linear quantile mixed models at multiple performance percentiles; subgroup analyses by sex, age, and performance level



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


### R

- `data.table`, `lme4`, `lmerTest`, `lqmm`, `splines`, `ggplot2`, `dplyr`, `parallel`



## Data Availability

Raw runner results were obtained from publicly available marathon results databases. Pollution data are from the Copernicus Atmosphere Monitoring Service (CAMS) Global Reanalysis (EAC4). Weather data are from the Open-Meteo historical weather API.

## License

This code is made available for the purposes of transparency and reproducibility of the associated research study.

