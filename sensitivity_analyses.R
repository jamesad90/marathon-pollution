#!/usr/bin/env Rscript
################################################################################
# SENSITIVITY ANALYSES
#
# Mentioned in methods but results not reported:
# 1. Alternative heat index specifications (quadratic, splines 3-7 df)
# 2. Exclusion of extreme pollution events (>95th or <5th percentile)
# 3. Alternative exposure windows (race start vs median finish time)
################################################################################

library(data.table)
library(lme4)
library(lmerTest)
library(splines)

cat(rep("=", 80), "\n")
cat("SENSITIVITY ANALYSES\n")
cat(rep("=", 80), "\n\n")

################################################################################
# 1. LOAD AND PREPARE DATA
################################################################################

cat("Loading data...\n")
df <- fread("merged_marathon_dataset_with_tokyo.csv")
df <- df[time_minutes >= 120 & time_minutes <= 480]
df <- df[!(marathon == "london" & year == 2020)]
df[, event := paste(marathon, year, sep = "_")]
df[, time_seconds := time_minutes * 60]
df[, sex := as.factor(sex)]
df[, year_factor := as.factor(year)]

calculate_heat_index <- function(temp_c, rh) {
  T <- temp_c * 9/5 + 32
  RH <- rh
  HI <- 0.5 * (T + 61.0 + ((T - 68.0) * 1.2) + (RH * 0.094))
  high_temp <- (T >= 80)
  if (any(high_temp, na.rm = TRUE)) {
    HI_full <- -42.379 + 2.04901523 * T + 10.14333127 * RH - 0.22475541 * T * RH -
               6.83783e-3 * T^2 - 5.481717e-2 * RH^2 + 1.22874e-3 * T^2 * RH +
               8.5282e-4 * T * RH^2 - 1.99e-6 * T^2 * RH^2
    adj1 <- ifelse((RH < 13) & (T >= 80) & (T <= 112),
                   -((13 - RH) / 4) * sqrt((17 - abs(T - 95)) / 17), 0)
    adj2 <- ifelse((RH > 85) & (T >= 80) & (T <= 87),
                   ((RH - 85) / 10) * ((87 - T) / 5), 0)
    HI_full <- HI_full + adj1 + adj2
    HI[high_temp] <- HI_full[high_temp]
  }
  HI_c <- (HI - 32) * 5/9
  return(HI_c)
}

df[, heat_index := calculate_heat_index(temp_race_day, humidity_race_day)]

cat(sprintf("Total sample: %s runners\n\n", format(nrow(df), big.mark = ",")))

################################################################################
# HELPER: fit model and extract pollutant coefficient
################################################################################

extract_coef <- function(model, pollutant_col) {
  coef_val <- fixef(model)[pollutant_col]
  ci <- tryCatch(
    confint(model, parm = pollutant_col, method = "profile", quiet = TRUE),
    error = function(e) confint(model, parm = pollutant_col, method = "Wald", quiet = TRUE)
  )
  return(data.table(
    coefficient = coef_val,
    ci_lower = ci[1],
    ci_upper = ci[2]
  ))
}

################################################################################
# 2. MAIN MODEL (baseline for comparison)
################################################################################

cat(rep("=", 60), "\n")
cat("BASELINE MODEL (main analysis specification)\n")
cat(rep("=", 60), "\n\n")

all_results <- list()

for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")

  for (sex_val in c("M", "F")) {
    sex_label <- ifelse(sex_val == "M", "Male", "Female")

    data_sex <- df[sex == sex_val & !is.na(get(poll)) &
                   !is.na(heat_index) & !is.na(wind_speed_race_day)]

    mod <- lmer(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
      poll)), data = data_sex, REML = TRUE)

    res <- extract_coef(mod, poll)
    res[, `:=`(analysis = "Main model (ns df=5)", pollutant = poll_label,
               sex = sex_val, n = nrow(data_sex))]
    all_results[[paste("main", poll_label, sex_val)]] <- res

    cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
  }
}

################################################################################
# 3. SENSITIVITY 1: Alternative heat index specifications
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("SENSITIVITY 1: Alternative heat index specifications\n")
cat(rep("=", 60), "\n\n", sep = "")

# 3a. Quadratic heat index
cat("--- Quadratic heat index ---\n")
for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")
  for (sex_val in c("M", "F")) {
    sex_label <- ifelse(sex_val == "M", "Male", "Female")
    data_sex <- df[sex == sex_val & !is.na(get(poll)) &
                   !is.na(heat_index) & !is.na(wind_speed_race_day)]

    mod <- lmer(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + heat_index + I(heat_index^2) + wind_speed_race_day + (1|event)",
      poll)), data = data_sex, REML = TRUE)

    res <- extract_coef(mod, poll)
    res[, `:=`(analysis = "Quadratic heat index", pollutant = poll_label,
               sex = sex_val, n = nrow(data_sex))]
    all_results[[paste("quad", poll_label, sex_val)]] <- res

    cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
  }
}

# 3b. Natural splines with varying df
for (spline_df in c(3, 4, 6, 7)) {
  cat(sprintf("\n--- Natural spline df=%d ---\n", spline_df))
  for (poll in c("no2_race_day", "pm25_race_day")) {
    poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")
    for (sex_val in c("M", "F")) {
      sex_label <- ifelse(sex_val == "M", "Male", "Female")
      data_sex <- df[sex == sex_val & !is.na(get(poll)) &
                     !is.na(heat_index) & !is.na(wind_speed_race_day)]

      mod <- lmer(as.formula(sprintf(
        "time_seconds ~ %s + year_factor + ns(heat_index, df=%d) + wind_speed_race_day + (1|event)",
        poll, spline_df)), data = data_sex, REML = TRUE)

      res <- extract_coef(mod, poll)
      res[, `:=`(analysis = sprintf("Spline df=%d", spline_df), pollutant = poll_label,
                 sex = sex_val, n = nrow(data_sex))]
      all_results[[paste("ns", spline_df, poll_label, sex_val)]] <- res

      cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                  poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
    }
  }
}

################################################################################
# 4. SENSITIVITY 2: Exclusion of extreme pollution
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("SENSITIVITY 2: Excluding extreme pollution events\n")
cat(rep("=", 60), "\n\n", sep = "")

for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")

  p5 <- quantile(df[[poll]], 0.05, na.rm = TRUE)
  p95 <- quantile(df[[poll]], 0.95, na.rm = TRUE)
  cat(sprintf("%s: excluding < %.2f (5th) or > %.2f (95th)\n", poll_label, p5, p95))

  df_trimmed <- df[get(poll) >= p5 & get(poll) <= p95]
  cat(sprintf("  Remaining: %s runners (removed %s)\n",
              format(nrow(df_trimmed), big.mark = ","),
              format(nrow(df) - nrow(df_trimmed), big.mark = ",")))

  for (sex_val in c("M", "F")) {
    sex_label <- ifelse(sex_val == "M", "Male", "Female")
    data_sex <- df_trimmed[sex == sex_val & !is.na(get(poll)) &
                           !is.na(heat_index) & !is.na(wind_speed_race_day)]

    mod <- lmer(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
      poll)), data = data_sex, REML = TRUE)

    res <- extract_coef(mod, poll)
    res[, `:=`(analysis = "Excl. extreme pollution (5-95th)", pollutant = poll_label,
               sex = sex_val, n = nrow(data_sex))]
    all_results[[paste("trim", poll_label, sex_val)]] <- res

    cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
  }
}

################################################################################
# 5. SENSITIVITY 3: Alternative exposure windows
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("SENSITIVITY 3: Alternative exposure windows\n")
cat(rep("=", 60), "\n\n", sep = "")

# Check what pollution columns are available
poll_cols <- names(df)[grepl("no2|pm25|pm2", names(df), ignore.case = TRUE)]
cat("Available pollution columns:\n")
for (col in poll_cols) {
  n_valid <- sum(!is.na(df[[col]]))
  if (n_valid > 0) {
    cat(sprintf("  %s: %s valid values, mean=%.2f\n", col,
                format(n_valid, big.mark = ","), mean(df[[col]], na.rm = TRUE)))
  }
}

# Check what weather columns are available for alternative windows
weather_cols <- names(df)[grepl("temp|humid|wind|weather", names(df), ignore.case = TRUE)]
cat("\nAvailable weather columns:\n")
for (col in weather_cols) {
  n_valid <- sum(!is.na(df[[col]]))
  if (n_valid > 0) {
    cat(sprintf("  %s: %s valid values\n", col, format(n_valid, big.mark = ",")))
  }
}

# If alternative exposure columns exist, run models with them
alt_no2_cols <- poll_cols[grepl("no2", poll_cols) & poll_cols != "no2_race_day"]
alt_pm25_cols <- poll_cols[grepl("pm25|pm2\\.5", poll_cols) & poll_cols != "pm25_race_day"]

if (length(alt_no2_cols) > 0 || length(alt_pm25_cols) > 0) {
  cat("\nRunning models with alternative exposure windows...\n")

  for (alt_col in c(alt_no2_cols, alt_pm25_cols)) {
    poll_label <- ifelse(grepl("no2", alt_col), "NO2", "PM2.5")
    cat(sprintf("\n--- %s (column: %s) ---\n", poll_label, alt_col))

    for (sex_val in c("M", "F")) {
      sex_label <- ifelse(sex_val == "M", "Male", "Female")
      data_sex <- df[sex == sex_val & !is.na(get(alt_col)) &
                     !is.na(heat_index) & !is.na(wind_speed_race_day)]

      if (nrow(data_sex) < 1000) {
        cat(sprintf("  %s: insufficient data\n", sex_label))
        next
      }

      mod <- lmer(as.formula(sprintf(
        "time_seconds ~ %s + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
        alt_col)), data = data_sex, REML = TRUE)

      res <- extract_coef(mod, alt_col)
      res[, `:=`(analysis = sprintf("Alt exposure: %s", alt_col), pollutant = poll_label,
                 sex = sex_val, n = nrow(data_sex))]
      all_results[[paste("alt", alt_col, sex_val)]] <- res

      cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                  poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
    }
  }
} else {
  cat("\nNo alternative exposure columns found in dataset.\n")
  cat("Note: To run alternative exposure window analyses, the dataset would need\n")
  cat("pollution values for different time windows (e.g., race start hour, median finish hour).\n")
  cat("This could be added by re-querying CAMS for different hourly windows.\n")
}

################################################################################
# 6. SENSITIVITY 4: Linear heat index (no transformation)
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("SENSITIVITY 4: Linear heat index (no transformation)\n")
cat(rep("=", 60), "\n\n", sep = "")

for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")
  for (sex_val in c("M", "F")) {
    sex_label <- ifelse(sex_val == "M", "Male", "Female")
    data_sex <- df[sex == sex_val & !is.na(get(poll)) &
                   !is.na(heat_index) & !is.na(wind_speed_race_day)]

    mod <- lmer(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + heat_index + wind_speed_race_day + (1|event)",
      poll)), data = data_sex, REML = TRUE)

    res <- extract_coef(mod, poll)
    res[, `:=`(analysis = "Linear heat index", pollutant = poll_label,
               sex = sex_val, n = nrow(data_sex))]
    all_results[[paste("linear", poll_label, sex_val)]] <- res

    cat(sprintf("  %s %s: %.2f [%.2f, %.2f]\n",
                poll_label, sex_label, res$coefficient, res$ci_lower, res$ci_upper))
  }
}

################################################################################
# 7. SAVE ALL RESULTS
################################################################################

results_all <- rbindlist(all_results, fill = TRUE)
results_all <- results_all[, .(analysis, pollutant, sex, coefficient, ci_lower, ci_upper, n)]
fwrite(results_all, "results/sensitivity_analyses.csv")

cat("\n\nSaved: results/sensitivity_analyses.csv\n")

# Summary comparison table
cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY: NO2 coefficient across all sensitivity analyses\n")
cat(rep("=", 60), "\n\n", sep = "")

no2_summary <- results_all[pollutant == "NO2"]
for (s in c("M", "F")) {
  sex_label <- ifelse(s == "M", "Male", "Female")
  cat(sprintf("\n%s:\n", sex_label))
  d <- no2_summary[sex == s]
  for (i in 1:nrow(d)) {
    cat(sprintf("  %-40s %.2f [%.2f, %.2f]\n",
                d$analysis[i], d$coefficient[i], d$ci_lower[i], d$ci_upper[i]))
  }
}

cat("\nDONE.\n")
