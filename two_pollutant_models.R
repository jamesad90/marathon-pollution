#!/usr/bin/env Rscript
################################################################################
# TWO-POLLUTANT MODELS: NO2 + PM2.5 IN SAME MODEL
#
# Requested by reviewer:
# 1. Report NO2-PM2.5 correlation
# 2. Run two-pollutant LMMs (both pollutants in same model)
# 3. Compare effect sizes with single-pollutant models (attenuation?)
################################################################################

library(data.table)
library(lme4)
library(lmerTest)
library(splines)

cat(rep("=", 80), "\n")
cat("TWO-POLLUTANT MODEL ANALYSIS\n")
cat(rep("=", 80), "\n\n")

################################################################################
# 1. LOAD AND PREPARE DATA
################################################################################

cat("Loading data...\n")
df <- fread("merged_marathon_dataset_with_tokyo.csv")

# Filter valid times
df <- df[time_minutes >= 120 & time_minutes <= 480]
df <- df[!(marathon == "london" & year == 2020)]

# Create event identifier
df[, event := paste(marathon, year, sep = "_")]
df[, time_seconds := time_minutes * 60]
df[, sex := as.factor(sex)]
df[, year_factor := as.factor(year)]

# Heat index calculation (same as main analysis)
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

# Complete cases for both pollutants
df_both <- df[!is.na(no2_race_day) & !is.na(pm25_race_day) &
              !is.na(heat_index) & !is.na(wind_speed_race_day)]

cat(sprintf("Complete cases (both pollutants): %s runners\n",
            format(nrow(df_both), big.mark = ",")))

################################################################################
# 2. POLLUTANT CORRELATION
################################################################################

cat("\n", rep("-", 60), "\n", sep = "")
cat("POLLUTANT CORRELATION\n")
cat(rep("-", 60), "\n\n", sep = "")

# Overall correlation
cor_overall <- cor(df_both$no2_race_day, df_both$pm25_race_day, method = "pearson")
cor_spearman <- cor(df_both$no2_race_day, df_both$pm25_race_day, method = "spearman")
cor_test <- cor.test(df_both$no2_race_day, df_both$pm25_race_day)

cat(sprintf("Pearson correlation (NO2 vs PM2.5): r = %.3f (p %s)\n",
            cor_overall,
            ifelse(cor_test$p.value < 0.001, "< 0.001", sprintf("= %.3f", cor_test$p.value))))
cat(sprintf("Spearman correlation: rho = %.3f\n", cor_spearman))

# Correlation by marathon
cat("\nCorrelation by marathon (event-level):\n")
event_means <- df_both[, .(mean_no2 = mean(no2_race_day),
                            mean_pm25 = mean(pm25_race_day)), by = .(marathon, year)]
for (m in sort(unique(event_means$marathon))) {
  d <- event_means[marathon == m]
  if (nrow(d) >= 3) {
    r <- cor(d$mean_no2, d$mean_pm25)
    cat(sprintf("  %s: r = %.3f (n = %d events)\n", m, r, nrow(d)))
  }
}
cor_event <- cor(event_means$mean_no2, event_means$mean_pm25)
cat(sprintf("\n  Overall (event-level): r = %.3f (n = %d events)\n",
            cor_event, nrow(event_means)))

################################################################################
# 3. TWO-POLLUTANT MODELS
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("TWO-POLLUTANT LINEAR MIXED MODELS\n")
cat(rep("=", 60), "\n\n", sep = "")

results_list <- list()

for (sex_val in c("M", "F")) {
  sex_label <- ifelse(sex_val == "M", "Male", "Female")
  cat(sprintf("\n--- %s runners ---\n", sex_label))

  data_sex <- df_both[sex == sex_val]
  cat(sprintf("  N = %s\n", format(nrow(data_sex), big.mark = ",")))

  # Single-pollutant model: NO2 only
  cat("  Fitting single-pollutant NO2 model...\n")
  mod_no2_single <- lmer(time_seconds ~ no2_race_day + year_factor +
                         ns(heat_index, df = 5) + wind_speed_race_day + (1|event),
                         data = data_sex, REML = TRUE)
  coef_no2_single <- fixef(mod_no2_single)["no2_race_day"]
  ci_no2_single <- confint(mod_no2_single, parm = "no2_race_day", method = "profile", quiet = TRUE)

  # Single-pollutant model: PM2.5 only
  cat("  Fitting single-pollutant PM2.5 model...\n")
  mod_pm25_single <- lmer(time_seconds ~ pm25_race_day + year_factor +
                          ns(heat_index, df = 5) + wind_speed_race_day + (1|event),
                          data = data_sex, REML = TRUE)
  coef_pm25_single <- fixef(mod_pm25_single)["pm25_race_day"]
  ci_pm25_single <- confint(mod_pm25_single, parm = "pm25_race_day", method = "profile", quiet = TRUE)

  # Two-pollutant model: NO2 + PM2.5
  cat("  Fitting two-pollutant model...\n")
  mod_two <- lmer(time_seconds ~ no2_race_day + pm25_race_day + year_factor +
                  ns(heat_index, df = 5) + wind_speed_race_day + (1|event),
                  data = data_sex, REML = TRUE)
  coef_no2_two <- fixef(mod_two)["no2_race_day"]
  coef_pm25_two <- fixef(mod_two)["pm25_race_day"]
  ci_two <- confint(mod_two, parm = c("no2_race_day", "pm25_race_day"),
                    method = "profile", quiet = TRUE)

  # Calculate attenuation
  atten_no2 <- (1 - coef_no2_two / coef_no2_single) * 100
  atten_pm25 <- (1 - coef_pm25_two / coef_pm25_single) * 100

  cat(sprintf("\n  SINGLE-POLLUTANT RESULTS (on complete-case sample):\n"))
  cat(sprintf("    NO2:  %.2f sec/µg/m³ [%.2f, %.2f]\n",
              coef_no2_single, ci_no2_single[1], ci_no2_single[2]))
  cat(sprintf("    PM2.5: %.2f sec/µg/m³ [%.2f, %.2f]\n",
              coef_pm25_single, ci_pm25_single[1], ci_pm25_single[2]))

  cat(sprintf("\n  TWO-POLLUTANT RESULTS:\n"))
  cat(sprintf("    NO2:  %.2f sec/µg/m³ [%.2f, %.2f]\n",
              coef_no2_two, ci_two["no2_race_day", 1], ci_two["no2_race_day", 2]))
  cat(sprintf("    PM2.5: %.2f sec/µg/m³ [%.2f, %.2f]\n",
              coef_pm25_two, ci_two["pm25_race_day", 1], ci_two["pm25_race_day", 2]))

  cat(sprintf("\n  ATTENUATION:\n"))
  cat(sprintf("    NO2:  %.1f%% change (positive = attenuated)\n", atten_no2))
  cat(sprintf("    PM2.5: %.1f%% change (positive = attenuated)\n", atten_pm25))

  # Store results
  results_list[[paste0("single_no2_", sex_val)]] <- data.table(
    model = "single_pollutant", pollutant = "NO2", sex = sex_val,
    coefficient = coef_no2_single,
    ci_lower = ci_no2_single[1], ci_upper = ci_no2_single[2],
    n = nrow(data_sex)
  )
  results_list[[paste0("single_pm25_", sex_val)]] <- data.table(
    model = "single_pollutant", pollutant = "PM2.5", sex = sex_val,
    coefficient = coef_pm25_single,
    ci_lower = ci_pm25_single[1], ci_upper = ci_pm25_single[2],
    n = nrow(data_sex)
  )
  results_list[[paste0("two_no2_", sex_val)]] <- data.table(
    model = "two_pollutant", pollutant = "NO2", sex = sex_val,
    coefficient = coef_no2_two,
    ci_lower = ci_two["no2_race_day", 1], ci_upper = ci_two["no2_race_day", 2],
    n = nrow(data_sex)
  )
  results_list[[paste0("two_pm25_", sex_val)]] <- data.table(
    model = "two_pollutant", pollutant = "PM2.5", sex = sex_val,
    coefficient = coef_pm25_two,
    ci_lower = ci_two["pm25_race_day", 1], ci_upper = ci_two["pm25_race_day", 2],
    n = nrow(data_sex)
  )
}

################################################################################
# 4. SAVE RESULTS
################################################################################

results_all <- rbindlist(results_list)
fwrite(results_all, "results/two_pollutant_model_results.csv")

# Save correlation info
cor_data <- data.table(
  level = c("individual", "event"),
  pearson_r = c(cor_overall, cor_event),
  spearman_rho = c(cor_spearman, NA),
  n = c(nrow(df_both), nrow(event_means))
)
fwrite(cor_data, "results/pollutant_correlation.csv")

cat("\n\nFiles saved:\n")
cat("  results/two_pollutant_model_results.csv\n")
cat("  results/pollutant_correlation.csv\n")
cat("\nDONE.\n")
