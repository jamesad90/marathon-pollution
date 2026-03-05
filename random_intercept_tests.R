#!/usr/bin/env Rscript
################################################################################
# RANDOM INTERCEPT SIGNIFICANCE & MARATHON-SPECIFIC ANALYSES
#
# Requested by reviewer:
# 1. Test if random intercept term is significant (LRT)
# 2. If yes, run marathon-specific analyses
# 3. Compare absolute times of top performers across locations
################################################################################

library(data.table)
library(lme4)
library(lmerTest)
library(splines)

cat(rep("=", 80), "\n")
cat("RANDOM INTERCEPT TESTS & MARATHON-SPECIFIC ANALYSES\n")
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

################################################################################
# 2. RANDOM INTERCEPT SIGNIFICANCE (Likelihood Ratio Test)
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("RANDOM INTERCEPT SIGNIFICANCE TEST\n")
cat(rep("=", 60), "\n\n", sep = "")

lrt_results <- list()

for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")

  for (sex_val in c("M", "F")) {
    sex_label <- ifelse(sex_val == "M", "Male", "Female")
    cat(sprintf("\n--- %s, %s ---\n", poll_label, sex_label))

    data_sex <- df[sex == sex_val & !is.na(get(poll)) &
                   !is.na(heat_index) & !is.na(wind_speed_race_day)]

    # Model WITH random intercept
    mod_re <- lmer(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
      poll)), data = data_sex, REML = FALSE)  # ML for LRT

    # Model WITHOUT random intercept (fixed effects only)
    mod_fe <- lm(as.formula(sprintf(
      "time_seconds ~ %s + year_factor + ns(heat_index, df=5) + wind_speed_race_day",
      poll)), data = data_sex)

    # LRT comparing models
    loglik_re <- logLik(mod_re)
    loglik_fe <- logLik(mod_fe)
    lrt_stat <- as.numeric(-2 * (loglik_fe - loglik_re))
    # Chi-sq with 1 df (variance component), but boundary test so use mixture
    p_value <- 0.5 * pchisq(lrt_stat, df = 1, lower.tail = FALSE)

    # Random intercept variance
    vc <- as.data.frame(VarCorr(mod_re))
    re_var <- vc[vc$grp == "event", "vcov"]
    re_sd <- vc[vc$grp == "event", "sdcor"]
    resid_var <- vc[vc$grp == "Residual", "vcov"]
    icc <- re_var / (re_var + resid_var)

    cat(sprintf("  Random intercept SD: %.2f seconds (%.2f minutes)\n", re_sd, re_sd / 60))
    cat(sprintf("  ICC: %.4f (%.2f%% of variance explained by event)\n", icc, icc * 100))
    cat(sprintf("  LRT statistic: %.2f\n", lrt_stat))
    cat(sprintf("  p-value: %s\n", ifelse(p_value < 0.001, "< 0.001", sprintf("%.4f", p_value))))
    cat(sprintf("  Significant: %s\n", ifelse(p_value < 0.05, "YES", "NO")))

    lrt_results[[paste(poll_label, sex_val)]] <- data.table(
      pollutant = poll_label, sex = sex_val,
      re_sd_seconds = re_sd, re_sd_minutes = re_sd / 60,
      icc = icc, lrt_statistic = lrt_stat, p_value = p_value,
      significant = p_value < 0.05
    )
  }
}

lrt_df <- rbindlist(lrt_results)
fwrite(lrt_df, "results/random_intercept_lrt.csv")
cat("\nSaved: results/random_intercept_lrt.csv\n")

################################################################################
# 3. MARATHON-SPECIFIC ANALYSES
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("MARATHON-SPECIFIC ANALYSES\n")
cat(rep("=", 60), "\n\n", sep = "")

marathon_results <- list()

for (poll in c("no2_race_day", "pm25_race_day")) {
  poll_label <- ifelse(poll == "no2_race_day", "NO2", "PM2.5")
  cat(sprintf("\n=== %s ===\n", poll_label))

  for (m in sort(unique(df$marathon))) {
    data_m <- df[marathon == m & !is.na(get(poll)) &
                 !is.na(heat_index) & !is.na(wind_speed_race_day)]

    n_events <- length(unique(data_m$event))
    cat(sprintf("\n  %s: N=%s runners, %d events\n",
                m, format(nrow(data_m), big.mark = ","), n_events))

    if (n_events < 3) {
      cat("    Skipping - insufficient events\n")
      next
    }

    tryCatch({
      # Pooled model (both sexes)
      mod <- lmer(as.formula(sprintf(
        "time_seconds ~ %s + sex + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
        poll)), data = data_m, REML = TRUE)

      coef_val <- fixef(mod)[poll]
      ci <- confint(mod, parm = poll, method = "profile", quiet = TRUE)

      cat(sprintf("    Coefficient: %.2f sec/µg/m³ [%.2f, %.2f]\n",
                  coef_val, ci[1], ci[2]))

      marathon_results[[paste(poll_label, m)]] <- data.table(
        pollutant = poll_label, marathon = m,
        coefficient = coef_val,
        ci_lower = ci[1], ci_upper = ci[2],
        n_runners = nrow(data_m), n_events = n_events
      )
    }, error = function(e) {
      cat(sprintf("    ERROR: %s\n", e$message))
    })
  }
}

marathon_df <- rbindlist(marathon_results)
fwrite(marathon_df, "results/marathon_specific_lmm.csv")
cat("\nSaved: results/marathon_specific_lmm.csv\n")

################################################################################
# 4. ABSOLUTE TIMES BY MARATHON FOR TOP PERFORMERS
################################################################################

cat("\n", rep("=", 60), "\n", sep = "")
cat("ABSOLUTE TIMES BY MARATHON (TOP PERFORMERS)\n")
cat(rep("=", 60), "\n\n", sep = "")

# Performance percentiles within each event
df[, perf_pct := percent_rank(time_seconds), by = event]

# Top 10% (elite)
elite <- df[perf_pct <= 0.10]

elite_summary <- elite[, .(
  mean_time_min = mean(time_minutes),
  median_time_min = median(time_minutes),
  sd_time_min = sd(time_minutes),
  mean_no2 = mean(no2_race_day, na.rm = TRUE),
  mean_pm25 = mean(pm25_race_day, na.rm = TRUE),
  n_runners = .N,
  n_events = length(unique(event))
), by = marathon]

cat("Elite runners (top 10%) — mean finish time by marathon:\n")
print(elite_summary[order(mean_time_min)])

fwrite(elite_summary, "results/elite_times_by_marathon.csv")
cat("\nSaved: results/elite_times_by_marathon.csv\n")

# Also top 5% for a tighter elite group
top5 <- df[perf_pct <= 0.05]
top5_summary <- top5[, .(
  mean_time_min = mean(time_minutes),
  median_time_min = median(time_minutes),
  mean_no2 = mean(no2_race_day, na.rm = TRUE),
  mean_pm25 = mean(pm25_race_day, na.rm = TRUE),
  n_runners = .N
), by = marathon]

cat("\nTop 5% runners — mean finish time by marathon:\n")
print(top5_summary[order(mean_time_min)])

fwrite(top5_summary, "results/top5_times_by_marathon.csv")

cat("\nDONE.\n")
