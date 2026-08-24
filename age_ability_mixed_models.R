#!/usr/bin/env Rscript
################################################################################
# AGE AND ABILITY SUBGROUP ANALYSIS WITH LINEAR MIXED MODELS
#
# Explores heterogeneous pollution effects by:
# - Age groups (18-29, 30-39, 40-49, 50-59, 60+)
# - Performance levels (elite, sub-elite, mid-pack, recreational, back-of-pack)
# - Age × Performance interactions
#
# Uses same mixed model framework as main analysis:
# - Random intercepts for marathon-year events
# - Heat index control (natural splines, 5 df)
# - Wind speed control (linear)
# - Year fixed effects
################################################################################

# Load required libraries
library(data.table)
library(lme4)
library(splines)
library(ggplot2)
library(dplyr)

cat(rep("=", 80), "\n")
cat("AGE AND ABILITY SUBGROUP ANALYSIS\n")
cat(rep("=", 80), "\n\n")

################################################################################
# 1. LOAD AND PREPARE DATA
################################################################################

cat("Loading data...\n")
df <- fread("merged_marathon_dataset_with_tokyo.csv")

# Filter valid times and ages
df <- df[time_minutes >= 120 & time_minutes <= 480]
df <- df[!(marathon == "london" & year == 2020)]
df <- df[age >= 18 & age <= 100]  # Valid ages only

cat(sprintf("Valid data: %s runners\n", format(nrow(df), big.mark=",")))

################################################################################
# 2. CREATE AGE GROUPS AND PERFORMANCE LEVELS
################################################################################

cat("\nCreating age groups and performance levels...\n")

# Create event identifier first (needed for performance percentiles)
df[, event := paste(marathon, year, sep="_")]
df[, time_seconds := time_minutes * 60]

# Age groups (standard marathon categories)
df[, age_group := fcase(
  age >= 18 & age < 30, "18-29",
  age >= 30 & age < 40, "30-39",
  age >= 40 & age < 50, "40-49",
  age >= 50 & age < 60, "50-59",
  age >= 60, "60+",
  default = NA_character_
)]

df[, age_group := factor(age_group, levels = c("18-29", "30-39", "40-49", "50-59", "60+"))]

# Performance levels (within-race percentiles)
df[, performance_percentile := percent_rank(time_seconds), by = .(event)]

df[, performance_group := fcase(
  performance_percentile <= 0.10, "Elite (top 10%)",
  performance_percentile <= 0.25, "Sub-Elite (10-25%)",
  performance_percentile <= 0.50, "Upper Mid (25-50%)",
  performance_percentile <= 0.75, "Mid-Pack (50-75%)",
  performance_percentile <= 0.90, "Recreational (75-90%)",
  default = "Back of Pack (90-100%)"
)]

df[, performance_group := factor(performance_group,
                                  levels = c("Elite (top 10%)",
                                           "Sub-Elite (10-25%)",
                                           "Upper Mid (25-50%)",
                                           "Mid-Pack (50-75%)",
                                           "Recreational (75-90%)",
                                           "Back of Pack (90-100%)"))]

# Sample size check
cat("\nSample sizes by age group:\n")
print(df[, .N, by = age_group][order(age_group)])

cat("\nSample sizes by performance group:\n")
print(df[, .N, by = performance_group][order(performance_group)])

cat("\nSample sizes by age × performance (showing smallest cells):\n")
age_perf_counts <- df[, .N, by = .(age_group, performance_group)]
print(age_perf_counts[order(N)][1:10])

################################################################################
# 3. CALCULATE HEAT INDEX
################################################################################

cat("\nCalculating heat index...\n")

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
df[, sex := as.factor(sex)]
df[, year_factor := as.factor(year)]

################################################################################
# 4. AGE-STRATIFIED MIXED MODELS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("AGE-STRATIFIED LINEAR MIXED MODELS\n")
cat(rep("=", 80), "\n\n", sep="")

# Function to fit age-stratified models
fit_age_stratified <- function(data, pollutant_col, pollutant_name) {
  cat(sprintf("\n--- %s ---\n", pollutant_name))

  results_list <- list()

  for (age_grp in levels(data$age_group)) {
    cat(sprintf("\nAge group: %s\n", age_grp))

    data_age <- data[age_group == age_grp & !is.na(get(pollutant_col)) &
                     !is.na(heat_index) & !is.na(wind_speed_race_day) & !is.na(sex)]

    if (nrow(data_age) < 1000) {
      cat(sprintf("  Skipping - insufficient sample size (N=%d)\n", nrow(data_age)))
      next
    }

    cat(sprintf("  N = %s runners\n", format(nrow(data_age), big.mark=",")))

    # Pooled model (both sexes)
    formula_str <- sprintf("time_seconds ~ %s + sex + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
                          pollutant_col)

    tryCatch({
      model <- lmer(as.formula(formula_str), data = data_age, REML = TRUE)
      coef_val <- fixef(model)[pollutant_col]
      ci <- confint(model, parm = pollutant_col, method = "profile", quiet = TRUE)

      results_list[[age_grp]] <- data.table(
        pollutant = pollutant_name,
        age_group = age_grp,
        coefficient = coef_val,
        ci_lower = ci[1],
        ci_upper = ci[2],
        n = nrow(data_age)
      )

      cat(sprintf("  Coefficient: %.2f sec/µg/m³ [95%% CI: %.2f, %.2f]\n",
                  coef_val, ci[1], ci[2]))
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
    })
  }

  results_df <- rbindlist(results_list, fill = TRUE)
  return(results_df)
}

# Complete cases for each pollutant
df_no2 <- df[!is.na(no2_race_day) & !is.na(age_group)]
df_pm25 <- df[!is.na(pm25_race_day) & !is.na(age_group)]

# Fit age-stratified models
results_age_no2 <- fit_age_stratified(df_no2, "no2_race_day", "NO2")
results_age_pm25 <- fit_age_stratified(df_pm25, "pm25_race_day", "PM2.5")

# Combine and save
results_age_all <- rbind(results_age_no2, results_age_pm25)
fwrite(results_age_all, "results/age_stratified_lmm.csv")
cat("\nSaved: results/age_stratified_lmm.csv\n")

################################################################################
# 5. PERFORMANCE-STRATIFIED MIXED MODELS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("PERFORMANCE-STRATIFIED LINEAR MIXED MODELS\n")
cat(rep("=", 80), "\n\n", sep="")

# Function to fit performance-stratified models
fit_performance_stratified <- function(data, pollutant_col, pollutant_name) {
  cat(sprintf("\n--- %s ---\n", pollutant_name))

  results_list <- list()

  for (perf_grp in levels(data$performance_group)) {
    cat(sprintf("\nPerformance group: %s\n", perf_grp))

    data_perf <- data[performance_group == perf_grp & !is.na(get(pollutant_col)) &
                      !is.na(heat_index) & !is.na(wind_speed_race_day) & !is.na(sex)]

    if (nrow(data_perf) < 1000) {
      cat(sprintf("  Skipping - insufficient sample size (N=%d)\n", nrow(data_perf)))
      next
    }

    cat(sprintf("  N = %s runners\n", format(nrow(data_perf), big.mark=",")))

    # Pooled model (both sexes)
    formula_str <- sprintf("time_seconds ~ %s + sex + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
                          pollutant_col)

    tryCatch({
      model <- lmer(as.formula(formula_str), data = data_perf, REML = TRUE)
      coef_val <- fixef(model)[pollutant_col]
      ci <- confint(model, parm = pollutant_col, method = "profile", quiet = TRUE)

      results_list[[perf_grp]] <- data.table(
        pollutant = pollutant_name,
        performance_group = perf_grp,
        coefficient = coef_val,
        ci_lower = ci[1],
        ci_upper = ci[2],
        n = nrow(data_perf)
      )

      cat(sprintf("  Coefficient: %.2f sec/µg/m³ [95%% CI: %.2f, %.2f]\n",
                  coef_val, ci[1], ci[2]))
    }, error = function(e) {
      cat(sprintf("  ERROR: %s\n", e$message))
    })
  }

  results_df <- rbindlist(results_list, fill = TRUE)
  return(results_df)
}

# Complete cases
df_no2_perf <- df[!is.na(no2_race_day) & !is.na(performance_group)]
df_pm25_perf <- df[!is.na(pm25_race_day) & !is.na(performance_group)]

# Fit performance-stratified models
results_perf_no2 <- fit_performance_stratified(df_no2_perf, "no2_race_day", "NO2")
results_perf_pm25 <- fit_performance_stratified(df_pm25_perf, "pm25_race_day", "PM2.5")

# Combine and save
results_perf_all <- rbind(results_perf_no2, results_perf_pm25)
fwrite(results_perf_all, "results/performance_stratified_lmm.csv")
cat("\nSaved: results/performance_stratified_lmm.csv\n")

################################################################################
# 6. AGE × PERFORMANCE STRATIFIED MODELS (2-way)
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("AGE × PERFORMANCE STRATIFIED LINEAR MIXED MODELS\n")
cat(rep("=", 80), "\n\n", sep="")

# Function to fit age × performance stratified models
fit_age_performance_stratified <- function(data, pollutant_col, pollutant_name) {
  cat(sprintf("\n--- %s ---\n", pollutant_name))

  results_list <- list()

  for (age_grp in levels(data$age_group)) {
    for (perf_grp in levels(data$performance_group)) {

      data_subset <- data[age_group == age_grp & performance_group == perf_grp &
                          !is.na(get(pollutant_col)) & !is.na(heat_index) & !is.na(wind_speed_race_day) & !is.na(sex)]

      # Only fit if sufficient sample size
      if (nrow(data_subset) < 500) next

      cat(sprintf("\n  %s × %s: N=%s\n", age_grp, perf_grp,
                  format(nrow(data_subset), big.mark=",")))

      # Pooled model (both sexes)
      formula_str <- sprintf("time_seconds ~ %s + sex + year_factor + ns(heat_index, df=5) + wind_speed_race_day + (1|event)",
                            pollutant_col)

      tryCatch({
        model <- lmer(as.formula(formula_str), data = data_subset, REML = TRUE)
        coef_val <- fixef(model)[pollutant_col]
        ci <- confint(model, parm = pollutant_col, method = "profile", quiet = TRUE)

        results_list[[paste(age_grp, perf_grp, sep="_")]] <- data.table(
          pollutant = pollutant_name,
          age_group = age_grp,
          performance_group = perf_grp,
          coefficient = coef_val,
          ci_lower = ci[1],
          ci_upper = ci[2],
          n = nrow(data_subset)
        )

        cat(sprintf("    Coefficient: %.2f sec/µg/m³ [%.2f, %.2f]\n",
                    coef_val, ci[1], ci[2]))
      }, error = function(e) {
        cat(sprintf("    ERROR: %s\n", e$message))
      })
    }
  }

  results_df <- rbindlist(results_list, fill = TRUE)
  return(results_df)
}

# Complete cases
df_no2_ap <- df[!is.na(no2_race_day) & !is.na(age_group) & !is.na(performance_group)]
df_pm25_ap <- df[!is.na(pm25_race_day) & !is.na(age_group) & !is.na(performance_group)]

cat("\nNote: Only fitting cells with N >= 500 to ensure stable estimates\n")

# Fit age × performance stratified models
results_ap_no2 <- fit_age_performance_stratified(df_no2_ap, "no2_race_day", "NO2")
results_ap_pm25 <- fit_age_performance_stratified(df_pm25_ap, "pm25_race_day", "PM2.5")

# Combine and save
results_ap_all <- rbind(results_ap_no2, results_ap_pm25)
fwrite(results_ap_all, "results/age_performance_stratified_lmm.csv")
cat("\nSaved: results/age_performance_stratified_lmm.csv\n")

################################################################################
# 7. VISUALIZATIONS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("Creating visualizations...\n")
cat(rep("=", 80), "\n\n", sep="")

# 1. Age-stratified forest plot
if (nrow(results_age_all) > 0) {
  results_age_all[, coef_min := coefficient / 60]
  results_age_all[, ci_lower_min := ci_lower / 60]
  results_age_all[, ci_upper_min := ci_upper / 60]

  p1 <- ggplot(results_age_all, aes(x = age_group, y = coef_min, color = pollutant)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = ci_lower_min, ymax = ci_upper_min),
                  width = 0.3, linewidth = 1, position = position_dodge(0.5)) +
    geom_point(size = 3, position = position_dodge(0.5)) +
    scale_color_manual(values = c("NO2" = "#d62728", "PM2.5" = "#ff7f0e")) +
    labs(
      title = "Pollution Effects on Marathon Performance by Age Group",
      subtitle = "Linear mixed model coefficients with 95% confidence intervals",
      x = "Age Group",
      y = "Coefficient (minutes per µg/m³)",
      color = "Pollutant"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave("plots/age_stratified_forest.png", p1, width = 9, height = 6, dpi = 300)
  cat("Saved: plots/age_stratified_forest.png\n")
}

# 2. Performance-stratified forest plot
if (nrow(results_perf_all) > 0) {
  results_perf_all[, coef_min := coefficient / 60]
  results_perf_all[, ci_lower_min := ci_lower / 60]
  results_perf_all[, ci_upper_min := ci_upper / 60]

  p2 <- ggplot(results_perf_all, aes(x = performance_group, y = coef_min, color = pollutant)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbar(aes(ymin = ci_lower_min, ymax = ci_upper_min),
                  width = 0.3, linewidth = 1, position = position_dodge(0.5)) +
    geom_point(size = 3, position = position_dodge(0.5)) +
    scale_color_manual(values = c("NO2" = "#d62728", "PM2.5" = "#ff7f0e")) +
    labs(
      title = "Pollution Effects on Marathon Performance by Performance Level",
      subtitle = "Linear mixed model coefficients with 95% confidence intervals",
      x = "Performance Level",
      y = "Coefficient (minutes per µg/m³)",
      color = "Pollutant"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  ggsave("plots/performance_stratified_forest.png", p2, width = 10, height = 6, dpi = 300)
  cat("Saved: plots/performance_stratified_forest.png\n")
}

# 3. Age × Performance heatmap
if (nrow(results_ap_all) > 0) {
  results_ap_all[, coef_min := coefficient / 60]

  for (poll in unique(results_ap_all$pollutant)) {
    data_plot <- results_ap_all[pollutant == poll]

    p3 <- ggplot(data_plot, aes(x = performance_group, y = age_group, fill = coef_min)) +
      geom_tile(color = "white", linewidth = 1) +
      geom_text(aes(label = sprintf("%.2f", coef_min)), color = "white", fontface = "bold") +
      scale_fill_gradient2(
        low = "#2166ac", mid = "white", high = "#b2182b",
        midpoint = 0, name = "Coefficient\n(min/µg/m³)"
      ) +
      labs(
        title = sprintf("%s Effects: Age × Performance Interaction", poll),
        subtitle = "Linear mixed model coefficients (minutes per µg/m³)",
        x = "Performance Level",
        y = "Age Group"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    ggsave(sprintf("plots/age_performance_heatmap_%s.png", tolower(poll)),
           p3, width = 11, height = 7, dpi = 300)
    cat(sprintf("Saved: plots/age_performance_heatmap_%s.png\n", tolower(poll)))
  }
}

################################################################################
# 8. SUMMARY STATISTICS
################################################################################

cat("\n", rep("=", 80), "\n", sep="")
cat("SUMMARY STATISTICS\n")
cat(rep("=", 80), "\n\n", sep="")

# Age group summary
cat("Age-stratified results (NO2):\n")
if (nrow(results_age_no2) > 0) {
  print(results_age_no2[, .(age_group, coefficient, ci_lower, ci_upper, n)])
}

cat("\nPerformance-stratified results (NO2):\n")
if (nrow(results_perf_no2) > 0) {
  print(results_perf_no2[, .(performance_group, coefficient, ci_lower, ci_upper, n)])
}

cat("\nAge × Performance results (NO2, top 10 largest effects):\n")
if (nrow(results_ap_no2) > 0) {
  top_effects <- results_ap_no2[order(-coefficient)][1:min(10, .N)]
  print(top_effects[, .(age_group, performance_group, coefficient, ci_lower, ci_upper, n)])
}

################################################################################
# 9. GENERATE SUMMARY REPORT
################################################################################

sink("results/age_ability_subgroup_report.txt")

cat(rep("=", 100), "\n")
cat("AGE AND ABILITY SUBGROUP ANALYSIS REPORT\n")
cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 100), "\n\n")

cat("OVERVIEW\n")
cat(rep("-", 100), "\n")
cat("This analysis explores heterogeneous pollution effects by age and performance level\n")
cat("using the same linear mixed model framework as the main analysis.\n\n")

cat("SAMPLE SIZE BY SUBGROUP\n")
cat(rep("-", 100), "\n")
cat("\nAge groups:\n")
print(df[, .N, by = age_group][order(age_group)])
cat("\nPerformance levels:\n")
print(df[, .N, by = performance_group][order(performance_group)])

cat("\n\nAGE-STRATIFIED RESULTS\n")
cat(rep("-", 100), "\n")
if (nrow(results_age_all) > 0) {
  print(results_age_all)
} else {
  cat("No results available\n")
}

cat("\n\nPERFORMANCE-STRATIFIED RESULTS\n")
cat(rep("-", 100), "\n")
if (nrow(results_perf_all) > 0) {
  print(results_perf_all)
} else {
  cat("No results available\n")
}

cat("\n\nAGE × PERFORMANCE STRATIFIED RESULTS (TOP 20 EFFECTS)\n")
cat(rep("-", 100), "\n")
if (nrow(results_ap_all) > 0) {
  top20 <- results_ap_all[order(-coefficient)][1:min(20, .N)]
  print(top20)
} else {
  cat("No results available\n")
}

cat("\n\nKEY FINDINGS\n")
cat(rep("-", 100), "\n")
cat("1. Age effects: See age-stratified results table above\n")
cat("2. Performance effects: See performance-stratified results table above\n")
cat("3. Interaction effects: See age × performance heatmaps\n")
cat("4. All models control for heat index, year, and event-level clustering\n")

cat("\n\nFILES GENERATED\n")
cat(rep("-", 100), "\n")
cat("Results:\n")
cat("  - results/age_stratified_lmm.csv\n")
cat("  - results/performance_stratified_lmm.csv\n")
cat("  - results/age_performance_stratified_lmm.csv\n")
cat("  - results/age_ability_subgroup_report.txt (this file)\n\n")
cat("Plots:\n")
cat("  - plots/age_stratified_forest.png\n")
cat("  - plots/performance_stratified_forest.png\n")
cat("  - plots/age_performance_heatmap_no2.png\n")
cat("  - plots/age_performance_heatmap_pm2.5.png\n")

cat("\n", rep("=", 100), "\n")
cat("END OF REPORT\n")
cat(rep("=", 100), "\n")

sink()

cat("\nSummary report saved to: results/age_ability_subgroup_report.txt\n")

cat("\n", rep("=", 80), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 80), "\n", sep="")
cat("\nResults saved to results/ directory\n")
cat("Plots saved to plots/ directory\n")
