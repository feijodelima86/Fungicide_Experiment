library(readr)
library(lme4)
library(lmerTest)  # For p-values in mixed models
library(emmeans)   # For post-hoc tests
library(ggplot2)   # For improved plotting
library(gridExtra) # For arranging multiple plots
#library(AICcmodavg) # For model comparison
library(sjPlot)    # For model visualization

# ==== Model Comparison Framework for Standard vs Mixed Effects ====
# This script creates and compares both model types for the fungicide experiment

# Read the data
data <- read.csv("0_Data/Data_Merged_OLRM_renamed.csv")


# Map group numbers to treatment names
treatment_map <- c(
  "1" = "Control1",    # Nothing Added
  "2" = "Control2",    # 100% acetone, No fungicide
  "3" = "Treatment1",  # 0.01% acetone and 0.01% fungicide
  "4" = "Treatment2",  # 0.1% acetone and 0.1% fungicide
  "5" = "Treatment3",  # 1% acetone and 1% fungicide
  "6" = "Treatment4",  # 10% acetone and 10% fungicide
  "7" = "Treatment5"   # 100% acetone and 100% fungicide
)

# Add a treatment column
data$Treatment <- treatment_map[as.character(data$Group)]

# Convert Day to numeric to ensure proper treatment in models
data$Day <- as.numeric(as.character(data$Day))

# Create a unique Sample identifier based on the bucket ID
# If Sample already identifies the unique buckets, we can use it directly
if(!"Sample.ID" %in% names(data)) {
  # Create a unique identifier that distinguishes each bucket
  data$Sample.ID <- paste(data$Sample, sep="_")
}

# Explicitly set treatment as a factor with Control1 as reference level
data$Treatment <- factor(data$Treatment, 
                         levels = c("Control1", "Control2", 
                                    "Treatment1", "Treatment2", 
                                    "Treatment3", "Treatment4", 
                                    "Treatment5"))

# Print basic dataset information
cat("===== Dataset Information =====\n")
print(paste("Number of observations:", nrow(data)))
print(paste("Number of unique samples:", length(unique(data$Sample.ID))))
print("Sample distribution by Treatment and Day:")
print(table(data$Treatment, data$Day))

# Check for missing values in the response variable
na_count <- sum(is.na(data$AI))
if(na_count > 0) {
  cat("\nWARNING:", na_count, "missing values found in AI\n")
  # Optionally filter out NA values
  data <- data[!is.na(data$AI),]
  cat("Filtered dataset now has", nrow(data), "observations\n")
}

# ===== PART 1: Estimate Acetone Effect =====
cat("\n===== PART 1: Estimate Acetone Effect =====\n")

# Filter to control groups
controls_data <- data[data$Treatment %in% c("Control1", "Control2"), ]
controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)

# 1A. Standard Linear Model for Acetone Effect
std_acetone_model <- lm(AI ~ Day + is_acetone, data = controls_data)
std_acetone_effect <- coef(std_acetone_model)["is_acetone"]

# 1B. Mixed Effects Model for Acetone Effect
mixed_acetone_model <- lmer(AI ~ Day + is_acetone + (1|Sample.ID), data = controls_data)
mixed_acetone_effect <- fixef(mixed_acetone_model)["is_acetone"]

# Print model summaries
cat("\n-- Standard Linear Model for Acetone Effect --\n")
print(summary(std_acetone_model))
cat("\nAcetone Effect (Standard Model):", std_acetone_effect, "\n")

cat("\n-- Mixed Effects Model for Acetone Effect --\n")
print(summary(mixed_acetone_model))
cat("\nAcetone Effect (Mixed Model):", mixed_acetone_effect, "\n")

# ===== PART 2: Create Adjusted Datasets =====
cat("\n===== PART 2: Creating Adjusted Datasets =====\n")

# The acetone concentrations for each treatment
acetone_conc <- c(
  "Control1" = 0,
  "Control2" = 100,
  "Treatment1" = 0.01,
  "Treatment2" = 0.1,
  "Treatment3" = 1, 
  "Treatment4" = 10,
  "Treatment5" = 100
)

# Calculate expected acetone effects under both models
data$std_expected_acetone_effect <- (acetone_conc[data$Treatment] / 100) * std_acetone_effect
data$mixed_expected_acetone_effect <- (acetone_conc[data$Treatment] / 100) * mixed_acetone_effect

# Create adjusted response variables
data$std_adjusted_AI <- data$AI - data$std_expected_acetone_effect
data$mixed_adjusted_AI <- data$AI - data$mixed_expected_acetone_effect

# Compare the adjustments
cat("\nComparison of Acetone Effect Adjustments by Treatment:\n")
adj_comparison <- aggregate(cbind(std_expected_acetone_effect, mixed_expected_acetone_effect) ~ Treatment, 
                            data = data, mean)
print(adj_comparison)

# ===== PART 3: Full Treatment Models =====
cat("\n===== PART 3: Full Treatment Models =====\n")

# 3A. Standard Linear Models
std_original_model <- lm(AI ~ Day + Treatment, data = data)
std_adjusted_model <- lm(std_adjusted_AI ~ Day + Treatment, data = data)

# 3B. Mixed Effects Models
mixed_original_model <- lmer(AI ~ Day + Treatment + (1|Sample.ID), data = data)
mixed_adjusted_model <- lmer(mixed_adjusted_AI ~ Day + Treatment + (1|Sample.ID), data = data)

# Print model summaries
cat("\n-- Standard Linear Model (Original AI) --\n")
print(summary(std_original_model))

cat("\n-- Standard Linear Model (Acetone-Adjusted AI) --\n")
print(summary(std_adjusted_model))

cat("\n-- Mixed Effects Model (Original AI) --\n")
print(summary(mixed_original_model))

cat("\n-- Mixed Effects Model (Acetone-Adjusted AI) --\n")
print(summary(mixed_adjusted_model))

# ===== PART 4: Model Comparison =====
cat("\n===== PART 4: Model Comparison =====\n")

# 4A. Compare information criteria
models <- list(
  StdOriginal = std_original_model,
  StdAdjusted = std_adjusted_model,
  MixedOriginal = mixed_original_model,
  MixedAdjusted = mixed_adjusted_model
)

# AIC and BIC comparison
aic_values <- sapply(models, AIC)
bic_values <- sapply(models, BIC)

model_comparison <- data.frame(
  Model = names(models),
  AIC = aic_values,
  BIC = bic_values
)
print(model_comparison)

# 4B. Compare ANOVA results
cat("\n-- ANOVA Results from Standard Model (Adjusted AI) --\n")
std_anova <- anova(std_adjusted_model)
print(std_anova)

cat("\n-- ANOVA Results from Mixed Model (Adjusted AI) --\n")
mixed_anova <- anova(mixed_adjusted_model)
print(mixed_anova)

# 4C. Compare treatment effects
# Extract coefficients and standard errors
std_coefs <- coef(summary(std_adjusted_model))
mixed_coefs <- cbind(fixef(mixed_adjusted_model), 
                     sqrt(diag(vcov(mixed_adjusted_model))))

# Create comparison table
coef_comparison <- data.frame(
  Term = rownames(std_coefs),
  Std_Estimate = std_coefs[, "Estimate"],
  Std_StdError = std_coefs[, "Std. Error"],
  Std_tvalue = std_coefs[, "t value"],
  Std_pvalue = std_coefs[, "Pr(>|t|)"]
)

mixed_terms <- names(fixef(mixed_adjusted_model))
coef_comparison$Mixed_Estimate <- NA
coef_comparison$Mixed_StdError <- NA

for(i in 1:length(mixed_terms)) {
  idx <- which(coef_comparison$Term == mixed_terms[i])
  if(length(idx) > 0) {
    coef_comparison$Mixed_Estimate[idx] <- mixed_coefs[i, 1]
    coef_comparison$Mixed_StdError[idx] <- mixed_coefs[i, 2]
  }
}

cat("\n-- Coefficient Comparison between Standard and Mixed Models --\n")
print(coef_comparison)

# 4D. Compare post-hoc tests
cat("\n-- Post-Hoc Comparisons: Standard Model --\n")
std_emmeans <- emmeans(std_adjusted_model, ~ Treatment)
std_pairs <- pairs(std_emmeans)
print(std_pairs)

cat("\n-- Post-Hoc Comparisons: Mixed Model --\n")
mixed_emmeans <- emmeans(mixed_adjusted_model, ~ Treatment)
mixed_pairs <- pairs(mixed_emmeans)
print(mixed_pairs)

# ===== PART 5: Data Visualization =====
cat("\n===== PART 5: Data Visualization =====\n")

# Function to calculate means and SEs for a given variable
calculate_stats <- function(data, variable, use_mixed = FALSE) {
  # For standard calculation
  if(!use_mixed) {
    means <- aggregate(data[[variable]] ~ Treatment + Day, data = data, FUN = mean)
    names(means)[3] <- "mean"
    
    ses <- aggregate(data[[variable]] ~ Treatment + Day, data = data, 
                     FUN = function(x) sd(x) / sqrt(length(x)))
    names(ses)[3] <- "se"
  } else {
    # For calculation accounting for repeated measures
    means <- aggregate(data[[variable]] ~ Treatment + Day, data = data, FUN = mean)
    names(means)[3] <- "mean"
    
    # Calculate SEs accounting for repeated measures
    ses <- aggregate(data[[variable]] ~ Treatment + Day, data = data, 
                     FUN = function(x) {
                       # Get the unique samples for this subset
                       samples <- unique(data$Sample.ID[data[[variable]] %in% x])
                       # Calculate SE based on number of unique samples
                       return(sd(x) / sqrt(length(samples)))
                     })
    names(ses)[3] <- "se"
  }
  
  result <- merge(means, ses, by = c("Treatment", "Day"))
  return(result)
}

# Calculate statistics for original and both adjusted AI values
std_mean_AI <- calculate_stats(data, "AI", FALSE)
mixed_mean_AI <- calculate_stats(data, "AI", TRUE)
std_mean_adjusted <- calculate_stats(data, "std_adjusted_AI", FALSE)
mixed_mean_adjusted <- calculate_stats(data, "mixed_adjusted_AI", TRUE)

# For ggplot - ensure treatments are ordered correctly
std_mean_AI$Treatment <- factor(std_mean_AI$Treatment, 
                                  levels = c("Control1", "Control2", 
                                             "Treatment1", "Treatment2", 
                                             "Treatment3", "Treatment4", 
                                             "Treatment5"))

mixed_mean_AI$Treatment <- factor(mixed_mean_AI$Treatment,
                                    levels = c("Control1", "Control2", 
                                               "Treatment1", "Treatment2", 
                                               "Treatment3", "Treatment4", 
                                               "Treatment5"))

std_mean_adjusted$Treatment <- factor(std_mean_adjusted$Treatment,
                                      levels = c("Control1", "Control2", 
                                                 "Treatment1", "Treatment2", 
                                                 "Treatment3", "Treatment4", 
                                                 "Treatment5"))

mixed_mean_adjusted$Treatment <- factor(mixed_mean_adjusted$Treatment,
                                        levels = c("Control1", "Control2", 
                                                   "Treatment1", "Treatment2", 
                                                   "Treatment3", "Treatment4", 
                                                   "Treatment5"))

# Create ggplot objects
plot_std_original <- ggplot(std_mean_AI, aes(x = Day, y = mean, color = Treatment, shape = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Original AI (Standard SE)", y = "AI") +
  theme_bw() +
  theme(legend.position = "right")

plot_mixed_original <- ggplot(mixed_mean_AI, aes(x = Day, y = mean, color = Treatment, shape = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Original AI (Mixed Model SE)", y = "AI") +
  theme_bw() +
  theme(legend.position = "right")

plot_std_adjusted <- ggplot(std_mean_adjusted, aes(x = Day, y = mean, color = Treatment, shape = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Acetone-Adjusted AI (Standard Model)", y = "Adjusted AI") +
  theme_bw() +
  theme(legend.position = "right")

plot_mixed_adjusted <- ggplot(mixed_mean_adjusted, aes(x = Day, y = mean, color = Treatment, shape = Treatment)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  labs(title = "Acetone-Adjusted AI (Mixed Model)", y = "Adjusted AI") +
  theme_bw() +
  theme(legend.position = "right")

# Arrange and display plots
grid.arrange(plot_std_original, plot_mixed_original, 
             plot_std_adjusted, plot_mixed_adjusted,
             nrow = 2, ncol = 2)

# Save plots if needed
ggsave("model_comparison_plots.pdf", 
       grid.arrange(plot_std_original, plot_mixed_original, 
                    plot_std_adjusted, plot_mixed_adjusted,
                    nrow = 2, ncol = 2),
       width = 12, height = 10)

# ===== PART 6: Random Effects Visualization =====
cat("\n===== PART 6: Random Effects Analysis =====\n")

# Extract random effects
random_effects <- ranef(mixed_adjusted_model)$Sample.ID
random_effects_df <- data.frame(
  Sample.ID = rownames(random_effects),
  RandomEffect = random_effects[,1]
)

# Calculate percentiles for interpretation
quantiles <- quantile(random_effects_df$RandomEffect, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

cat("\nRandom Effects Summary:\n")
print(summary(random_effects_df$RandomEffect))
cat("\nRandom Effects Quantiles:\n")
print(quantiles)

# Visualize random effects
ggplot(random_effects_df, aes(x = RandomEffect)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Random Effects by Sample", 
       x = "Random Effect Value", y = "Count") +
  theme_bw()
ggsave("random_effects_distribution.pdf", width = 8, height = 6)

# Visualize random effects by model predictions
plot_ranef(mixed_adjusted_model)
ggsave("random_effects_visualization.pdf", width = 10, height = 8)

# ===== PART 7: Diagnostics =====
cat("\n===== PART 7: Model Diagnostics =====\n")

# Standard model diagnostics
par(mfrow = c(2, 2))
plot(std_adjusted_model, main = "Standard Model Diagnostics")

# Mixed model diagnostics
# Create diagnostic plots for the mixed model
plot_model(mixed_adjusted_model, type = "diag")
ggsave("mixed_model_diagnostics.pdf", width = 10, height = 8)

# QQ plots of residuals
std_resid <- residuals(std_adjusted_model)
mixed_resid <- residuals(mixed_adjusted_model)

par(mfrow = c(1, 2))
qqnorm(std_resid, main = "QQ Plot - Standard Model")
qqline(std_resid)
qqnorm(mixed_resid, main = "QQ Plot - Mixed Model")
qqline(mixed_resid)

# ===== PART 8: Summary Report =====
cat("\n===== PART 8: Summary Report =====\n")

# Summarize key findings
cat("\nKey Statistical Findings:\n")
cat("1. Acetone Effect Estimates:\n")
cat("   - Standard Model:", std_acetone_effect, "\n")
cat("   - Mixed Model:", mixed_acetone_effect, "\n")
cat("   - Difference:", mixed_acetone_effect - std_acetone_effect, "\n")

cat("\n2. Model Fit Comparison (AIC, lower is better):\n")
cat("   - Standard Original Model:", AIC(std_original_model), "\n")
cat("   - Standard Adjusted Model:", AIC(std_adjusted_model), "\n") 
cat("   - Mixed Original Model:", AIC(mixed_original_model), "\n")
cat("   - Mixed Adjusted Model:", AIC(mixed_adjusted_model), "\n")

cat("\n3. Treatment Effects (p-values):\n")
cat("   - Standard Model:", std_anova$`Pr(>F)`[2], "\n")
cat("   - Mixed Model:", mixed_anova$`Pr(>F)`[2], "\n")

cat("\n4. Day Effects (p-values):\n")
cat("   - Standard Model:", std_anova$`Pr(>F)`[1], "\n")
cat("   - Mixed Model:", mixed_anova$`Pr(>F)`[1], "\n")

cat("\n5. Residual Standard Error:\n")
cat("   - Standard Model:", sigma(std_adjusted_model), "\n")
cat("   - Mixed Model (Residual):", attr(VarCorr(mixed_adjusted_model), "sc"), "\n")

cat("\n6. Random Effect Standard Deviation:\n")
cat("   - Mixed Model:", as.data.frame(VarCorr(mixed_adjusted_model))$sdcor[1], "\n")

cat("\n7. Significant Pairwise Differences (p < 0.05):\n")
std_sig <- summary(std_pairs)
std_sig_pairs <- std_sig[std_sig$p.value < 0.05, ]
mixed_sig <- summary(mixed_pairs)
mixed_sig_pairs <- mixed_sig[mixed_sig$p.value < 0.05, ]

cat("   - Standard Model Significant Pairs:", nrow(std_sig_pairs), "\n")
cat("   - Mixed Model Significant Pairs:", nrow(mixed_sig_pairs), "\n")

cat("\n===== Conclusion =====\n")
cat("This analysis compared standard linear models and mixed-effects models\n")
cat("for analyzing the fungicide experiment data. The mixed-effects approach\n")
cat("properly accounts for the repeated measurements from the same buckets over time,\n")
cat("which can lead to more accurate estimates of treatment effects and standard errors.\n")

cat("\nThe key differences between the approaches are:\n")
cat("1. Standard Errors: Mixed models account for correlation within samples, often\n")
cat("   resulting in larger, more conservative standard errors\n")
cat("2. Model Fit: Mixed models typically show better fit according to AIC/BIC\n")
cat("3. P-values and Statistical Significance: Differences in significance can\n")
cat("   occur between the two approaches due to the accounting for repeated measures\n")

cat("\nBased on the results, the recommended model for this data is the\n")
if(AIC(mixed_adjusted_model) < AIC(std_adjusted_model)) {
  cat("mixed-effects model with acetone adjustment, which shows better fit and\n")
  cat("accounts for the experimental design with repeated measurements.\n")
} else {
  cat("standard model with acetone adjustment, which shows better fit despite\n")
  cat("not accounting for repeated measurements.\n")
}

