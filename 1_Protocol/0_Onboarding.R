library(readr)


# Analysis of C/P Ratio with Acetone Adjustment

# Read the data
data <- read.csv(text = 'Sample,Group,Day,Fungi_Est,Acet_Est,CP_ratio
a1,1,3,0,0,0.133333
a2,1,3,0,0,0.2
a3,1,3,0,0,0.1875
b1,2,3,0,100,0.157895
b7,2,3,0,100,0.190476
c4,2,3,0,100,0.166667
b2,3,3,0.01,0.01,0.166667
b8,3,3,0.01,0.01,0.217391
c5,3,3,0.01,0.01,0.153846
b3,4,3,0.1,0.1,0.15
b9,4,3,0.1,0.1,0.190476
c6,4,3,0.1,0.1,0.176471
b4,5,3,1,1,0.166667
c1,5,3,1,1,0.157895
c7,5,3,1,1,0.2
b5,6,3,10,10,0.142857
c2,6,3,10,10,0.25
c8,6,3,10,10,0.157895
b6,7,3,100,100,0.136364
c3,7,3,100,100,0.166667
c9,7,3,100,100,0.176471
a1,1,5,0,0,0.15
a2,1,5,0,0,0.181818
a3,1,5,0,0,0.111111
b1,2,5,0,100,0.2
b7,2,5,0,100,0.235294
c4,2,5,0,100,0.263158
b2,3,5,0.01,0.01,0.227273
b8,3,5,0.01,0.01,0.190476
c5,3,5,0.01,0.01,0.222222
b3,4,5,0.1,0.1,0.176471
b9,4,5,0.1,0.1,0.208333
c6,4,5,0.1,0.1,0.185185
b4,5,5,1,1,0.176471
c1,5,5,1,1,0.2
c7,5,5,1,1,0.181818
b5,6,5,10,10,0.157895
c2,6,5,10,10,0.16
c8,6,5,10,10,0.176471
b6,7,5,100,100,0.133333
c3,7,5,100,100,0.125
c9,7,5,100,100,0.210526
a1,1,7,0,0,0.25
a2,1,7,0,0,0.25
a3,1,7,0,0,0.238095
b1,2,7,0,100,0.266667
b7,2,7,0,100,0.3
c4,2,7,0,100,0.333333
b2,3,7,0.01,0.01,0.272727
b8,3,7,0.01,0.01,0.25
c5,3,7,0.01,0.01,0.294118
b3,4,7,0.1,0.1,0.166667
b9,4,7,0.1,0.1,0.1875
c6,4,7,0.1,0.1,0.2
b4,5,7,1,1,0.2
c1,5,7,1,1,0.157895
c7,5,7,1,1,0.166667
b5,6,7,10,10,0.166667
c2,6,7,10,10,0.173913
c8,6,7,10,10,0.166667
b6,7,7,100,100,0.133333
c3,7,7,100,100,0.157895
c9,7,7,100,100,0.181818
a1,1,10,0,0,0.3
a2,1,10,0,0,0.333333
a3,1,10,0,0,0.285714
b1,2,10,0,100,0.384615
b7,2,10,0,100,0.333333
c4,2,10,0,100,0.357143
b2,3,10,0.01,0.01,0.285714
b8,3,10,0.01,0.01,0.285714
c5,3,10,0.01,0.01,0.307692
b3,4,10,0.1,0.1,0.285714
b9,4,10,0.1,0.1,0.285714
c6,4,10,0.1,0.1,0.277778
b4,5,10,1,1,0.363636
c1,5,10,1,1,0.25
c7,5,10,1,1,0.181818
b5,6,10,10,10,0.307692
c2,6,10,10,10,0.181818
c8,6,10,10,10,0.263158
b6,7,10,100,100,0.384615
c3,7,10,100,100,0.272727
c9,7,10,100,100,0.266667')


data

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

# Verify the loaded data
head(data)
cat("\nData structure:\n")
str(data)
cat("\nSummary of C/P ratio:\n")
summary(data$CP_ratio)

# Check the number of samples per treatment and day
cat("\nSamples per treatment and day:\n")
print(table(data$Treatment, data$Day))

# Step 1: Estimate acetone effect from controls
# Run GLM with time as covariate comparing Control 1 vs Control 2
controls_data <- data[data$Treatment %in% c("Control1", "Control2"), ]
controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)

acetone_model <- glm(CP_ratio ~ Day + is_acetone, 
                     data = controls_data, 
                     family = gaussian())

plot(acetone_model)

cat("\nAcetone effect model summary:\n")
print(summary(acetone_model))

# Extract the acetone effect (per 100% concentration)
acetone_effect_estimate <- coef(acetone_model)["is_acetone"]
cat("\nEstimated effect of 100% acetone:", acetone_effect_estimate, "\n")

# Step 2: Create proportional acetone effects
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

# Calculate expected acetone effect for each treatment based on concentration
data$expected_acetone_effect <- (acetone_conc[data$Treatment] / 100) * acetone_effect_estimate
data$acetone_adjusted_CP_ratio <- data$CP_ratio - data$expected_acetone_effect



# Step 3: Run analysis on acetone-adjusted responses
treatment_model <- glm(acetone_adjusted_CP_ratio ~ Day + Treatment + Acet_Est, 
                       data = data, 
                       family = gaussian())

cat("\nTreatment model (with acetone adjustment) summary:\n")
print(summary(treatment_model))

# Calculate means and standard errors for original and adjusted CP ratios
mean_CP_ratio <- aggregate(CP_ratio ~ Treatment + Day, data = data, 
                           FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_CP_ratio <- do.call(data.frame, mean_CP_ratio)
names(mean_CP_ratio)[3:4] <- c("mean_CP_ratio", "se_CP_ratio")

mean_adjusted <- aggregate(acetone_adjusted_CP_ratio ~ Treatment + Day, data = data,
                           FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_adjusted <- do.call(data.frame, mean_adjusted)
names(mean_adjusted)[3:4] <- c("mean_adjusted", "se_adjusted")

# View the means
cat("\nOriginal CP Ratio means (sample):\n")
print(head(mean_CP_ratio))
cat("\nAcetone-adjusted CP Ratio means (sample):\n")
print(head(mean_adjusted))

# For plotting - ensure treatments are ordered correctly
mean_CP_ratio$Treatment <- factor(mean_CP_ratio$Treatment, 
                                  levels = c("Control1", "Control2", 
                                             "Treatment1", "Treatment2", 
                                             "Treatment3", "Treatment4", 
                                             "Treatment5"))

mean_adjusted$Treatment <- factor(mean_adjusted$Treatment,
                                  levels = c("Control1", "Control2", 
                                             "Treatment1", "Treatment2", 
                                             "Treatment3", "Treatment4", 
                                             "Treatment5"))

# Base R plotting code for both plots
# Define a color palette for the treatments
treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)  # Different point symbols

# Set up plotting layout - two plots stacked vertically
par(mfrow=c(2,1))

# Function to create the plot with error bars
create_treatment_plot <- function(data, y_var, se_var, y_label, main_title) {
  # Get the unique treatments and days
  treatments <- levels(data$Treatment)
  days <- sort(unique(data$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data[[y_var]] - data[[se_var]], na.rm = TRUE)
  y_max <- max(data[[y_var]] + data[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(y_min, y_max * 1.1),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Plot each treatment
  for (i in 1:length(treatments)) {
    treatment <- treatments[i]
    treatment_data <- data[data$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = 2)
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = treatment_colors[i], lwd = 1.5)
    }
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.3, 0), 
         legend = treatments, 
         col = treatment_colors[1:length(treatments)],
         pch = treatment_pch[1:length(treatments)],
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Plot 1: Original CP Ratio
create_treatment_plot(mean_CP_ratio, "mean_CP_ratio", "se_CP_ratio", 
                      "C/P Ratio", "Original C/P Ratio by Treatment and Day")

# Plot 2: Acetone-Adjusted CP Ratio
create_treatment_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                      "Adjusted C/P Ratio", "Acetone-Adjusted C/P Ratio by Treatment and Day")



# Step 4: Compare to traditional ANCOVA and analyze dose-response relationship

# Run traditional ANCOVA (without acetone adjustment)
ancova_model <- glm(CP_ratio ~ Day + Treatment+ Acet_Est, 
                    data = data, 
                    family = gaussian())

cat("\nTraditional ANCOVA model summary:\n")
print(summary(ancova_model))

# Compare AIC between approaches
cat("\nModel comparison:\n")
cat("AIC for traditional ANCOVA:", AIC(ancova_model), "\n")
cat("AIC for treatment model with acetone adjustment:", AIC(treatment_model), "\n")

# Optional: Analyze dose-response relationship for fungicide after acetone adjustment
# First define fungicide concentrations for each treatment
fungicide_conc <- c(
  "Control1" = 0,
  "Control2" = 0,
  "Treatment1" = 0.01,
  "Treatment2" = 0.1,
  "Treatment3" = 1, 
  "Treatment4" = 10,
  "Treatment5" = 100
)

# Add fungicide concentration to the data
data$Fungi_Conc <- fungicide_conc[data$Treatment]

# Filter out only treatments with fungicide (not controls)
fungicide_data <- data[data$Fungi_Conc > 0, ]

# Log-transform fungicide concentration for dose-response analysis
fungicide_data$log_fungi <- log10(fungicide_data$Fungi_Conc)

# Model fungicide effect with log-transformed concentration
fungicide_model <- glm(acetone_adjusted_CP_ratio ~ Day + log_fungi,
                       data = fungicide_data,
                       family = gaussian())

cat("\nFungicide dose-response model (log-transformed) summary:\n")
print(summary(fungicide_model))

# Create dose-response visualization for the final day
final_day_data <- fungicide_data[fungicide_data$Day == max(fungicide_data$Day), ]

# Calculate means and standard errors by fungicide concentration for the final day
dose_response <- aggregate(acetone_adjusted_CP_ratio ~ Fungi_Conc, data = final_day_data,
                           FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
dose_response <- do.call(data.frame, dose_response)
names(dose_response)[2:3] <- c("mean_response", "se_response")

# Add log-transformed concentration for plotting
dose_response$log_conc <- log10(dose_response$Fungi_Conc)

# Plot dose-response curve for the final day
par(mfrow=c(1,1))  # Reset to single plot
par(mar=c(5, 5, 4, 2))  # Adjust margins

# Create the dose-response plot
plot(dose_response$log_conc, dose_response$mean_response,
     type = "o", pch = 16, col = "blue", lwd = 2,
     xlab = "Fungicide Concentration (log10 %)", 
     ylab = "Adjusted C/P Ratio",
     main = paste("Dose-Response Curve (Day", max(fungicide_data$Day), ")"),
     xaxt = "n")  # Custom x-axis labels

# Create custom x-axis with actual concentration values
axis(1, at = log10(unique(fungicide_data$Fungi_Conc)), 
     labels = unique(fungicide_data$Fungi_Conc))

# Add error bars
arrows(dose_response$log_conc, 
       dose_response$mean_response - dose_response$se_response,
       dose_response$log_conc, 
       dose_response$mean_response + dose_response$se_response,
       angle = 90, code = 3, length = 0.05, col = "blue", lwd = 1.5)

# Add grid lines
grid(lty = 2, col = "gray80")

# Optional: Add regression line if the fungicide effect is significant
if (summary(fungicide_model)$coefficients["log_fungi", "Pr(>|t|)"] < 0.05) {
  # Create sequence of x values for prediction
  x_seq <- seq(min(dose_response$log_conc), max(dose_response$log_conc), length.out = 100)
  
  # Create prediction data frame
  pred_data <- data.frame(
    log_fungi = x_seq,
    Day = max(fungicide_data$Day)  # Use the final day for prediction
  )
  
  # Make predictions
  pred_y <- predict(fungicide_model, newdata = pred_data, type = "response")
  
  # Add regression line
  lines(x_seq, pred_y, col = "red", lwd = 2, lty = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Observed", "Predicted"), 
         col = c("blue", "red"),
         pch = c(16, NA),
         lty = c(1, 2),
         lwd = 2,
         bty = "n")
}


# Step 5: Advanced Analysis and Visualization

# 1. Time-trend analysis for each treatment
# This examines how each treatment's response changes over time

# Function to perform time-trend analysis for a specific treatment
analyze_time_trend <- function(treatment_name) {
  # Subset data for this treatment
  treat_data <- data[data$Treatment == treatment_name, ]
  
  # Run linear model of CP_ratio over time
  time_model <- lm(CP_ratio ~ Day, data = treat_data)
  
  # Run model for acetone-adjusted values
  if (treatment_name != "Control1") {  # Control1 has no acetone adjustment
    adj_model <- lm(acetone_adjusted_CP_ratio ~ Day, data = treat_data)
  } else {
    adj_model <- time_model  # Same for Control1
  }
  
  # Return results
  return(list(
    treatment = treatment_name,
    time_model = time_model,
    adj_model = adj_model,
    orig_slope = coef(time_model)["Day"],
    orig_p = summary(time_model)$coefficients["Day", "Pr(>|t|)"],
    adj_slope = coef(adj_model)["Day"],
    adj_p = summary(adj_model)$coefficients["Day", "Pr(>|t|)"]
  ))
}

# Apply function to each treatment
treatments <- unique(data$Treatment)
time_trend_results <- lapply(treatments, analyze_time_trend)

# Compile results into a data frame
time_trend_summary <- data.frame(
  Treatment = sapply(time_trend_results, function(x) x$treatment),
  Original_Slope = sapply(time_trend_results, function(x) x$orig_slope),
  Original_P = sapply(time_trend_results, function(x) x$orig_p),
  Adjusted_Slope = sapply(time_trend_results, function(x) x$adj_slope),
  Adjusted_P = sapply(time_trend_results, function(x) x$adj_p)
)

# Add significance indicators
time_trend_summary$Original_Sig <- ifelse(time_trend_summary$Original_P < 0.05, "*", "")
time_trend_summary$Adjusted_Sig <- ifelse(time_trend_summary$Adjusted_P < 0.05, "*", "")

# Display time trend results
cat("\nTime Trend Analysis by Treatment:\n")
print(time_trend_summary)

# 2. Create visualization of time trends
par(mfrow=c(1,1))
par(mar=c(5, 5, 4, 8))  # Bottom, left, top, right margins

# Set up color palette
treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")

# Create empty plot
plot(NULL, xlim = range(data$Day), 
     ylim = range(c(time_trend_summary$Original_Slope, time_trend_summary$Adjusted_Slope)) * 1.2,
     xlab = "Treatment", ylab = "Time Slope (Change in C/P Ratio per Day)",
     main = "Time Trends by Treatment",
     xaxt = "n")

# Add grid
grid(nx = NA, ny = NULL, lty = 2, col = "gray90")

# Add custom x-axis
axis(1, at = 1:length(treatments), labels = treatments)

# Add horizontal line at y=0 (no change over time)
abline(h = 0, lty = 2)

# Plot original slopes
points(1:length(treatments), time_trend_summary$Original_Slope, 
       pch = 16, col = treatment_colors, cex = 1.5)

# Plot adjusted slopes
points(1:length(treatments) + 0.25, time_trend_summary$Adjusted_Slope, 
       pch = 17, col = treatment_colors, cex = 1.5)

# Connect pairs with lines
for (i in 1:length(treatments)) {
  lines(c(i, i + 0.25), 
        c(time_trend_summary$Original_Slope[i], time_trend_summary$Adjusted_Slope[i]),
        col = treatment_colors[i], lwd = 1.5)
  
  # Add significance indicators for original slopes
  if (time_trend_summary$Original_Sig[i] == "*") {
    text(i, time_trend_summary$Original_Slope[i], "*", pos = 3, cex = 1.2)
  }
  
  # Add significance indicators for adjusted slopes
  if (time_trend_summary$Adjusted_Sig[i] == "*") {
    text(i + 0.25, time_trend_summary$Adjusted_Slope[i], "*", pos = 3, cex = 1.2)
  }
}

# Add legend
legend("topright", inset = c(-0.25, 0), 
       legend = c("Original", "Acetone-adjusted"), 
       pch = c(16, 17),
       col = "black",
       cex = 0.9, 
       bty = "n",
       xpd = TRUE)

# 3. Analyze fungicide effect after controlling for time and acetone
# This directly models the effect of fungicide concentration

# Create log-transformed fungicide concentration
data$log_fungi <- ifelse(data$Fungi_Est > 0, log10(data$Fungi_Est), NA)

# Filter to treatments with fungicide
fungi_data <- data[!is.na(data$log_fungi), ]

# Run model
fungi_effect_model <- glm(acetone_adjusted_CP_ratio ~ Day + log_fungi, 
                          data = fungi_data, 
                          family = gaussian())

cat("\nFungicide Effect Model (controlling for time):\n")
print(summary(fungi_effect_model))

# 4. Create 3D visualization of the treatment-time-response relationship
# Set up data for 3D plotting with base R
treatments_numeric <- as.numeric(factor(data$Treatment, 
                                        levels = c("Control1", "Control2", 
                                                   "Treatment1", "Treatment2", 
                                                   "Treatment3", "Treatment4", 
                                                   "Treatment5")))

# Calculate the prediction surface
# Create grid of treatment and day values
treatment_grid <- rep(1:7, each = 50)
day_grid <- rep(seq(min(data$Day), max(data$Day), length.out = 50), 7)
grid_data <- data.frame(Treatment = factor(treatment_grid, 
                                           levels = 1:7,
                                           labels = c("Control1", "Control2", 
                                                      "Treatment1", "Treatment2", 
                                                      "Treatment3", "Treatment4", 
                                                      "Treatment5")),
                        Day = day_grid)

# Add all required columns that exist in the original data and are used by the model
# Add Acet_Est and Fungi_Est based on treatment
grid_data$Acet_Est <- c(
  "Control1" = 0,
  "Control2" = 100,
  "Treatment1" = 0.01,
  "Treatment2" = 0.1,
  "Treatment3" = 1, 
  "Treatment4" = 10,
  "Treatment5" = 100
)[grid_data$Treatment]

grid_data$Fungi_Est <- c(
  "Control1" = 0,
  "Control2" = 0,
  "Treatment1" = 0.01,
  "Treatment2" = 0.1,
  "Treatment3" = 1, 
  "Treatment4" = 10,
  "Treatment5" = 100
)[grid_data$Treatment]

# Add Group column if needed by the model
grid_data$Group <- as.numeric(factor(grid_data$Treatment, 
                                     levels = c("Control1", "Control2", 
                                                "Treatment1", "Treatment2", 
                                                "Treatment3", "Treatment4", 
                                                "Treatment5")))

# Add expected acetone effect for each treatment
grid_data$expected_acetone_effect <- (grid_data$Acet_Est / 100) * acetone_effect_estimate

# Now predict from ANCOVA model
grid_data$predicted <- predict(ancova_model, newdata = grid_data)

# Adjust for acetone effect
grid_data$adjusted_predicted <- grid_data$predicted - grid_data$expected_acetone_effect

# Create 3D perspective plot
# Convert to matrices for persp
x <- unique(day_grid)
y <- 1:7
z_matrix <- matrix(grid_data$adjusted_predicted, nrow = 50, ncol = 7)

# Create 3D plot
persp(x, y, z_matrix, 
      theta = 30, phi = 20, expand = 0.6, col = "lightblue",
      xlab = "Day", ylab = "Treatment", zlab = "Adjusted C/P Ratio",
      main = "3D Surface of Treatment-Time-Response Relationship",
      ticktype = "detailed")

# 5. Output final summary of findings
cat("\n\n=== SUMMARY OF FINDINGS ===\n")
cat("1. Acetone Effect Estimate:", acetone_effect_estimate, "\n")

# Summarize time trends
cat("2. Time Trends:\n")
for (i in 1:nrow(time_trend_summary)) {
  cat("   - ", time_trend_summary$Treatment[i], ": ", 
      round(time_trend_summary$Adjusted_Slope[i], 4), 
      " per day", 
      ifelse(time_trend_summary$Adjusted_P[i] < 0.05, " (significant)", ""), 
      "\n", sep = "")
}

# Summarize fungicide effect
fungi_effect <- coef(fungi_effect_model)["log_fungi"]
fungi_p <- summary(fungi_effect_model)$coefficients["log_fungi", "Pr(>|t|)"]
cat("3. Fungicide Effect (log scale):", round(fungi_effect, 4), 
    ifelse(fungi_p < 0.05, " (significant)", ""), "\n")

# Comparison of models
cat("4. Model Comparison:\n")
cat("   - AIC of traditional ANCOVA:", round(AIC(ancova_model), 2), "\n")
cat("   - AIC of acetone-adjusted model:", round(AIC(treatment_model), 2), "\n")
cat("   - Difference:", round(AIC(ancova_model) - AIC(treatment_model), 2), 
    " (", ifelse(AIC(ancova_model) > AIC(treatment_model), "better", "worse"), 
    " with acetone adjustment)\n", sep = "")

