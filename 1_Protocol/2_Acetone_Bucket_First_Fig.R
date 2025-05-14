library(readr)
# Fix Matrix package issue first
if(requireNamespace("Matrix", quietly = TRUE)) {
  # If Matrix is installed but outdated, try to update it
  update.packages("Matrix", repos = "https://cloud.r-project.org")
} else {
  # If Matrix isn't installed, install it
  install.packages("Matrix", repos = "https://cloud.r-project.org")
}

# Now install/load lme4 and dependencies
if(!requireNamespace("lme4", quietly = TRUE)) {
  install.packages("lme4", repos = "https://cloud.r-project.org", dependencies = TRUE)
}
library(lme4)  # For mixed-effects models

# Optional: adds p-values to lmer models
if(!requireNamespace("lmerTest", quietly = TRUE)) {
  install.packages("lmerTest", repos = "https://cloud.r-project.org")
}
library(lmerTest)

# Analysis of C/P Ratio with Acetone Adjustment

# Read the data
data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

names(data)

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
str(data)
summary(data$SiO2)

# Make sure Sample is a factor for use as random effect
if(!"Sample" %in% names(data)) {
  # If Sample column doesn't exist, create a unique identifier
  data$Sample <- factor(paste(data$Treatment, seq_along(data$Treatment), sep="_"))
} else {
  # Ensure Sample is a factor
  data$Sample <- factor(data$Sample)
}

# Check the number of samples per treatment and day
print(table(data$Treatment, data$Day))

# Step 1: Estimate acetone effect from controls
controls_data <- data[data$Treatment %in% c("Control1", "Control2"), ]
controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)

# Try mixed model first, but fall back to simpler model if needed
tryCatch({
  # Attempt to use linear mixed-effects model with random intercept for Sample
  acetone_model2 <- lmer(SiO2 ~ Day + is_acetone + (1|Sample), 
                        data = controls_data)
  
  # Extract the acetone effect from the mixed model
  acetone_effect_estimate <- fixef(acetone_model)["is_acetone"]
  
  # Check model
  cat("Successfully fit mixed effects model\n")
  print(summary(acetone_model))
  
}, error = function(e) {
  # If the mixed model fails, fall back to a regular linear model
  message("Mixed model failed with error: ", e$message)
  message("Falling back to standard linear model")
  
  acetone_model <- lm(SiO2 ~ Day + is_acetone, data = controls_data)
  
  # Extract the acetone effect from the linear model
  acetone_effect_estimate <<- coef(acetone_model)["is_acetone"]
  
  # Check model
  print(summary(acetone_model))
})

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
data$acetone_adjusted_SiO2 <- data$SiO2 - data$expected_acetone_effect

# Step 3: Run analysis on acetone-adjusted responses
# Try mixed model first, but fall back to simpler model if needed
tryCatch({
  # Include the Acet_Est variable if it exists in your dataset
  if("Acet_Est" %in% names(data)) {
    treatment_model <- lmer(acetone_adjusted_SiO2 ~ Day + Treatment + Acet_Est + (1|Sample), 
                            data = data)
  } else {
    treatment_model <- lmer(acetone_adjusted_SiO2 ~ Day + Treatment + (1|Sample), 
                            data = data)
  }
  
  cat("Successfully fit mixed effects model for treatment analysis\n")
  print(summary(treatment_model))
  
}, error = function(e) {
  # If the mixed model fails, fall back to a regular linear model
  message("Mixed model failed with error: ", e$message)
  message("Falling back to standard linear model for treatment analysis")
  
  if("Acet_Est" %in% names(data)) {
    treatment_model <<- lm(acetone_adjusted_SiO2 ~ Day + Treatment + Acet_Est, 
                           data = data)
  } else {
    treatment_model <<- lm(acetone_adjusted_SiO2 ~ Day + Treatment, 
                           data = data)
  }
  
  print(summary(treatment_model))
})



# Calculate means and standard errors for original and adjusted SiO2
mean_SiO2 <- aggregate(SiO2 ~ Treatment + Day, data = data, 
                       FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_SiO2 <- do.call(data.frame, mean_SiO2)
names(mean_SiO2)[3:4] <- c("mean_SiO2", "se_SiO2")

mean_adjusted <- aggregate(acetone_adjusted_SiO2 ~ Treatment + Day, data = data,
                           FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_adjusted <- do.call(data.frame, mean_adjusted)
names(mean_adjusted)[3:4] <- c("mean_adjusted", "se_adjusted")

# View the means
print(head(mean_SiO2))
print(head(mean_adjusted))

# For plotting - ensure treatments are ordered correctly
mean_SiO2$Treatment <- factor(mean_SiO2$Treatment, 
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

# Plot 1: Original SiO2 values
create_treatment_plot(mean_SiO2, "mean_SiO2", "se_SiO2", 
                      "SiO2", "Original SiO2 Values")

# Plot 2: Acetone-Adjusted SiO2 values
create_treatment_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                      "Adjusted SiO2", "Acetone-Adjusted SiO2 Values")

