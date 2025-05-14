library(readr)

# Analysis of C/P Ratio with Acetone Adjustment

# Read the data
data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

# Here you can search the entire document and replace the name of the variable of interest. 

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
summary(data$AI)

# Check the number of samples per treatment and day
print(table(data$Treatment, data$Day))

# Step 1: Estimate acetone effect from controls
# Run GLM with time as covariate comparing Control 1 vs Control 2
controls_data <- data[data$Treatment %in% c("Control1", "Control2"), ]
controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)

acetone_model <- glm(AI ~ Day + is_acetone, 
                     data = controls_data, 
                     family = gaussian())

plot(acetone_model)

print(summary(acetone_model))

# Extract the acetone effect (per 100% concentration)
acetone_effect_estimate <- coef(acetone_model)["is_acetone"]

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
data$acetone_adjusted_AI <- data$AI - data$expected_acetone_effect



# Step 3: Run analysis on acetone-adjusted responses
treatment_model <- glm(acetone_adjusted_AI ~ Day + Treatment + Acet_Est, 
                       data = data, 
                       family = gaussian())

print(summary(treatment_model))

# Calculate means and standard errors for original and adjusted CP ratios
mean_AI <- aggregate(AI ~ Treatment + Day, data = data, 
                            FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_AI <- do.call(data.frame, mean_AI)
names(mean_AI)[3:4] <- c("mean_AI", "se_AI")

mean_adjusted <- aggregate(acetone_adjusted_AI ~ Treatment + Day, data = data,
                           FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
mean_adjusted <- do.call(data.frame, mean_adjusted)
names(mean_adjusted)[3:4] <- c("mean_adjusted", "se_adjusted")

# View the means
print(head(mean_AI))
print(head(mean_adjusted))

# For plotting - ensure treatments are ordered correctly
mean_AI$Treatment <- factor(mean_AI$Treatment, 
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
create_treatment_plot <- function(data, y_var, se_var, y_label, mAIn_title) {
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
       xlab = "Day", ylab = y_label, mAIn = mAIn_title,
       xaxt = "n", # we'll add custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.mAIn = 1.3)
  
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
create_treatment_plot(mean_AI, "mean_AI", "se_AI", 
                      "AI", "")

# Plot 2: Acetone-Adjusted CP Ratio
create_treatment_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                      "Adjusted AI", "")

