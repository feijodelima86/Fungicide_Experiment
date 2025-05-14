
library(readr)

# Analysis of C/P Ratio with Acetone Adjustment

# Read the data
data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

# Here you can search the entire document and replace the name of the variable of interest. 

names(data)

data$AI<-log10(data$AI)

# Map group numbers to treatment names
treatment_map <- c(
  "1" = "Control1",    # Nothing Added
  "2" = "Control2",    # Acetone only, No fungicide
  "3" = "Treatment1",  # Lowest concentrations
  "4" = "Treatment2",  
  "5" = "Treatment3",  
  "6" = "Treatment4",  
  "7" = "Treatment5"   # Highest concentrations
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

# Map treatments to their actual fungicide concentrations (μg/L)
fungicide_conc <- c(
  "Control1" = 0,      # No fungicide
  "Control2" = 0,      # No fungicide, only acetone
  "Treatment1" = 0.054, # 0.054 μg/L (factor of 10 from 0.54)
  "Treatment2" = 0.54,  # 0.54 μg/L (factor of 10 from 5.4)
  "Treatment3" = 5.4,   # 5.4 μg/L (factor of 10 from 54)
  "Treatment4" = 54,    # 54 μg/L (factor of 10 from 540)
  "Treatment5" = 540    # 540 μg/L (highest concentration)
)

# Map treatments to their actual acetone concentrations (μg/L)
acetone_conc <- c(
  "Control1" = 0,      # No acetone
  "Control2" = 59,     # Highest acetone concentration, No fungicide
  "Treatment1" = 0.0059, # 0.0059 μg/L (factor of 10 from 0.059)
  "Treatment2" = 0.059,  # 0.059 μg/L (factor of 10 from 0.59)
  "Treatment3" = 0.59,   # 0.59 μg/L (factor of 10 from 5.9)
  "Treatment4" = 5.9,    # 5.9 μg/L (factor of 10 from 59)
  "Treatment5" = 59      # 59 μg/L (highest concentration)
)

# Step 1: Estimate acetone effect from controls
# Run GLM with time as covariate comparing Control 1 vs Control 2
controls_data <- data[data$Treatment %in% c("Control1", "Control2"), ]
controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)

acetone_model <- glm(AI ~ Day + is_acetone, 
                     data = controls_data, 
                     family = gaussian())

plot(acetone_model)

print(summary(acetone_model))

# Extract the acetone effect (per 59 μg/L concentration - the highest level)
acetone_effect_estimate <- coef(acetone_model)["is_acetone"]

# Step 2: Create proportional acetone effects
# Calculate expected acetone effect for each treatment based on concentration
data$expected_acetone_effect <- (acetone_conc[data$Treatment] / 59) * acetone_effect_estimate
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

# Plot 1: Original AI
create_treatment_plot(mean_AI, "mean_AI", "se_AI", 
                      "AI", "")

# Plot 2: Acetone-Adjusted AI
create_treatment_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                      "Adjusted AI", "")



# Function for the selective plotting from paste-2.txt
# Function to create the plot with error bars and selected treatments
create_treatment_plot_selective <- function(data, y_var, se_var, y_label, main_title, selected_treatments) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  treatments <- levels(droplevels(data_subset$Treatment))
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
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
  
  # Get treatment indices for color and point shape matching
  all_treatments <- levels(data$Treatment)
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data[data$Treatment == treatment, ]
    
    # Get the index of this treatment in the full treatment list for consistent colors
    idx <- treatment_indices[i]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[idx], lwd = 2)
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[idx], pch = treatment_pch[idx], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = treatment_colors[idx], lwd = 1.5)
    }
  }
  
  # Add legend outside the plot with only the selected treatments
  legend("topright", inset = c(-0.3, 0), 
         legend = selected_treatments, 
         col = treatment_colors[treatment_indices],
         pch = treatment_pch[treatment_indices],
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Set up plotting layout for selective plots
par(mfrow=c(2,1))

# Plot 1: Control 1, 2 and treatment 7 (Control1, Control2, Treatment5)
first_plot_treatments <- c("Control1", "Control2", "Treatment5")
create_treatment_plot_selective(mean_AI, "mean_AI", "se_AI", 
                                "AI", "Controls and Highest Treatment", 
                                first_plot_treatments)

# Plot 2: Treatments 1-5 (Treatment1 through Treatment5)
second_plot_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
create_treatment_plot_selective(mean_adjusted, "mean_adjusted", "se_adjusted", 
                                "Adjusted AI", "All Fungicide Treatments (Acetone-Adjusted)", 
                                second_plot_treatments)

# Concentration plot function with significance indicators
create_concentration_plot <- function(data, y_var, se_var, y_label, main_title = "", 
                                      significant_points = NULL) {
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- data[data$Treatment %in% treatments_to_include, ]
  
  # Get unique days
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[y_var]] - plot_data[[se_var]], na.rm = TRUE)
  y_max <- max(plot_data[[y_var]] + plot_data[[se_var]], na.rm = TRUE)
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(fungicide_conc[treatments_to_include])
  max_conc <- max(fungicide_conc[treatments_to_include])
  
  # Create empty plot
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(y_min, y_max * 1.2),  # Extra space for significance markers
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = y_label, 
       main = main_title,
       log = "x",  # Log scale for x-axis
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add x-axis with treatment concentrations
  conc_for_axis <- c(0.054, 0.54, 5.4, 54, 540)
  axis_labels <- c("0.054", "0.54", "5.4", "54", "540")
  axis(1, at = conc_for_axis, labels = axis_labels)
  
  # Add grid
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray80")
  
  # Plot each day separately
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    
    # Sort by concentration
    day_data <- day_data[order(day_data$Concentration), ]
    if(nrow(day_data) == 0) next
    
    # Plot the line
    lines(day_data$Concentration, day_data[[y_var]], 
          col = day_colors[i], lwd = 2)
    
    # Plot points and add significance markers if needed
    for (j in 1:nrow(day_data)) {
      # Check if this point is significant
      is_sig <- FALSE
      if (!is.null(significant_points)) {
        for (sig_point in significant_points) {
          if (sig_point$day == day && 
              sig_point$concentration == day_data$Concentration[j]) {
            is_sig <- TRUE
            break
          }
        }
      }
      
      # Draw the point
      points(day_data$Concentration[j], day_data[[y_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
      
      # Add significance marker if needed
      if (is_sig) {
        text(day_data$Concentration[j], 
             day_data[[y_var]][j] + 1.3*day_data[[se_var]][j], 
             "*", cex = 1.5, font = 2)  # Bold asterisk
      }
    }
    
    # Add error bars
    for (j in 1:nrow(day_data)) {
      arrows(day_data$Concentration[j], 
             day_data[[y_var]][j] - day_data[[se_var]][j],
             day_data$Concentration[j], 
             day_data[[y_var]][j] + day_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = day_colors[i], lwd = 1.5)
    }
  }
  
  # Add legend for days
  legend("topright", inset = c(-0.3, 0), 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
  
  # Add significance legend if needed
  if (!is.null(significant_points)) {
    legend("topleft", legend = "* p < 0.05", bty = "n", cex = 0.9)
  }
}

# Example of how to define significant points
# Replace with your actual significant findings
significant_points <- list(
  list(day = 5, concentration = 54),
  list(day = 10, concentration = 5.4),
  list(day = 15, concentration = 0.54)
)

# Set up for two stacked plots
par(mfrow = c(2,1))

# Plot 1: Original AI by fungicide concentration
create_concentration_plot(mean_AI, "mean_AI", "se_AI", 
                          "AI", "AI by Fungicide Concentration",
                          significant_points)

# Plot 2: Adjusted AI by fungicide concentration
create_concentration_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                          "Adjusted AI", "Adjusted AI by Fungicide Concentration",
                          significant_points)
