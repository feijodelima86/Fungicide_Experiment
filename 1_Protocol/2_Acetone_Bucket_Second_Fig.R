# Code to generate separate plots for specific treatment groups

# Define a color palette for the treatments
treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
treatment_pch <- c(15, 16, 21, 17, 18, 19, 20)  # Different point symbols

# Function to create the plot with error bars
create_treatment_plot <- function(data, treatments_to_include, y_var, se_var, y_label, main_title) {
  # Filter the data for the specified treatments
  plot_data <- data[data$Treatment %in% treatments_to_include, ]
  
  # Get the unique treatments and days
  treatments <- treatments_to_include[treatments_to_include %in% levels(data$Treatment)]
  days <- sort(unique(plot_data$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(plot_data[[y_var]] - plot_data[[se_var]], na.rm = TRUE)
  y_max <- max(plot_data[[y_var]] + plot_data[[se_var]], na.rm = TRUE)
  
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
    treatment_data <- plot_data[plot_data$Treatment == treatment, ]
    
    # Get the index of this treatment in the original treatment list (for consistent colors)
    original_index <- which(levels(data$Treatment) == treatment)
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[original_index], lwd = 2)
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[original_index], pch = treatment_pch[original_index], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = treatment_colors[original_index], lwd = 1.5)
    }
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.3, 0), 
         legend = treatments, 
         col = treatment_colors[match(treatments, levels(data$Treatment))],
         pch = treatment_pch[match(treatments, levels(data$Treatment))],
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Set up plotting parameters
# Create a new graphics device for these plots
# Set up plotting layout - two plots stacked vertically
par(mfrow=c(2,1))

# Define the treatment groups
controls_and_t5 <- c("Control1", "Control2", "Treatment5")
treatments_1_to_5 <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")

# PLOT 1: Control 1, Control 2, and Treatment 5 (UNCORRECTED)
create_treatment_plot(
  mean_SiO2, 
  controls_and_t5,
  "mean_SiO2", 
  "se_SiO2", 
  "SiO2", 
  "Control Groups and Treatment 5 (Uncorrected)"
)

# PLOT 2: Treatments 1-5 (CORRECTED for acetone)
create_treatment_plot(
  mean_adjusted, 
  treatments_1_to_5,
  "mean_adjusted", 
  "se_adjusted", 
  "Acetone-Adjusted SiO2", 
  "Treatments 1-5 (Acetone-Adjusted)"
)