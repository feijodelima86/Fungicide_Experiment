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