# Function to create plot with fungicide concentration on x-axis and days as different colors
create_concentration_plot <- function(data, y_var, se_var, y_label, main_title = "") {
  # Define colors for days - shades of blue from light to dark
  # Get unique days
  days <- sort(unique(data$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Get unique treatments and map to concentrations
  treatments <- levels(data$Treatment)
  
  # Map treatments to their fungicide concentrations (not acetone concentrations)
  fungicide_conc <- c(
    "Control1" = 0,
    "Control2" = 0,     # No fungicide, only acetone
    "Treatment1" = 0.01,
    "Treatment2" = 0.1,
    "Treatment3" = 1, 
    "Treatment4" = 10,
    "Treatment5" = 100
  )
  
  # Convert treatment to concentration values for the x-axis
  # We'll create a new data frame with concentrations
  plot_data <- data
  plot_data$Concentration <- fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine the y-axis limits including error bars
  y_min <- min(plot_data[[y_var]] - plot_data[[se_var]], na.rm = TRUE)
  y_max <- max(plot_data[[y_var]] + plot_data[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Use log scale for concentration (x-axis)
  # Add a small value to zero concentrations to avoid log(0)
  log_conc <- unique(plot_data$Concentration)
  log_conc[log_conc == 0] <- 0.001  # Small value to represent zero on log scale
  
  # Create empty plot with appropriate axes
  plot(NULL, 
       xlim = c(min(log_conc), max(fungicide_conc)), 
       ylim = c(y_min, y_max * 1.1),
       xlab = "Fungicide Concentration (%)", 
       ylab = y_label, 
       main = main_title,
       log = "x",  # Use log scale for x-axis
       xaxt = "n", # we'll add custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with specific concentrations
  conc_for_axis <- c(0.001, 0.01, 0.1, 1, 10, 100)
  axis_labels <- c("0", "0.01", "0.1", "1", "10", "100")
  axis(1, at = conc_for_axis, labels = axis_labels)
  
  # Add grid for readability
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray80")
  
  # Plot each day separately
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    
    # Sort by concentration to ensure lines connect points in order
    day_data <- day_data[order(day_data$Concentration), ]
    
    # Skip if no data for this day
    if(nrow(day_data) == 0) next
    
    # Replace zero with small value for log scale
    plot_conc <- day_data$Concentration
    plot_conc[plot_conc == 0] <- 0.001
    
    # Plot the line
    lines(plot_conc, day_data[[y_var]], 
          col = day_colors[i], lwd = 2)
    
    # Plot the points
    points(plot_conc, day_data[[y_var]], 
           col = day_colors[i], pch = day_pch[i], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(day_data)) {
      arrows(plot_conc[j], 
             day_data[[y_var]][j] - day_data[[se_var]][j],
             plot_conc[j], 
             day_data[[y_var]][j] + day_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = day_colors[i], lwd = 1.5)
    }
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.3, 0), 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Set up plotting layout - two plots stacked vertically
par(mfrow=c(2,1))

# Plot 1: Original CP Ratio
create_concentration_plot(mean_AI, "mean_AI", "se_AI", 
                          "AI")

# Plot 2: Acetone-Adjusted CP Ratio
create_concentration_plot(mean_adjusted, "mean_adjusted", "se_adjusted", 
                          "Adjusted AI")