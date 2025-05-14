library(readr)

# Function for analyzing response variables to fungicide with acetone adjustment
analyze_fungicide_effect <- function(data, response_var, log_transform = FALSE, alpha = 0.05) {
  # Make a copy of the data to avoid modifying the original
  analysis_data <- data
  
  # Check for required column
  if (!"Group" %in% names(analysis_data) && "treatment_group" %in% names(analysis_data)) {
    analysis_data$Group <- analysis_data$treatment_group
  }
  
  # Log transform if requested
  if (log_transform) {
    analysis_data[[response_var]] <- log10(analysis_data[[response_var]])
    cat("Applied log10 transform to", response_var, "\n")
  }
  
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
  analysis_data$Treatment <- treatment_map[as.character(analysis_data$Group)]
  
  # Convert Day to numeric to ensure proper treatment in models
  analysis_data$Day <- as.numeric(as.character(analysis_data$Day))
  
  # Check the number of samples per treatment and day
  print(table(analysis_data$Treatment, analysis_data$Day))
  
  # Map treatments to their actual fungicide concentrations (μg/L)
  fungicide_conc <- c(
    "Control1" = 0,      # No fungicide
    "Control2" = 0,      # No fungicide, only acetone
    "Treatment1" = 0.054, # 0.054 μg/L
    "Treatment2" = 0.54,  # 0.54 μg/L
    "Treatment3" = 5.4,   # 5.4 μg/L
    "Treatment4" = 54,    # 54 μg/L
    "Treatment5" = 540    # 540 μg/L
  )
  
  # Map treatments to their actual acetone concentrations (μg/L)
  acetone_conc <- c(
    "Control1" = 0,      # No acetone
    "Control2" = 59,     # Highest acetone concentration, No fungicide
    "Treatment1" = 0.0059, # 0.0059 μg/L
    "Treatment2" = 0.059,  # 0.059 μg/L
    "Treatment3" = 0.59,   # 0.59 μg/L
    "Treatment4" = 5.9,    # 5.9 μg/L
    "Treatment5" = 59      # 59 μg/L
  )
  
  # Step 1: Estimate acetone effect from controls
  # Run GLM with time as covariate comparing Control 1 vs Control 2
  controls_data <- analysis_data[analysis_data$Treatment %in% c("Control1", "Control2"), ]
  controls_data$is_acetone <- ifelse(controls_data$Treatment == "Control2", 1, 0)
  
  if (log_transform) {
    family_type <- gaussian(link = "identity")
  } else {
    family_type <- gaussian()
  }
  
  acetone_model <- glm(as.formula(paste(response_var, "~ Day + is_acetone")), 
                       data = controls_data, 
                       family = family_type)
  
  print(summary(acetone_model))
  
  # Extract the acetone effect (per 59 μg/L concentration - the highest level)
  acetone_effect_estimate <- coef(acetone_model)["is_acetone"]
  
  # Step 2: Create proportional acetone effects
  # Calculate expected acetone effect for each treatment based on concentration
  analysis_data$expected_acetone_effect <- (acetone_conc[analysis_data$Treatment] / 59) * acetone_effect_estimate
  adj_var_name <- paste0("acetone_adjusted_", response_var)
  analysis_data[[adj_var_name]] <- analysis_data[[response_var]] - analysis_data$expected_acetone_effect
  
  # Step 3: Run analysis on acetone-adjusted responses
  # Check if 'Acet_Est' exists in the data
  acet_est_term <- if("Acet_Est" %in% names(analysis_data)) "+ Acet_Est" else ""
  
  treatment_model <- glm(as.formula(paste(adj_var_name, "~ Day + Treatment", acet_est_term)), 
                         data = analysis_data, 
                         family = family_type)
  
  print(summary(treatment_model))
  
  # Find significant differences between treatments and Control2 after acetone adjustment
  significant_points <- list()
  days <- sort(unique(analysis_data$Day))
  
  for (day in days) {
    day_data <- analysis_data[analysis_data$Day == day, ]
    
    # Compare each treatment to Control2 (acetone control) using the acetone-adjusted values
    treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
    for (trt in treatments_to_include) {
      trt_data <- day_data[day_data$Treatment %in% c("Control2", trt), ]
      if (nrow(trt_data) < 4) next  # Skip if too few observations
      
      trt_data$is_treatment <- ifelse(trt_data$Treatment == trt, 1, 0)
      model <- try(glm(as.formula(paste(adj_var_name, "~ is_treatment")), 
                       data = trt_data, family = family_type), silent = TRUE)
      
      if (!inherits(model, "try-error")) {
        p_val <- summary(model)$coefficients["is_treatment", "Pr(>|t|)"]
        if (!is.na(p_val) && p_val < alpha) {
          significant_points <- c(significant_points, 
                                  list(list(day = day, 
                                            treatment = trt,
                                            concentration = fungicide_conc[trt],
                                            p_value = p_val)))
        }
      }
    }
  }
  
  cat("\nSignificant differences from Control2 (acetone only) after acetone adjustment (p <", alpha, "):\n")
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      cat("Day:", point$day, "Treatment:", point$treatment, 
          "Concentration:", point$concentration, "p-value:", point$p_value, "\n")
    }
  } else {
    cat("No significant differences found.\n")
  }
  
  # Calculate means and standard errors for original and adjusted response variables
  mean_original <- aggregate(as.formula(paste(response_var, "~ Treatment + Day")), 
                             data = analysis_data, 
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                 se = sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))))
  mean_original <- do.call(data.frame, mean_original)
  colnames(mean_original)[3:4] <- c(paste0("mean_", response_var), paste0("se_", response_var))
  
  mean_adjusted <- aggregate(as.formula(paste(adj_var_name, "~ Treatment + Day")), 
                             data = analysis_data,
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                 se = sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))))
  mean_adjusted <- do.call(data.frame, mean_adjusted)
  colnames(mean_adjusted)[3:4] <- c(paste0("mean_adjusted_", response_var), 
                                    paste0("se_adjusted_", response_var))
  
  # For plotting - ensure treatments are ordered correctly
  mean_original$Treatment <- factor(mean_original$Treatment, 
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
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Generate all visualizations
  # Set up plotting layout - two plots stacked vertically
  par(mfrow=c(2,1))
  
  # VISUALIZATIONS: SECTION 1 - Original vs Adjusted plots
  # Original variable plot
  create_plot(mean_original, 
              paste0("mean_", response_var), 
              paste0("se_", response_var), 
              response_var, 
              "Original Data")
  
  # Acetone-Adjusted variable plot
  create_plot(mean_adjusted, 
              paste0("mean_adjusted_", response_var), 
              paste0("se_adjusted_", response_var), 
              paste0("Adjusted ", response_var), 
              "Acetone-Adjusted Data")
  
  # VISUALIZATIONS: SECTION 2 - Selective plots
  # Set up plotting layout for selective plots
  par(mfrow=c(2,1))
  
  # Plot 1: Control 1, 2 and treatment 5 (highest)
  first_plot_treatments <- c("Control1", "Control2", "Treatment5")
  create_selective_plot(mean_original, 
                        paste0("mean_", response_var), 
                        paste0("se_", response_var), 
                        response_var, 
                        "Controls and Highest Treatment", 
                        first_plot_treatments)
  
  # Plot 2: Treatments 1-5 (Treatment1 through Treatment5)
  second_plot_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  create_selective_plot(mean_adjusted, 
                        paste0("mean_adjusted_", response_var), 
                        paste0("se_adjusted_", response_var), 
                        paste0("Adjusted ", response_var), 
                        "All Fungicide Treatments (Acetone-Adjusted)", 
                        second_plot_treatments)
  
  # VISUALIZATIONS: SECTION 3 - Concentration plots
  # Set up for concentration plots
  par(mfrow = c(2,1))
  
  # Plot 1: Original variable by fungicide concentration
  create_concentration_plot_adaptive(mean_original, 
                            paste0("mean_", response_var), 
                            paste0("se_", response_var), 
                            response_var, 
                            "Original Data by Fungicide Concentration",
                            fungicide_conc,
                            NULL)  # No significance testing for original data
  
  # Plot 2: Adjusted variable by fungicide concentration with significance
  create_concentration_plot_adaptive(mean_adjusted, 
                            paste0("mean_adjusted_", response_var), 
                            paste0("se_adjusted_", response_var), 
                            paste0("Adjusted ", response_var), 
                            "Acetone-Adjusted Data by Fungicide Concentration\nAsterisks (*) indicate significant differences from Control2",
                            fungicide_conc,
                            significant_points)
  
  # Return list of results
  invisible(list(
    model_acetone = acetone_model,
    model_treatment = treatment_model,
    data_adjusted = analysis_data,
    means_original = mean_original,
    means_adjusted = mean_adjusted,
    acetone_effect = acetone_effect_estimate,
    significant_points = significant_points
  ))
}

# Function to create the treatment plot with error bars
create_plot <- function(data, y_var, se_var, y_label, main_title) {
  # Get the unique treatments and days
  treatments <- levels(data$Treatment)
  days <- sort(unique(data$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data[[y_var]] - data[[se_var]], na.rm = TRUE)
  y_max <- max(data[[y_var]] + data[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
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

# Function to create the selective plot with error bars
create_selective_plot <- function(data, y_var, se_var, y_label, main_title, selected_treatments) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
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
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
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

# Modified concentration plot function with adaptive y-axis scaling
create_concentration_plot_adaptive <- function(data, y_var, se_var, y_label, main_title, 
                                               fungicide_conc, significant_points = NULL, log_transformed = FALSE) {
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
  
  # For log-transformed data, adjust y-axis padding
  if (log_transformed) {
    # For log-transformed data, use less padding
    y_padding <- 0.3  # About 30% padding for log-scale
    y_min_adjusted <- y_min - (y_padding * abs(y_min))
    y_max_adjusted <- y_max + (y_padding * abs(y_max))
  } else {
    # For untransformed data, need more padding especially for zero-bounded data
    y_range <- y_max - y_min
    y_min_adjusted <- max(0, y_min - (0.1 * y_range))  # Don't go below 0 for natural AI values
    y_max_adjusted <- y_max + (0.2 * y_range)  # Add 20% padding at the top
  }
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(fungicide_conc[treatments_to_include])
  max_conc <- max(fungicide_conc[treatments_to_include])
  
  # Create empty plot
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(y_min_adjusted, y_max_adjusted),
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
      if (!is.null(significant_points) && length(significant_points) > 0) {
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
        # For log-transformed data, position asterisks with less padding
        if (log_transformed) {
          text(day_data$Concentration[j], 
               day_data[[y_var]][j] + 1.1*day_data[[se_var]][j], 
               "*", cex = 1.5, font = 2)  # Bold asterisk
        } else {
          text(day_data$Concentration[j], 
               day_data[[y_var]][j] + 1.3*day_data[[se_var]][j], 
               "*", cex = 1.5, font = 2)  # Bold asterisk
        }
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
  if (!is.null(significant_points) && length(significant_points) > 0) {
    legend("topleft", 
           legend = "* p < 0.05 (different from Control2)", 
           bty = "n", cex = 0.9)
  }
}

# Example usage:
# data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")
 

 # # Analyze AI (Autotrophic Index) with log transform
# analyze_fungicide_effect(data, "AI", log_transform = TRUE)

 
# # Analyze C.P.ratio without log transform
# analyze_fungicide_effect(data, "C.P.ratio", log_transform = FALSE)
# 
# # Analyze another variable with custom alpha
# analyze_fungicide_effect(data, "Chlorophyll.a", log_transform = FALSE, alpha = 0.1)