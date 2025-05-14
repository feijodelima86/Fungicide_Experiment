# Modified analyze_fungicide_effect function to test for differences across treatments at specific days
analyze_fungicide_effect <- function(data, response_var, log_transform = FALSE, alpha = 0.05) {
  # Make a copy of the data to avoid modifying the original
  analysis_data <- data
  
  # Check for required column
  if (!"Group" %in% names(analysis_data) && "treatment_group" %in% names(analysis_data)) {
    analysis_data$Group <- analysis_data$treatment_group
  }
  
  # Original response variable name for labeling
  original_var_name <- response_var
  
  # Check for negative values if log transform is requested
  if (log_transform) {
    min_value <- min(analysis_data[[response_var]], na.rm = TRUE)
    if (min_value <= 0) {
      warning(paste("Response variable contains values ≤ 0. Cannot apply log10 transform directly.",
                    "Consider adding a constant or using a different transformation."))
      
      # Ask if user wants to add a constant for log transformation
      cat("Minimum value is", min_value, "\n")
      cat("Would you like to add a constant to make all values positive for log transformation? (Y/N): ")
      user_choice <- readline()
      
      if (tolower(user_choice) == "y") {
        # Add a constant to make all values positive
        offset <- abs(min_value) + 1  # Add 1 to ensure all values are > 0
        analysis_data[[response_var]] <- analysis_data[[response_var]] + offset
        cat("Added constant of", offset, "to all values before log transform\n")
        analysis_data[[response_var]] <- log10(analysis_data[[response_var]])
        cat("Applied log10 transform to", response_var, "\n")
      } else {
        # Skip log transformation
        cat("Skipping log transformation as requested\n")
        log_transform <- FALSE
      }
    } else {
      # Apply log transform normally
      analysis_data[[response_var]] <- log10(analysis_data[[response_var]])
      cat("Applied log10 transform to", response_var, "\n")
    }
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
  analysis_data$expected_acetone_effect <- (acetone_conc[analysis_data$Treatment] / 59) * acetone_effect_estimate
  adj_var_name <- paste0("acetone_adjusted_", response_var)
  analysis_data[[adj_var_name]] <- analysis_data[[response_var]] - analysis_data$expected_acetone_effect
  
  # Step 3: Run analysis on acetone-adjusted responses
  acet_est_term <- if("Acet_Est" %in% names(analysis_data)) "+ Acet_Est" else ""
  
  treatment_model <- glm(as.formula(paste(adj_var_name, "~ Day + Treatment", acet_est_term)), 
                         data = analysis_data, 
                         family = family_type)
  
  print(summary(treatment_model))
  
  # NEW: Test for differences across all treatments at each specific day
  significant_points <- list()
  days <- sort(unique(analysis_data$Day))
  
  for (day in days) {
    cat("\n=== Day", day, "Treatment Comparisons ===\n")
    
    # Filter data for the current day
    day_data <- analysis_data[analysis_data$Day == day, ]
    
    # Only perform ANOVA if we have enough treatments with data
    if (length(unique(day_data$Treatment)) >= 2) {
      # Run ANOVA for all treatments on this day
      anova_formula <- as.formula(paste(adj_var_name, "~ Treatment"))
      anova_model <- try(aov(anova_formula, data = day_data), silent = TRUE)
      
      if (!inherits(anova_model, "try-error")) {
        # Print the ANOVA summary
        anova_result <- summary(anova_model)
        print(anova_result)
        
        # Check if there's a significant treatment effect
        p_value <- anova_result[[1]]["Treatment", "Pr(>F)"]
        if (!is.na(p_value) && p_value < alpha) {
          cat("Significant treatment effect on Day", day, "(p =", p_value, ")\n")
          
          # Perform Tukey's HSD post-hoc test to find specific differences
          tukey_result <- TukeyHSD(anova_model, "Treatment", conf.level = 1-alpha)
          print(tukey_result)
          
          # Extract significant comparisons
          tukey_p_values <- tukey_result$Treatment[, "p adj"]
          significant_pairs <- which(tukey_p_values < alpha)
          
          if (length(significant_pairs) > 0) {
            for (i in significant_pairs) {
              pair_name <- names(tukey_p_values)[i]
              treatments <- strsplit(pair_name, "-")[[1]]
              difference <- tukey_result$Treatment[pair_name, "diff"]
              p_val <- tukey_p_values[i]
              
              # Add to significant points
              significant_points <- c(significant_points, 
                                      list(list(day = day,
                                                treatment1 = treatments[1],
                                                treatment2 = treatments[2],
                                                difference = difference,
                                                p_value = p_val)))
            }
          }
        } else {
          cat("No significant treatment effect on Day", day, "\n")
        }
      } else {
        cat("Error in ANOVA for Day", day, "\n")
      }
    } else {
      cat("Not enough treatment groups on Day", day, "for comparison\n")
    }
  }
  
  # Print summary of significant differences
  cat("\n=== Summary of Significant Treatment Differences (p <", alpha, ") ===\n")
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      sign_char <- ifelse(point$difference > 0, ">", "<")
      cat("Day:", point$day, ":", point$treatment1, sign_char, point$treatment2, 
          "Difference:", round(point$difference, 4), 
          "p-value:", round(point$p_value, 4), "\n")
    }
  } else {
    cat("No significant differences found between treatments.\n")
  }
  
  # Create safer variable name suffixes for aggregation
  mean_var <- paste0("mean_", response_var)
  se_var <- paste0("se_", response_var)
  mean_adj_var <- paste0("mean_adjusted_", response_var)
  se_adj_var <- paste0("se_adjusted_", response_var)
  
  # Calculate means and standard errors for original and adjusted response variables
  mean_original <- aggregate(as.formula(paste(response_var, "~ Treatment + Day")), 
                             data = analysis_data, 
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                 se = sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))))
  mean_original <- do.call(data.frame, mean_original)
  colnames(mean_original)[3:4] <- c(mean_var, se_var)
  
  mean_adjusted <- aggregate(as.formula(paste(adj_var_name, "~ Treatment + Day")), 
                             data = analysis_data,
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                                 se = sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))))
  mean_adjusted <- do.call(data.frame, mean_adjusted)
  colnames(mean_adjusted)[3:4] <- c(mean_adj_var, se_adj_var)
  
  # Ensure treatments are ordered correctly
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
  
  # Return list of results
  result <- list(
    model_acetone = acetone_model,
    model_treatment = treatment_model,
    data_adjusted = analysis_data,
    means_original = mean_original,
    means_adjusted = mean_adjusted,
    acetone_effect = acetone_effect_estimate,
    significant_points = significant_points,
    # Store column names for easier reference
    mean_var = mean_var,
    se_var = se_var,
    mean_adj_var = mean_adj_var,
    se_adj_var = se_adj_var,
    response_var = response_var,
    original_var_name = original_var_name,
    log_transform = log_transform,
    fungicide_conc = fungicide_conc
  )
  
  return(result)
}
# Function to create the treatment plot with error bars
create_plot <- function(data, y_var, se_var, y_label, main_title, log_scale = FALSE) {
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
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max * 1.1),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
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
         legend = c("Controls:", treatments[1:2], "Treatments:", treatments[3:7]), 
         col = c(NA, treatment_colors[1:2], NA, treatment_colors[3:7]),
         pch = c(NA, treatment_pch[1:2], NA, treatment_pch[3:7]),
         lty = c(NA, 1, 1, NA, 1, 1, 1, 1, 1),
         lwd = c(NA, 2, 2, NA, 2, 2, 2, 2, 2),
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}
# Function to create the selective plot with error bars
create_selective_plot <- function(data, y_var, se_var, y_label, main_title, selected_treatments, log_scale = FALSE) {
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
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max * 1.1),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
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
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}

# Function to create selective plot with significance indicators for specific days
create_selective_plot_sig_days <- function(data, y_var, se_var, y_label, main_title, 
                                           selected_treatments, significant_points = NULL,
                                           selected_day = NULL, log_scale = FALSE) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Filter for specific day if provided
  if (!is.null(selected_day)) {
    data_subset <- data_subset[data_subset$Day == selected_day, ]
    if (nrow(data_subset) == 0) {
      stop(paste("No data available for day", selected_day))
    }
    days <- selected_day
  } else {
    # Get unique days
    days <- sort(unique(data_subset$Day))
  }
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max * 1.2),  # Extra space for significance markers
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_subset[data_subset$Treatment == treatment, ]
    
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
  
  # Add significance markers if provided
  if (!is.null(significant_points) && length(significant_points) > 0) {
    # Get the significant comparisons involving the selected treatments
    filtered_sig_points <- Filter(function(x) {
      x$treatment1 %in% selected_treatments && 
        x$treatment2 %in% selected_treatments &&
        (is.null(selected_day) || x$day == selected_day)
    }, significant_points)
    
    if (length(filtered_sig_points) > 0) {
      # Group significance markers by day
      days_with_sig <- unique(sapply(filtered_sig_points, function(x) x$day))
      
      for (day in days_with_sig) {
        # Get comparisons for this day
        day_comparisons <- Filter(function(x) x$day == day, filtered_sig_points)
        
        if (length(day_comparisons) > 0) {
          # Calculate y positions for significance markers
          y_start <- y_max * 1.05
          y_spacing <- (y_max * 0.1) / (length(day_comparisons) + 1)
          
          for (k in 1:length(day_comparisons)) {
            comp <- day_comparisons[[k]]
            
            # Add asterisks based on p-value
            p_val <- comp$p_value
            asterisks <- ifelse(p_val < 0.001, "***", 
                                ifelse(p_val < 0.01, "**", 
                                       ifelse(p_val < 0.05, "*", "")))
            
            # If we have two specific data points, draw a line between them
            y_pos <- y_start + (k * y_spacing)
            text(day, y_pos, asterisks, cex = 1.2)
            
            # Add small ticks to indicate which treatments are being compared
            trt1_idx <- match(comp$treatment1, selected_treatments)
            trt2_idx <- match(comp$treatment2, selected_treatments)
            
            if (!is.na(trt1_idx) && !is.na(trt2_idx)) {
              # Find the data points for these treatments on this day
              trt1_data <- data_subset[data_subset$Treatment == comp$treatment1 & data_subset$Day == day, ]
              trt2_data <- data_subset[data_subset$Treatment == comp$treatment2 & data_subset$Day == day, ]
              
              if (nrow(trt1_data) > 0 && nrow(trt2_data) > 0) {
                # Draw vertical lines from data points to significance line
                segments(day, trt1_data[[y_var]], day, y_pos, 
                         col = treatment_colors[treatment_indices[trt1_idx]], 
                         lty = 2, lwd = 1)
                segments(day, trt2_data[[y_var]], day, y_pos, 
                         col = treatment_colors[treatment_indices[trt2_idx]], 
                         lty = 2, lwd = 1)
              }
            }
          }
        }
      }
      
      # Add legend for significance symbols
      legend("topleft", 
             legend = c("* p < 0.05", "** p < 0.01", "*** p < 0.001"), 
             bty = "n", cex = 0.8)
    }
  }
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}
# Function to create selective plot with significance groups (A, B, C notation)
create_selective_plot_sig_days <- function(data, y_var, se_var, y_label, main_title, 
                                           selected_treatments, significant_points = NULL,
                                           selected_day = NULL, log_scale = FALSE) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Filter for specific day if provided
  if (!is.null(selected_day)) {
    data_subset <- data_subset[data_subset$Day == selected_day, ]
    if (nrow(data_subset) == 0) {
      stop(paste("No data available for day", selected_day))
    }
    days <- selected_day
  } else {
    # Get unique days
    days <- sort(unique(data_subset$Day))
  }
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max * 1.1),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
  # Calculate the statistical groupings (A, B, C, etc.) for each day
  group_labels <- list()
  
  if (!is.null(significant_points) && length(significant_points) > 0) {
    # Process each day
    for (day in days) {
      # Get significant comparisons for this day among selected treatments
      day_comparisons <- Filter(function(x) {
        x$day == day && 
          x$treatment1 %in% selected_treatments && 
          x$treatment2 %in% selected_treatments
      }, significant_points)
      
      if (length(day_comparisons) > 0) {
        # Initialize all treatments with group "A"
        day_groups <- setNames(rep("A", length(selected_treatments)), selected_treatments)
        
        # Add additional group letters based on significant differences
        current_letter <- "A"
        
        # Sort treatments by their mean values (descending) for this day
        treatment_means <- sapply(selected_treatments, function(trt) {
          trt_data <- data_subset[data_subset$Treatment == trt & data_subset$Day == day, ]
          if (nrow(trt_data) > 0) trt_data[[y_var]] else NA
        })
        
        sorted_treatments <- names(sort(treatment_means, decreasing = TRUE))
        
        # First pass: assign initial groups based on highest value
        for (i in 2:length(sorted_treatments)) {
          trt_i <- sorted_treatments[i]
          
          # Check if this treatment differs from any treatment with higher mean
          for (j in 1:(i-1)) {
            trt_j <- sorted_treatments[j]
            
            # Look for this comparison in the significant differences
            is_different <- any(sapply(day_comparisons, function(comp) {
              (comp$treatment1 == trt_i && comp$treatment2 == trt_j) || 
                (comp$treatment1 == trt_j && comp$treatment2 == trt_i)
            }))
            
            if (is_different) {
              # This treatment is in a different group
              next_letter <- LETTERS[which(LETTERS == current_letter) + 1]
              day_groups[trt_i] <- next_letter
              current_letter <- next_letter
              break  # Only increase letter once per treatment
            }
          }
        }
        
        # Second pass: refine groups to make them transitive (A, AB, B, etc.)
        for (i in 1:length(sorted_treatments)) {
          trt_i <- sorted_treatments[i]
          
          for (j in 1:length(sorted_treatments)) {
            if (i == j) next
            trt_j <- sorted_treatments[j]
            
            # Look for this comparison in the significant differences
            is_different <- any(sapply(day_comparisons, function(comp) {
              (comp$treatment1 == trt_i && comp$treatment2 == trt_j) || 
                (comp$treatment1 == trt_j && comp$treatment2 == trt_i)
            }))
            
            if (!is_different && day_groups[trt_i] != day_groups[trt_j]) {
              # These treatments are not significantly different but have different groups
              # Create an intermediate group (like AB)
              group_i <- day_groups[trt_i]
              group_j <- day_groups[trt_j]
              
              # Combine groups if they're not already combined
              if (nchar(group_i) == 1 && nchar(group_j) == 1) {
                combined_group <- paste0(min(group_i, group_j), max(group_i, group_j))
                if (group_i < group_j) {
                  day_groups[trt_j] <- combined_group
                } else {
                  day_groups[trt_i] <- combined_group
                }
              }
            }
          }
        }
        
        # Store the group labels for this day
        group_labels[[as.character(day)]] <- day_groups
      }
    }
  }
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_subset[data_subset$Treatment == treatment, ]
    
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
      
      # Add group label if available for this day
      day <- treatment_data$Day[j]
      if (as.character(day) %in% names(group_labels)) {
        day_group <- group_labels[[as.character(day)]]
        if (treatment %in% names(day_group)) {
          text(day + 0.15, treatment_data[[y_var]][j], 
               day_group[treatment], 
               col = treatment_colors[idx], 
               cex = 0.8, pos = 4)  # pos=4 means position text to the right
        }
      }
    }
  }
  
  # Add legend for explaining group labels
  if (length(group_labels) > 0) {
    legend("topleft", 
           legend = c("Treatments labeled with the same letter", "are not significantly different (p ≥ 0.05)"), 
           bty = "n", cex = 0.8)
  }
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}
# Function to create selective plot with black bold letters for significance groups (A, B, C notation)
create_selective_plot_sig_days <- function(data, y_var, se_var, y_label, main_title, 
                                           selected_treatments, significant_points = NULL,
                                           selected_day = NULL, log_scale = FALSE) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Filter for specific day if provided
  if (!is.null(selected_day)) {
    data_subset <- data_subset[data_subset$Day == selected_day, ]
    if (nrow(data_subset) == 0) {
      stop(paste("No data available for day", selected_day))
    }
    days <- selected_day
  } else {
    # Get unique days
    days <- sort(unique(data_subset$Day))
  }
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max * 1.1),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
  # Calculate the statistical groupings (A, B, C, etc.) for each day separately
  # This ensures differences across treatments are only calculated within each day
  group_labels <- list()
  
  if (!is.null(significant_points) && length(significant_points) > 0) {
    # Process each day individually
    for (day in days) {
      # Get significant comparisons ONLY for this specific day among selected treatments
      day_comparisons <- Filter(function(x) {
        x$day == day && 
          x$treatment1 %in% selected_treatments && 
          x$treatment2 %in% selected_treatments
      }, significant_points)
      
      if (length(day_comparisons) > 0) {
        # Initialize all treatments with group "A"
        day_groups <- setNames(rep("A", length(selected_treatments)), selected_treatments)
        
        # Add additional group letters based on significant differences
        current_letter <- "A"
        
        # Sort treatments by their mean values (descending) for this day
        treatment_means <- sapply(selected_treatments, function(trt) {
          trt_data <- data_subset[data_subset$Treatment == trt & data_subset$Day == day, ]
          if (nrow(trt_data) > 0) trt_data[[y_var]] else NA
        })
        
        sorted_treatments <- names(sort(treatment_means, decreasing = TRUE))
        
        # First pass: assign initial groups based on highest value
        for (i in 2:length(sorted_treatments)) {
          trt_i <- sorted_treatments[i]
          
          # Check if this treatment differs from any treatment with higher mean
          for (j in 1:(i-1)) {
            trt_j <- sorted_treatments[j]
            
            # Look for this comparison in the significant differences
            is_different <- any(sapply(day_comparisons, function(comp) {
              (comp$treatment1 == trt_i && comp$treatment2 == trt_j) || 
                (comp$treatment1 == trt_j && comp$treatment2 == trt_i)
            }))
            
            if (is_different) {
              # This treatment is in a different group
              next_letter <- LETTERS[which(LETTERS == current_letter) + 1]
              day_groups[trt_i] <- next_letter
              current_letter <- next_letter
              break  # Only increase letter once per treatment
            }
          }
        }
        
        # Second pass: refine groups to make them transitive (A, AB, B, etc.)
        for (i in 1:length(sorted_treatments)) {
          trt_i <- sorted_treatments[i]
          
          for (j in 1:length(sorted_treatments)) {
            if (i == j) next
            trt_j <- sorted_treatments[j]
            
            # Look for this comparison in the significant differences
            is_different <- any(sapply(day_comparisons, function(comp) {
              (comp$treatment1 == trt_i && comp$treatment2 == trt_j) || 
                (comp$treatment1 == trt_j && comp$treatment2 == trt_i)
            }))
            
            if (!is_different && day_groups[trt_i] != day_groups[trt_j]) {
              # These treatments are not significantly different but have different groups
              # Create an intermediate group (like AB)
              group_i <- day_groups[trt_i]
              group_j <- day_groups[trt_j]
              
              # Combine groups if they're not already combined
              if (nchar(group_i) == 1 && nchar(group_j) == 1) {
                combined_group <- paste0(min(group_i, group_j), max(group_i, group_j))
                if (group_i < group_j) {
                  day_groups[trt_j] <- combined_group
                } else {
                  day_groups[trt_i] <- combined_group
                }
              }
            }
          }
        }
        
        # Store the group labels for this day
        group_labels[[as.character(day)]] <- day_groups
      }
    }
  }
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_subset[data_subset$Treatment == treatment, ]
    
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
      
      # Add group label if available for this day
      day <- treatment_data$Day[j]
      if (as.character(day) %in% names(group_labels)) {
        day_group <- group_labels[[as.character(day)]]
        if (treatment %in% names(day_group)) {
          # Use black and bold text for group labels
          text(day + 0.15, treatment_data[[y_var]][j], 
               day_group[treatment], 
               col = "black",  # Always black letters
               font = 2,       # Bold font
               cex = 0.9, pos = 4)  # pos=4 means position text to the right
        }
      }
    }
  }
  
  # Add legend for explaining group labels
  if (length(group_labels) > 0) {
    legend("topleft", 
           legend = c("Treatments labeled with the same letter", "within a day are not significantly different (p ≥ 0.05)"), 
           bty = "n", cex = 0.8)
  }
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}
# Function to create selective plot with black bold letters for significance groups (A, B, C notation)
create_selective_plot_sig_days <- function(data, y_var, se_var, y_label, main_title, 
                                           selected_treatments, significant_points = NULL,
                                           selected_day = NULL, log_scale = FALSE) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Filter for specific day if provided
  if (!is.null(selected_day)) {
    data_subset <- data_subset[data_subset$Day == selected_day, ]
    if (nrow(data_subset) == 0) {
      stop(paste("No data available for day", selected_day))
    }
    days <- selected_day
  } else {
    # Get unique days
    days <- sort(unique(data_subset$Day))
  }
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Adjusted y-max to be 1.2 * (max value + max standard error)
  max_value_with_se <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  y_max <- 1.2 * max_value_with_se
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
  # Calculate the statistical groupings (A, B, C, etc.) for each day separately
  # This ensures differences across treatments are only calculated within each day
  group_labels <- list()
  
  if (!is.null(significant_points) && length(significant_points) > 0) {
    # Process each day individually
    for (day in days) {
      # Get significant comparisons ONLY for this specific day among selected treatments
      day_comparisons <- Filter(function(x) {
        x$day == day && 
          x$treatment1 %in% selected_treatments && 
          x$treatment2 %in% selected_treatments
      }, significant_points)
      
      if (length(day_comparisons) > 0) {
        # For this day, get the selected treatments that have data
        day_treatments <- unique(data_subset$Treatment[data_subset$Day == day])
        
        # Create a matrix to track which treatments are significantly different
        n_treatments <- length(day_treatments)
        diff_matrix <- matrix(FALSE, nrow = n_treatments, ncol = n_treatments)
        rownames(diff_matrix) <- day_treatments
        colnames(diff_matrix) <- day_treatments
        
        # Fill in the difference matrix
        for (comp in day_comparisons) {
          if (comp$treatment1 %in% day_treatments && comp$treatment2 %in% day_treatments) {
            t1 <- comp$treatment1
            t2 <- comp$treatment2
            diff_matrix[t1, t2] <- TRUE
            diff_matrix[t2, t1] <- TRUE
          }
        }
        
        # Sort treatments by their mean values (descending) for this day
        treatment_means <- sapply(day_treatments, function(trt) {
          trt_data <- data_subset[data_subset$Treatment == trt & data_subset$Day == day, ]
          if (nrow(trt_data) > 0) trt_data[[y_var]] else NA
        })
        
        sorted_treatments <- names(sort(treatment_means, decreasing = TRUE))
        
        # Simple lettering approach: use distinct letters for distinct groups
        # First, initialize all to NA
        day_groups <- setNames(rep(NA, length(day_treatments)), day_treatments)
        
        # Start with letter A
        current_group <- "A"
        for (trt in sorted_treatments) {
          if (is.na(day_groups[trt])) {
            # Assign this treatment to the current group
            day_groups[trt] <- current_group
            
            # Find all treatments not significantly different from this one
            # and assign them to the same group if they don't have a group yet
            for (other_trt in day_treatments) {
              if (other_trt != trt && is.na(day_groups[other_trt]) && !diff_matrix[trt, other_trt]) {
                day_groups[other_trt] <- current_group
              }
            }
            
            # Move to next letter for next group
            current_group <- LETTERS[which(LETTERS == current_group) + 1]
          }
        }
        
        # Store the group labels for this day
        group_labels[[as.character(day)]] <- day_groups
      }
    }
  }
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_subset[data_subset$Treatment == treatment, ]
    
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
      
      # Add group label if available for this day
      day <- treatment_data$Day[j]
      if (as.character(day) %in% names(group_labels)) {
        day_group <- group_labels[[as.character(day)]]
        if (treatment %in% names(day_group)) {
          # Use black and bold text for group labels
          text(day + 0.15, treatment_data[[y_var]][j], 
               day_group[treatment], 
               col = "black",  # Always black letters
               font = 2,       # Bold font
               cex = 0.9, pos = 4)  # pos=4 means position text to the right
        }
      }
    }
  }
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}
# Function to create selective plot with black bold letters for significance groups
create_selective_plot_sig_days <- function(data, y_var, se_var, y_label, main_title, 
                                           selected_treatments, significant_points = NULL,
                                           selected_day = NULL, log_scale = FALSE) {
  # Filter data for selected treatments only
  data_subset <- data[data$Treatment %in% selected_treatments, ]
  
  # Filter for specific day if provided
  if (!is.null(selected_day)) {
    data_subset <- data_subset[data_subset$Day == selected_day, ]
    if (nrow(data_subset) == 0) {
      stop(paste("No data available for day", selected_day))
    }
    days <- selected_day
  } else {
    # Get unique days
    days <- sort(unique(data_subset$Day))
  }
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  max_value_with_se <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  y_max <- 1.1 * max_value_with_se  # More conservative adjustment to prevent empty space
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "green3", "blue", "cyan", "magenta", "orange")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  # For log scale, ensure all y values are positive
  log_param <- ""
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), ylim = c(max(y_min, 0.001), y_max),
       xlab = "Day", ylab = y_label, main = main_title,
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  all_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                      "Treatment3", "Treatment4", "Treatment5")
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Determine if we have controls and treatments
  has_controls <- any(selected_treatments %in% c("Control1", "Control2"))
  has_treatments <- any(selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5"))
  
  # Debugging: print significant_points to see what we're working with
  # cat("Number of significant points:", length(significant_points), "\n")
  # if (!is.null(significant_points) && length(significant_points) > 0) {
  #   for (i in 1:min(5, length(significant_points))) {
  #     cat("Point", i, ":", 
  #         "day =", significant_points[[i]]$day, 
  #         "treatments =", significant_points[[i]]$treatment1, "-", 
  #         significant_points[[i]]$treatment2, 
  #         "p-value =", significant_points[[i]]$p_value, "\n")
  #   }
  # }
  
  # Process significant points to generate letter groups for each day
  day_letter_groups <- list()
  
  # Step 1: For each day, extract the significant pairwise differences
  for (day in days) {
    # Get all combinations of treatments for this day
    day_data <- data_subset[data_subset$Day == day, ]
    day_treatments <- unique(day_data$Treatment)
    
    if (length(day_treatments) < 2) {
      next  # Skip days with fewer than 2 treatments
    }
    
    # Create a matrix to track significant differences (TRUE = significantly different)
    sig_diff_matrix <- matrix(FALSE, 
                              nrow = length(day_treatments), 
                              ncol = length(day_treatments),
                              dimnames = list(day_treatments, day_treatments))
    
    # Collect mean values for sorting
    means <- numeric(length(day_treatments))
    names(means) <- day_treatments
    
    for (t in day_treatments) {
      t_data <- day_data[day_data$Treatment == t, ]
      means[t] <- t_data[[y_var]]
    }
    
    # Fill matrix with significant differences from significant_points
    if (!is.null(significant_points) && length(significant_points) > 0) {
      for (point in significant_points) {
        # Only consider points for this specific day
        if (point$day == day) {
          t1 <- point$treatment1
          t2 <- point$treatment2
          
          # Check if both treatments are in our selected set
          if (t1 %in% day_treatments && t2 %in% day_treatments) {
            sig_diff_matrix[t1, t2] <- TRUE
            sig_diff_matrix[t2, t1] <- TRUE
          }
        }
      }
    }
    
    # Sort treatments by descending mean value
    sorted_treatments <- names(sort(means, decreasing = TRUE))
    
    # Start with all treatments in group "A"
    letter_groups <- rep("A", length(day_treatments))
    names(letter_groups) <- day_treatments
    
    # Assign letters A, B, C, etc. based on significant differences
    available_letters <- LETTERS
    current_letter_idx <- 1
    
    # Step 2: Assign letters to each treatment
    for (i in 1:length(sorted_treatments)) {
      treat_i <- sorted_treatments[i]
      
      # If this treatment already has a letter other than "A", skip it
      if (letter_groups[treat_i] != "A") {
        next
      }
      
      # Assign current letter
      letter_groups[treat_i] <- available_letters[current_letter_idx]
      
      # Find all treatments that are NOT significantly different from this one
      similar_treats <- c()
      for (j in 1:length(sorted_treatments)) {
        treat_j <- sorted_treatments[j]
        if (i != j && !sig_diff_matrix[treat_i, treat_j]) {
          similar_treats <- c(similar_treats, treat_j)
        }
      }
      
      # All similar treatments get the same letter
      if (length(similar_treats) > 0) {
        letter_groups[similar_treats] <- available_letters[current_letter_idx]
      }
      
      # Next letter for the next group
      current_letter_idx <- current_letter_idx + 1
    }
    
    # Store letter groups for this day
    day_letter_groups[[as.character(day)]] <- letter_groups
  }
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_subset[data_subset$Treatment == treatment, ]
    
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
      day <- treatment_data$Day[j]
      arrows(day, 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             day, 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.05, col = treatment_colors[idx], lwd = 1.5)
      
      # Add letter group label
      day_str <- as.character(day)
      if (day_str %in% names(day_letter_groups)) {
        day_group <- day_letter_groups[[day_str]]
        if (treatment %in% names(day_group)) {
          letter <- day_group[treatment]
          text(day + 0.2, treatment_data[[y_var]][j], 
               letter, 
               col = "black", # Always black
               font = 2,      # Bold
               cex = 0.9)
        }
      }
    }
  }
  
  # Add legend with organized controls and treatments sections
  if (has_controls && has_treatments) {
    # Separate controls and treatments in legend
    controls <- selected_treatments[selected_treatments %in% c("Control1", "Control2")]
    treatments <- selected_treatments[selected_treatments %in% c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")]
    
    control_indices <- match(controls, all_treatments)
    treatment_indices <- match(treatments, all_treatments)
    
    legend_items <- c("Controls:", controls, "Treatments:", treatments)
    legend_colors <- c(NA, treatment_colors[control_indices], NA, treatment_colors[treatment_indices])
    legend_pch <- c(NA, treatment_pch[control_indices], NA, treatment_pch[treatment_indices])
    legend_lty <- c(NA, rep(1, length(controls)), NA, rep(1, length(treatments)))
    legend_lwd <- c(NA, rep(2, length(controls)), NA, rep(2, length(treatments)))
    
    legend("topright", inset = c(-0.3, 0), 
           legend = legend_items, 
           col = legend_colors,
           pch = legend_pch,
           lty = legend_lty,
           lwd = legend_lwd,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  } else {
    # Just show the selected treatments without categories
    legend("topright", inset = c(-0.3, 0), 
           legend = selected_treatments, 
           col = treatment_colors[treatment_indices],
           pch = treatment_pch[treatment_indices],
           lty = 1, lwd = 2,
           cex = 0.9, 
           bty = "n",  # no box around legend
           xpd = TRUE) # allow plotting outside the plot region
  }
}
