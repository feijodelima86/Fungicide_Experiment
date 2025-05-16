analyze_fungicide_effect <- function(data, response_var, log_transform = FALSE, alpha = 0.05, log_base = 10) {
  # Make a copy of the data to avoid modifying the original
  analysis_data <- data
  
  # Check for required column
  if (!"Group" %in% names(analysis_data) && "treatment_group" %in% names(analysis_data)) {
    analysis_data$Group <- analysis_data$treatment_group
  }
  
  # Original response variable name for labeling
  original_var_name <- response_var
  
  # Handle log transformation with improved approach
  if (log_transform) {
    min_value <- min(analysis_data[[response_var]], na.rm = TRUE)
    
    # Automatically apply log+1 transformation for all data
    # This works for both positive and non-positive values
    if (log_base == 10) {
      analysis_data[[response_var]] <- log10(analysis_data[[response_var]] + abs(min_value) + 1)
      cat("Applied log10(x + |min| + 1) transform to", response_var, "\n")
      cat("Added constant of", abs(min_value) + 1, "to ensure all values are positive before log transform\n")
    } else if (log_base == "e" || log_base == "natural") {
      analysis_data[[response_var]] <- log(analysis_data[[response_var]] + abs(min_value) + 1)
      cat("Applied natural log(x + |min| + 1) transform to", response_var, "\n")
      cat("Added constant of", abs(min_value) + 1, "to ensure all values are positive before log transform\n")
    } else {
      # For any other base
      analysis_data[[response_var]] <- log(analysis_data[[response_var]] + abs(min_value) + 1, base = log_base)
      cat("Applied log base", log_base, "(x + |min| + 1) transform to", response_var, "\n")
      cat("Added constant of", abs(min_value) + 1, "to ensure all values are positive before log transform\n")
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
  
  # Always use gaussian with identity link for both transformed and untransformed data
  family_type <- gaussian(link = "identity")
  
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
    log_transform_constant = if(log_transform) abs(min_value) + 1 else NA,
    log_base = log_base,
    fungicide_conc = fungicide_conc
  )
  
  return(result)
}
# Simplified plotting functions for fungicide analysis results


# Function to generate just the data visualization plots (no scorecards)
plot_fungicide_data <- function(results, 
                                y_multiplier_linear = 1.1, 
                                y_multiplier_log = 1.05,
                                line_spacing_factor = 1.0) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Define treatment set - only treatments, no controls
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  
  # Create a side-by-side layout
  par(mfrow = c(1, 2))
  
  #------------------------------------------------------------------
  # PLOT 1: ALL TREATMENTS ONLY (ADJUSTED DATA)
  #------------------------------------------------------------------
  
  # Prepare data
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get unique days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset upper limit multiplier
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "", # No caption
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Define line styles for treatments only - darker grays and thicker lines
  # Treatment lines vary from T1 (dashed) to T5 (extra thick black)
  line_types <- c(2, 1, 1, 1, 1)  # T1 is dashed, others solid
  
  # Apply line spacing factor to line widths
  base_line_widths <- c(1.3, 1.6, 2.0, 2.5, 3.5)  # Original thick lines
  line_widths <- base_line_widths * line_spacing_factor
  
  # Treatment colors and point symbols
  treatment_colors <- c("gray65", "gray55", "gray35", "gray15", "black")
  treatment_pch <- c(17, 18, 19, 20, 21)  # Different point symbols
  
  # Plot each treatment in order of concentration (low to high)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line with custom styling
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = line_widths[i], lty = line_types[i])
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.3)
    
    # Add error bars (smaller and less prominent)
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.03, 
             col = adjustcolor(treatment_colors[i], alpha.f = 0.7), lwd = 0.8)
    }
  }
  
  # Add treatments as text labels at the end of each line
  # Get the last day's data for each treatment
  last_day <- max(days)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    # Extract the last day's data point
    last_point <- data_to_use[data_to_use$Treatment == treatment & data_to_use$Day == last_day, ]
    if (nrow(last_point) > 0) {
      # Add treatment label slightly to the right of the last point
      # Shorten treatment name to just T1, T2, etc.
      label <- paste0("T", substr(treatment, 10, 10))
      
      text(last_day + 0.3, last_point[[y_var]], label, 
           col = treatment_colors[i], cex = 0.9, adj = 0)
    }
  }
  
  # Add legend at top left
  legend("topleft",
         legend = c("T1", "T2", "T3", "T4", "T5"), 
         col = treatment_colors,
         pch = treatment_pch,
         lty = line_types,
         lwd = line_widths,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  #------------------------------------------------------------------
  # PLOT 2: ADJUSTED DATA BY CONCENTRATION
  #------------------------------------------------------------------
  
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- results$means_adjusted[results$means_adjusted$Treatment %in% treatments_to_include, ]
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Get unique concentrations
  concentrations <- sort(unique(plot_data$Concentration))
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_adj_var]] - plot_data[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_adj_var]] + plot_data[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up plotting area with adjusted margins
  par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(results$fungicide_conc[treatments_to_include])
  max_conc <- max(results$fungicide_conc[treatments_to_include])
  
  # Get unique days and set up colors
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Handle log scale appropriately
  # For x-axis (concentration), we always use log scale
  # For y-axis, check if values are appropriate for log scale
  log_param <- "x"  # Always log scale for x-axis (concentration)
  
  # Reset log parameters and multipliers for second plot
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use smaller multiplier for log scale
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with improved y-axis limits
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = paste0("Adjusted ", y_label), 
       main = "", # No caption
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1)
  
  # Add x-axis with treatment concentrations
  conc_for_axis <- c(0.054, 0.54, 5.4, 54, 540)
  axis_labels <- c("0.054", "0.54", "5.4", "54", "540")
  axis(1, at = conc_for_axis, labels = axis_labels)
  
  # Add grid
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray80")
  
  # Increase line widths to match the thickness in the first plot
  # Using maximum thickness from the first plot as a baseline
  enhanced_line_width <- max(line_widths)
  
  # Plot each day separately
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    
    # Sort by concentration
    day_data <- day_data[order(day_data$Concentration), ]
    if(nrow(day_data) == 0) next
    
    # Plot the line with enhanced thickness based on line spacing factor
    lines(day_data$Concentration, day_data[[results$mean_adj_var]], 
          col = day_colors[i], lwd = enhanced_line_width * line_spacing_factor)
    
    # Plot points
    for (j in 1:nrow(day_data)) {
      # Draw the point
      points(day_data$Concentration[j], day_data[[results$mean_adj_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
      
      # Add error bars
      arrows(day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
             day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
             angle = 90, code = 3, length = 0.05, col = day_colors[i], lwd = 1.5)
    }
  }
  
  # Add day labels at the right end of each line where possible
  # Find the rightmost point for each day
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    if(nrow(day_data) == 0) next
    
    # Get the point with maximum concentration
    max_conc_idx <- which.max(day_data$Concentration)
    if(length(max_conc_idx) > 0) {
      # Label just to the right of the last point
      text(day_data$Concentration[max_conc_idx] * 1.15, 
           day_data[[results$mean_adj_var]][max_conc_idx], 
           paste("Day", day), 
           col = day_colors[i], cex = 0.8, adj = 0)
    }
  }
  
  # Add legend at top left
  legend("topleft", 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         lty = 1, 
         lwd = 2,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  # Return color schemes for use in scorecards
  return(list(
    treatment_colors = treatment_colors,
    day_colors = day_colors,
    days = days,
    concentrations = concentrations
  ))
}

# Function to generate just the significance scorecards
plot_significance_scorecards <- function(results, 
                                         color_info = NULL,
                                         asterisk_size = 1.2) {
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # If color_info not provided, get defaults
  if (is.null(color_info)) {
    # Define treatment set for analysis
    selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
    
    # Get the unique days
    data_to_use <- results$means_adjusted
    days <- sort(unique(data_to_use$Day[data_to_use$Treatment %in% selected_treatments]))
    
    # Darker gray colors for all treatments
    treatment_colors <- c("gray65", "gray55", "gray35", "gray15", "black")
    
    # Create a gradient of blue colors from light to dark
    num_days <- length(days)
    day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
    
    # Get concentrations
    plot_data <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
    plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
    concentrations <- sort(unique(plot_data$Concentration))
  } else {
    # Use provided color info
    treatment_colors <- color_info$treatment_colors
    day_colors <- color_info$day_colors
    days <- color_info$days
    concentrations <- color_info$concentrations
  }
  
  # Define treatment set
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  
  # Create a side-by-side layout for scorecards
  par(mfrow = c(1, 2))
  
  #------------------------------------------------------------------
  # SCORECARD 1: TREATMENT SCORECARD (COMPACT VERSION)
  #------------------------------------------------------------------
  
  # Create margin for scorecard
  par(mar = c(2, 3, 3, 1))
  
  # Create an empty plot for the scorecard - tight and clean
  plot(0, 0, type = "n", 
       xlim = c(0, length(days) + 1), 
       ylim = c(0, length(selected_treatments) + 1),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  
  # Add title
  text(length(days)/2, length(selected_treatments) + 0.7, "Significance by Treatment and Day:", 
       adj = c(0.5, 0.5), font = 2, cex = 1.0)
  
  # Draw background grid for clarity
  for (i in 1:length(days)) {
    for (j in 1:length(selected_treatments)) {
      rect(i-0.4, j-0.4, i+0.4, j+0.4, border = "lightgray", col = NULL)
    }
  }
  
  # Draw day labels at the bottom
  for (i in 1:length(days)) {
    text(i, 0.3, paste0("D", days[i]), cex = 0.9)
  }
  
  # Draw treatment labels on the left
  for (j in 1:length(selected_treatments)) {
    text(0.3, j, paste0("T", j), cex = 0.9, col = treatment_colors[j], font = 2)
  }
  
  # Add significance asterisks in the grid
  for (i in 1:length(days)) {
    for (j in 1:length(selected_treatments)) {
      # Check if this combination is significant
      day <- days[i]
      treatment <- selected_treatments[j]
      lookup_key <- paste(day, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Get p-value for significance
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Add asterisk with treatment color
        text(i, j, asterisks, col = treatment_colors[j], cex = asterisk_size, font = 2)
      }
    }
  }
  
  #------------------------------------------------------------------
  # SCORECARD 2: CONCENTRATION SCORECARD (COMPACT VERSION)
  #------------------------------------------------------------------
  
  # Create margin for scorecard
  par(mar = c(2, 3, 3, 1))
  
  # Create an empty plot for the scorecard - tight and clean
  plot(0, 0, type = "n", 
       xlim = c(0, length(concentrations) + 1), 
       ylim = c(0, length(days) + 1),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
  
  # Add title
  text(length(concentrations)/2, length(days) + 0.7, "Significance by Day and Concentration:", 
       adj = c(0.5, 0.5), font = 2, cex = 1.0)
  
  # Draw background grid for clarity
  for (i in 1:length(concentrations)) {
    for (j in 1:length(days)) {
      rect(i-0.4, j-0.4, i+0.4, j+0.4, border = "lightgray", col = NULL)
    }
  }
  
  # Draw concentration labels at the bottom
  for (i in 1:length(concentrations)) {
    # Format concentration with fewer digits
    conc_label <- format(concentrations[i], scientific = FALSE, digits = 2)
    text(i, 0.3, conc_label, cex = 0.8, srt = 45, adj = c(1, 1))
  }
  
  # Draw day labels on the left
  for (j in 1:length(days)) {
    text(0.3, j, paste0("D", days[j]), cex = 0.9, col = day_colors[j], font = 2)
  }
  
  # Add significance asterisks in the grid
  for (i in 1:length(concentrations)) {
    conc <- concentrations[i]
    
    # Find the treatment that corresponds to this concentration
    treatment <- NULL
    for (trt in selected_treatments) {
      if (abs(results$fungicide_conc[trt] - conc) < 0.0001) {
        treatment <- trt
        break
      }
    }
    
    if (!is.null(treatment)) {
      # For each day, check if significant at this concentration
      for (j in 1:length(days)) {
        day <- days[j]
        
        # Check if this combination is significant
        lookup_key <- paste(day, treatment, sep = "_")
        
        if (!is.null(sig_lookup[[lookup_key]])) {
          # Get p-value for significance
          p_value <- sig_lookup[[lookup_key]]
          asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
          
          # Add asterisk with day color
          text(i, j, asterisks, col = day_colors[j], cex = asterisk_size, font = 2)
        }
      }
    }
  }
  
  # Reset the plotting parameters
  par(mfrow = c(1, 1))
}




# Plot treatment comparison - original vs adjusted data with improved y-axis scaling
plot_treatment_comparison <- function(results) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: ORIGINAL DATA ---
  
  # Get the unique treatments and days
  treatments <- levels(results$means_original$Treatment)
  days <- sort(unique(results$means_original$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(results$means_original[[results$mean_var]] - results$means_original[[results$se_var]], na.rm = TRUE)
  y_max <- max(results$means_original[[results$mean_var]] + results$means_original[[results$se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Adjust the upper limit multiplier based on scale type
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with appropriate axes
  # Use different multipliers for min and max
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = y_label, main = "Original Data",
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
    treatment_data <- results$means_original[results$means_original$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line
    lines(treatment_data$Day, treatment_data[[results$mean_var]], 
          col = treatment_colors[i], lwd = 2)
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[results$mean_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[results$mean_var]][j] - treatment_data[[results$se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[results$mean_var]][j] + treatment_data[[results$se_var]][j],
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
  
  # --- PLOT 2: ADJUSTED DATA ---
  
  # Get the unique treatments and days
  treatments <- levels(results$means_adjusted$Treatment)
  days <- sort(unique(results$means_adjusted$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(results$means_adjusted[[results$mean_adj_var]] - results$means_adjusted[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(results$means_adjusted[[results$mean_adj_var]] + results$means_adjusted[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset the upper limit multiplier
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with appropriate axes - using improved y-axis limits
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), main = "Acetone-Adjusted Data",
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
    treatment_data <- results$means_adjusted[results$means_adjusted$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line
    lines(treatment_data$Day, treatment_data[[results$mean_adj_var]], 
          col = treatment_colors[i], lwd = 2)
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[results$mean_adj_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.5)
    
    # Add error bars
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[results$mean_adj_var]][j] - treatment_data[[results$se_adj_var]][j],
             treatment_data$Day[j], 
             treatment_data[[results$mean_adj_var]][j] + treatment_data[[results$se_adj_var]][j],
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

# Plot selective treatment comparison with improved y-axis scaling
plot_selective_treatment_comparison <- function(results) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: CONTROLS AND HIGHEST TREATMENT ---
  
  # Define treatment set for first plot
  selected_treatments <- c("Control1", "Control2", "Treatment5")
  data_to_use <- results$means_original
  y_var <- results$mean_var
  se_var <- results$se_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Define upper limit multiplier based on scale type
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = y_label, main = "Controls and Highest Treatment",
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
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
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
      
      # If this is a significant point, add an asterisk (but not for control treatments)
      if (treatment != "Control1" && treatment != "Control2") {
        day_value <- treatment_data$Day[j]
        lookup_key <- paste(day_value, treatment, sep = "_")
        
        if (!is.null(sig_lookup[[lookup_key]])) {
          # Add bold asterisk next to the point
          p_value <- sig_lookup[[lookup_key]]
          asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
          
          # Adjust asterisk positioning based on the scale
          text_x <- treatment_data$Day[j] + 0.2
          
          # For log scale, use multiplicative offset, for linear use additive
          if (log_param == "y") {
            text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
          } else {
            text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.2
          }
          
          text(text_x, text_y, labels = asterisks, 
               col = treatment_colors[idx], font = 2, cex = 1.3)
        }
      }
    }
  }
  
  # Add legend
  legend("topright", inset = c(-0.3, 0), 
         legend = selected_treatments, 
         col = treatment_colors[treatment_indices],
         pch = treatment_pch[treatment_indices],
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
  
  # --- PLOT 2: ALL TREATMENTS (ADJUSTED DATA) ---
  
  # Define treatment set for second plot
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset upper limit multiplier
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "All Fungicide Treatments (Acetone-Adjusted)",
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Plot each treatment
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
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
      
      # If this is a significant point, add an asterisk
      day_value <- treatment_data$Day[j]
      lookup_key <- paste(day_value, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Add bold asterisk next to the point
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Adjust asterisk positioning based on the scale
        text_x <- treatment_data$Day[j] + 0.2
        
        # For log scale, use multiplicative offset, for linear use additive
        if (log_param == "y") {
          text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
        } else {
          text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.2
        }
        
        text(text_x, text_y, labels = asterisks, 
             col = treatment_colors[idx], font = 2, cex = 1.3)
      }
    }
  }
  
  # Add legend for treatments
  legend("topright", inset = c(-0.3, 0), 
         legend = selected_treatments, 
         col = treatment_colors[treatment_indices],
         pch = treatment_pch[treatment_indices],
         lty = 1, lwd = 2,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Plot concentration response with improved y-axis scaling
plot_concentration_response <- function(results) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: ORIGINAL DATA BY CONCENTRATION ---
  
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- results$means_original[results$means_original$Treatment %in% treatments_to_include, ]
  
  # Get unique days
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_var]] - plot_data[[results$se_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_var]] + plot_data[[results$se_var]], na.rm = TRUE)
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(results$fungicide_conc[treatments_to_include])
  max_conc <- max(results$fungicide_conc[treatments_to_include])
  
  # Handle log scale appropriately
  # For x-axis (concentration), we always use log scale
  # For y-axis, check if values are appropriate for log scale
  log_param <- "x"  # Always log scale for x-axis (concentration)
  
  # Define upper limit multiplier based on scale type
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with improved y-axis limits
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = y_label, 
       main = "Original Data by Fungicide Concentration",
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
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
    lines(day_data$Concentration, day_data[[results$mean_var]], 
          col = day_colors[i], lwd = 2)
    
    # Plot points
    for (j in 1:nrow(day_data)) {
      # Draw the point
      points(day_data$Concentration[j], day_data[[results$mean_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
    }
    
    # Add error bars
    for (j in 1:nrow(day_data)) {
      arrows(day_data$Concentration[j], 
             day_data[[results$mean_var]][j] - day_data[[results$se_var]][j],
             day_data$Concentration[j], 
             day_data[[results$mean_var]][j] + day_data[[results$se_var]][j],
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
  
  # --- PLOT 2: ADJUSTED DATA BY CONCENTRATION WITH SIGNIFICANCE ---
  
  # Filter out control groups - include only Treatment1-5
  data_filtered <- results$means_adjusted[results$means_adjusted$Treatment %in% treatments_to_include, ]
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_adj_var]] - plot_data[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_adj_var]] + plot_data[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Reset log parameters and multipliers for second plot
  log_param <- "x"  # Always log scale for x-axis (concentration)
  upper_limit_multiplier <- 1.1  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use smaller multiplier for log scale
      upper_limit_multiplier <- 1.05
    }
  }
  
  # Create empty plot with improved y-axis limits
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = paste0("Adjusted ", y_label), 
       main = "Acetone-Adjusted Data by Fungicide Concentration\nAsterisks (*) indicate significant differences from Control2",
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add x-axis with treatment concentrations
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
    lines(day_data$Concentration, day_data[[results$mean_adj_var]], 
          col = day_colors[i], lwd = 2)
    
    # Plot points and add significance markers if needed
    for (j in 1:nrow(day_data)) {
      # Check if this point is significant
      is_sig <- FALSE
      if (!is.null(results$significant_points) && length(results$significant_points) > 0) {
        for (sig_point in results$significant_points) {
          if (sig_point$day == day && 
              sig_point$concentration == day_data$Concentration[j]) {
            is_sig <- TRUE
            break
          }
        }
      }
      
      # Draw the point
      points(day_data$Concentration[j], day_data[[results$mean_adj_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
      
      # Add significance marker if needed with improved positioning
      if (is_sig) {
        # Adjust the position of the asterisk based on scale
        if (log_param == "xy") {
          # For log scale, use multiplicative offset
          text_y <- day_data[[results$mean_adj_var]][j] * 1.03 + day_data[[results$se_adj_var]][j] * 1.1
        } else {
          # For linear scale, use additive offset
          text_y <- day_data[[results$mean_adj_var]][j] + 1.3 * day_data[[results$se_adj_var]][j]
        }
        
        text(day_data$Concentration[j], text_y, "*", cex = 1.5, font = 2)  # Bold asterisk
      }
    }
    
    # Add error bars
    for (j in 1:nrow(day_data)) {
      arrows(day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
             day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
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
}

# Plot treatment comparison with jittered points and customizable y-axis scaling
plot_jittered_treatment_comparison <- function(results, jitter_width = 0.5, 
                                               y_multiplier_linear = 1.1, 
                                               y_multiplier_log = 1.01) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: CONTROLS AND HIGHEST TREATMENT ---
  
  # Define treatment set for first plot
  selected_treatments <- c("Control1", "Control2", "Treatment5")
  data_to_use <- results$means_original
  y_var <- results$mean_var
  se_var <- results$se_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Treatment colors and point symbols
  treatment_colors <- c("black", "red", "#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725")
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Use the custom multiplier parameter based on scale type
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use custom log multiplier parameter
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes and customized y-limits
  plot(NULL, xlim = range(days) + c(-0.5, 0.5), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = y_label, main = "Controls and Highest Treatment",
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
  
  # Set jitter offsets for each treatment
  jitter_offsets <- seq(-jitter_width, jitter_width, length.out = length(selected_treatments))
  
  # Plot each treatment (NO LINES, just points with jitter and wider error bars)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Get the index of this treatment in the full treatment list for consistent colors
    idx <- treatment_indices[i]
    
    # Sort by day
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Apply jitter to the x-coordinates
    jittered_x <- treatment_data$Day + jitter_offsets[i]
    
    # Plot the points with jitter
    points(jittered_x, treatment_data[[y_var]], 
           col = treatment_colors[idx], pch = treatment_pch[idx], cex = 1.7)
    
    # Add wider error bars (box-and-whisker style)
    for (j in 1:nrow(treatment_data)) {
      # Vertical error bar
      arrows(jittered_x[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             jittered_x[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.08, col = treatment_colors[idx], lwd = 2)
      
      # Add horizontal caps at ends of error bars (for box-and-whisker appearance)
      # Top cap
      segments(jittered_x[j] - 0.05, treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
               jittered_x[j] + 0.05, treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
               col = treatment_colors[idx], lwd = 2)
      
      # Bottom cap
      segments(jittered_x[j] - 0.05, treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
               jittered_x[j] + 0.05, treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
               col = treatment_colors[idx], lwd = 2)
      
      # If this is a significant point, add a bold black asterisk (but not for control treatments)
      if (treatment != "Control1" && treatment != "Control2") {
        day_value <- treatment_data$Day[j]
        lookup_key <- paste(day_value, treatment, sep = "_")
        
        if (!is.null(sig_lookup[[lookup_key]])) {
          # Add bold black asterisk with improved positioning
          p_value <- sig_lookup[[lookup_key]]
          asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
          
          # Position the asterisk based on the scale
          text_x <- jittered_x[j]
          
          # For log scale, use a smaller multiplicative offset
          if (log_param == "y") {
            text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
          } else {
            # For linear scale, use the original offset
            text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.4
          }
          
          text(text_x, text_y, labels = asterisks, 
               col = "black", font = 2, cex = 1.5)
        }
      }
    }
  }
  
  # Add legend
  legend("topright", inset = c(-0.3, 0), 
         legend = selected_treatments, 
         col = treatment_colors[treatment_indices],
         pch = treatment_pch[treatment_indices],
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
  
  # --- PLOT 2: ALL TREATMENTS (ADJUSTED DATA) ---
  
  # Define treatment set for second plot
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset the upper limit multiplier, using custom parameter
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use custom log multiplier parameter
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes and improved y-limits
  plot(NULL, xlim = range(days) + c(-0.5, 0.5), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "All Fungicide Treatments (Acetone-Adjusted)",
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Get treatment indices for color and point shape matching
  treatment_indices <- match(selected_treatments, all_treatments)
  
  # Set jitter offsets for each treatment
  jitter_offsets <- seq(-jitter_width, jitter_width, length.out = length(selected_treatments))
  
  # Plot each treatment (NO LINES, just points with jitter and wider error bars)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Get the index of this treatment in the full treatment list for consistent colors
    idx <- treatment_indices[i]
    
    # Sort by day
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Apply jitter to the x-coordinates
    jittered_x <- treatment_data$Day + jitter_offsets[i]
    
    # Plot the points with jitter
    points(jittered_x, treatment_data[[y_var]], 
           col = treatment_colors[idx], pch = treatment_pch[idx], cex = 1.7)
    
    # Add wider error bars (box-and-whisker style)
    for (j in 1:nrow(treatment_data)) {
      # Vertical error bar
      arrows(jittered_x[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             jittered_x[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.08, col = treatment_colors[idx], lwd = 2)
      
      # Add horizontal caps at ends of error bars (for box-and-whisker appearance)
      # Top cap
      segments(jittered_x[j] - 0.05, treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
               jittered_x[j] + 0.05, treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
               col = treatment_colors[idx], lwd = 2)
      
      # Bottom cap
      segments(jittered_x[j] - 0.05, treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
               jittered_x[j] + 0.05, treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
               col = treatment_colors[idx], lwd = 2)
      
      # If this is a significant point, add an asterisk with improved positioning
      day_value <- treatment_data$Day[j]
      lookup_key <- paste(day_value, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Add bold black asterisk
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Position the asterisk based on scale
        text_x <- jittered_x[j]
        
        # For log scale, use a smaller multiplicative offset
        if (log_param == "y") {
          text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
        } else {
          # For linear scale, use the original offset
          text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.4
        }
        
        text(text_x, text_y, labels = asterisks, 
             col = "black", font = 2, cex = 1.5)
      }
    }
  }
  
  # Add legend for treatments
  legend("topright", inset = c(-0.3, 0), 
         legend = selected_treatments, 
         col = treatment_colors[treatment_indices],
         pch = treatment_pch[treatment_indices],
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Plot jittered concentration response with customizable y-axis scaling
plot_jittered_concentration_response <- function(results, jitter_width = 0.5,
                                                 y_multiplier_linear = 1.1,
                                                 y_multiplier_log = 1.01) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: ORIGINAL DATA BY CONCENTRATION ---
  
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- results$means_original[results$means_original$Treatment %in% treatments_to_include, ]
  
  # Get unique days
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_var]] - plot_data[[results$se_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_var]] + plot_data[[results$se_var]], na.rm = TRUE)
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(results$fungicide_conc[treatments_to_include])
  max_conc <- max(results$fungicide_conc[treatments_to_include])
  
  # Handle log scale appropriately
  # For x-axis (concentration), we always use log scale
  # For y-axis, check if values are appropriate for log scale
  log_param <- "x"  # Always log scale for x-axis (concentration)
  
  # Use custom multiplier parameter based on scale type
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use custom log multiplier parameter
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with improved y-axis limits using the custom multiplier
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 2),  # Expanded to accommodate jitter
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = y_label, 
       main = "Original Data by Fungicide Concentration",
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
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
    
    # Jitter the x-coordinates (for log scale we need to use multiplication)
    jitter_factor <- 10^(jitter_width * (i - (num_days + 1)/2) / ((num_days + 1)/2))
    jittered_x <- day_data$Concentration * jitter_factor
    
    # Plot points with jitter (NO LINES)
    for (j in 1:nrow(day_data)) {
      # Draw the point
      points(jittered_x[j], day_data[[results$mean_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.7)
      
      # Add box-and-whisker style error bars
      # Vertical error bar
      arrows(jittered_x[j], 
             day_data[[results$mean_var]][j] - day_data[[results$se_var]][j],
             jittered_x[j], 
             day_data[[results$mean_var]][j] + day_data[[results$se_var]][j],
             angle = 90, code = 3, length = 0.08, col = day_colors[i], lwd = 2)
      
      # Due to log scale on x-axis, we need to calculate cap width proportionally
      cap_width <- jittered_x[j] * 0.1  # 10% of x value
      
      # Top cap
      segments(jittered_x[j] - cap_width, day_data[[results$mean_var]][j] + day_data[[results$se_var]][j],
               jittered_x[j] + cap_width, day_data[[results$mean_var]][j] + day_data[[results$se_var]][j],
               col = day_colors[i], lwd = 2)
      
      # Bottom cap
      segments(jittered_x[j] - cap_width, day_data[[results$mean_var]][j] - day_data[[results$se_var]][j],
               jittered_x[j] + cap_width, day_data[[results$mean_var]][j] - day_data[[results$se_var]][j],
               col = day_colors[i], lwd = 2)
    }
  }
  
  # Add legend for days
  legend("topright", inset = c(-0.3, 0), 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
  
  # --- PLOT 2: ADJUSTED DATA BY CONCENTRATION WITH SIGNIFICANCE ---
  
  # Filter out control groups - include only Treatment1-5
  data_filtered <- results$means_adjusted[results$means_adjusted$Treatment %in% treatments_to_include, ]
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_adj_var]] - plot_data[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_adj_var]] + plot_data[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up plotting area
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Reset log parameters and multipliers for second plot
  log_param <- "x"  # Always log scale for x-axis (concentration)
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use custom log multiplier parameter
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with improved y-axis limits using the custom multiplier
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 2),  # Expanded to accommodate jitter
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = paste0("Adjusted ", y_label), 
       main = "Acetone-Adjusted Data by Fungicide Concentration",
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add x-axis with treatment concentrations
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
    
    # Jitter the x-coordinates (for log scale we need to use multiplication)
    jitter_factor <- 10^(jitter_width * (i - (num_days + 1)/2) / ((num_days + 1)/2))
    jittered_x <- day_data$Concentration * jitter_factor
    
    # Plot points and add significance markers if needed (NO LINES)
    for (j in 1:nrow(day_data)) {
      # Draw the point
      points(jittered_x[j], day_data[[results$mean_adj_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.7)
      
      # Add box-and-whisker style error bars
      # Vertical error bar
      arrows(jittered_x[j], 
             day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
             jittered_x[j], 
             day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
             angle = 90, code = 3, length = 0.08, col = day_colors[i], lwd = 2)
      
      # Due to log scale on x-axis, we need to calculate cap width proportionally
      cap_width <- jittered_x[j] * 0.1  # 10% of x value
      
      # Top cap
      segments(jittered_x[j] - cap_width, day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
               jittered_x[j] + cap_width, day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
               col = day_colors[i], lwd = 2)
      
      # Bottom cap
      segments(jittered_x[j] - cap_width, day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
               jittered_x[j] + cap_width, day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
               col = day_colors[i], lwd = 2)
      
      # Check if this point is significant
      treatment <- day_data$Treatment[j]
      lookup_key <- paste(day, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Add significance marker with improved positioning
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Adjust asterisk position based on scale
        if (log_param == "xy") {
          # For log scale, use multiplicative offset
          text_y <- day_data[[results$mean_adj_var]][j] * 1.03 + day_data[[results$se_adj_var]][j] * 1.1
        } else {
          # For linear scale, use additive offset
          text_y <- day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j] * 1.5
        }
        
        text(jittered_x[j], text_y, asterisks, col = "black", font = 2, cex = 1.5)  # Bold black asterisk
      }
    }
  }
  
  # Add legend for days
  legend("topright", inset = c(-0.3, 0), 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         cex = 0.9, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Alternative plot_selective_treatment_comparison with concentration-based line styling
plot_selective_treatment_comparison_alt <- function(results, y_multiplier_linear = 1.1, 
                                                    y_multiplier_log = 1.05) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # Setup layout for two plots
  par(mfrow = c(2, 1))
  
  # --- PLOT 1: CONTROLS AND ALL TREATMENTS (ORIGINAL DATA) ---
  
  # Define treatment set for first plot - include all treatments
  selected_treatments <- c("Control1", "Control2", "Treatment1", "Treatment2", 
                           "Treatment3", "Treatment4", "Treatment5")
  data_to_use <- results$means_original
  y_var <- results$mean_var
  se_var <- results$se_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Define upper limit multiplier based on scale type
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = y_label, main = "Original Data - All Treatments",
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Define custom line types, widths, and colors based on treatment concentration
  # Controls keep original colors (black, red)
  # Treatments use darker grays and thicker lines overall
  line_types <- c(1, 1, 2, 1, 1, 1, 1)  # T1 is dashed (2), others solid (1)
  line_widths <- c(1, 1, 1.3, 1.6, 2.0, 2.5, 3.5)  # Overall thicker lines for treatments
  
  # Use original colors for controls, darker grayscale for treatments
  treatment_colors <- c("black", "red", "gray65", "gray55", "gray35", "gray15", "black")  
  treatment_pch <- c(15, 16, 17, 18, 19, 20, 21)  # Different point symbols for each treatment
  
  # Plot each treatment in order of concentration (low to high)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line with custom styling
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = line_widths[i], lty = line_types[i])
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.3)
    
    # Add error bars (smaller and less prominent)
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.03, 
             col = adjustcolor(treatment_colors[i], alpha.f = 0.7), lwd = 0.8)
      
      # If this is a significant point, add an asterisk (but not for control treatments)
      if (treatment != "Control1" && treatment != "Control2") {
        day_value <- treatment_data$Day[j]
        lookup_key <- paste(day_value, treatment, sep = "_")
        
        if (!is.null(sig_lookup[[lookup_key]])) {
          # Add asterisk next to the point
          p_value <- sig_lookup[[lookup_key]]
          asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
          
          # Adjust asterisk positioning based on the scale
          text_x <- treatment_data$Day[j] + 0.2
          
          # For log scale, use multiplicative offset, for linear use additive
          if (log_param == "y") {
            text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
          } else {
            text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.2
          }
          
          text(text_x, text_y, labels = asterisks, 
               col = treatment_colors[i], font = 2, cex = 1.2)
        }
      }
    }
  }
  
  # Add treatments as text labels at the end of each line
  # Get the last day's data for each treatment
  last_day <- max(days)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    # Extract the last day's data point
    last_point <- data_to_use[data_to_use$Treatment == treatment & data_to_use$Day == last_day, ]
    if (nrow(last_point) > 0) {
      # Add treatment label slightly to the right of the last point
      # Shorten treatment name to just T1, T2, etc. for treatments
      label <- treatment
      if (startsWith(treatment, "Treatment")) {
        label <- paste0("T", substr(treatment, 10, 10))
      }
      
      text(last_day + 0.3, last_point[[y_var]], label, 
           col = treatment_colors[i], cex = 0.9, adj = 0)
    }
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.3, 0), 
         legend = c("Controls:", "Control 1", "Control 2", "Treatments:", "T1 (lowest)", "T2", "T3", "T4", "T5 (highest)"), 
         col = c(NA, "black", "red", NA, "black", "black", "black", "black", "black"),
         pch = c(NA, 15, 16, NA, 17, 18, 19, 20, 21),
         lty = c(NA, 1, 1, NA, 1, 1, 1, 1, 1),
         lwd = c(NA, 1, 1, NA, 0.8, 1.2, 1.6, 2.0, 2.5),
         cex = 0.8, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
  
  # --- PLOT 2: ALL TREATMENTS ONLY (ADJUSTED DATA) ---
  
  # Define treatment set for second plot - only treatments, no controls
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 10))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset upper limit multiplier
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "Acetone-Adjusted Data - Treatments Only",
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Define line styles for treatments only - darker grays and thicker lines
  # Treatment lines vary from T1 (dashed) to T5 (extra thick black)
  line_types <- c(2, 1, 1, 1, 1)  # T1 is dashed, others solid
  line_widths <- c(1.3, 1.6, 2.0, 2.5, 3.5)  # Overall thicker lines
  
  # Darker gray colors for all treatments
  treatment_colors <- c("gray65", "gray55", "gray35", "gray15", "black")
  treatment_pch <- c(17, 18, 19, 20, 21)  # Different point symbols
  
  # Plot each treatment in order of concentration (low to high)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line with custom styling
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = line_widths[i], lty = line_types[i])
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.3)
    
    # Add error bars (smaller and less prominent)
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.03, 
             col = adjustcolor(treatment_colors[i], alpha.f = 0.7), lwd = 0.8)
      
      # If this is a significant point, add an asterisk
      day_value <- treatment_data$Day[j]
      lookup_key <- paste(day_value, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Add asterisk with improved positioning
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Adjust asterisk positioning based on the scale
        text_x <- treatment_data$Day[j] + 0.2
        
        # For log scale, use multiplicative offset, for linear use additive
        if (log_param == "y") {
          text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
        } else {
          text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.2
        }
        
        text(text_x, text_y, labels = asterisks, 
             col = treatment_colors[i], font = 2, cex = 1.2)
      }
    }
  }
  
  # Add treatments as text labels at the end of each line
  # Get the last day's data for each treatment
  last_day <- max(days)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    # Extract the last day's data point
    last_point <- data_to_use[data_to_use$Treatment == treatment & data_to_use$Day == last_day, ]
    if (nrow(last_point) > 0) {
      # Add treatment label slightly to the right of the last point
      # Shorten treatment name to just T1, T2, etc.
      label <- paste0("T", substr(treatment, 10, 10))
      
      text(last_day + 0.3, last_point[[y_var]], label, 
           col = treatment_colors[i], cex = 0.9, adj = 0)
    }
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.3, 0), 
         legend = c("Treatments:", "T1", "T2", "T3", "T4", "T5"), 
         col = c(NA, treatment_colors),
         pch = c(NA, treatment_pch),
         lty = c(NA, line_types),
         lwd = c(NA, line_widths),
         cex = 0.8, 
         bty = "n",  # no box around legend
         xpd = TRUE) # allow plotting outside the plot region
}

# Function to display side-by-side plots of treatments and concentration response
plot_treatment_and_concentration <- function(results, y_multiplier_linear = 1.1, 
                                             y_multiplier_log = 1.05) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Extract significant point information
  significant_points <- results$significant_points
  
  # Create lookup table for significant points
  sig_lookup <- list()
  if (length(significant_points) > 0) {
    for (point in significant_points) {
      key <- paste(point$day, point$treatment, sep = "_")
      sig_lookup[[key]] <- point$p_value
    }
  }
  
  # Setup layout for side-by-side plots
  par(mfrow = c(1, 2))
  
  #------------------------------------------------------------------
  # PLOT 1: ALL TREATMENTS ONLY (ADJUSTED DATA) - from plot_selective_treatment_comparison_alt
  #------------------------------------------------------------------
  
  # Define treatment set for first plot - only treatments, no controls
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get the unique treatments and days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 4, 2))  # c(bottom, left, top, right) - reduced right margin
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset upper limit multiplier
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "Acetone-Adjusted Data - Treatments Only",
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Define line styles for treatments only - darker grays and thicker lines
  # Treatment lines vary from T1 (dashed) to T5 (extra thick black)
  line_types <- c(2, 1, 1, 1, 1)  # T1 is dashed, others solid
  line_widths <- c(1.3, 1.6, 2.0, 2.5, 3.5)  # Overall thicker lines
  
  # Darker gray colors for all treatments
  treatment_colors <- c("gray65", "gray55", "gray35", "gray15", "black")
  treatment_pch <- c(17, 18, 19, 20, 21)  # Different point symbols
  
  # Plot each treatment in order of concentration (low to high)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line with custom styling
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = line_widths[i], lty = line_types[i])
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.3)
    
    # Add error bars (smaller and less prominent)
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.03, 
             col = adjustcolor(treatment_colors[i], alpha.f = 0.7), lwd = 0.8)
      
      # If this is a significant point, add an asterisk
      day_value <- treatment_data$Day[j]
      lookup_key <- paste(day_value, treatment, sep = "_")
      
      if (!is.null(sig_lookup[[lookup_key]])) {
        # Add asterisk with improved positioning
        p_value <- sig_lookup[[lookup_key]]
        asterisks <- if (p_value < 0.001) "***" else if (p_value < 0.01) "**" else "*"
        
        # Adjust asterisk positioning based on the scale
        text_x <- treatment_data$Day[j] + 0.2
        
        # For log scale, use multiplicative offset, for linear use additive
        if (log_param == "y") {
          text_y <- treatment_data[[y_var]][j] * 1.03 + treatment_data[[se_var]][j] * 1.1
        } else {
          text_y <- treatment_data[[y_var]][j] + treatment_data[[se_var]][j] * 1.2
        }
        
        text(text_x, text_y, labels = asterisks, 
             col = treatment_colors[i], font = 2, cex = 1.2)
      }
    }
  }
  
  # Add treatments as text labels at the end of each line
  # Get the last day's data for each treatment
  last_day <- max(days)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    # Extract the last day's data point
    last_point <- data_to_use[data_to_use$Treatment == treatment & data_to_use$Day == last_day, ]
    if (nrow(last_point) > 0) {
      # Add treatment label slightly to the right of the last point
      # Shorten treatment name to just T1, T2, etc.
      label <- paste0("T", substr(treatment, 10, 10))
      
      text(last_day + 0.3, last_point[[y_var]], label, 
           col = treatment_colors[i], cex = 0.9, adj = 0)
    }
  }
  
  # Add compact legend in the upper left
  legend("topleft",
         legend = c("T1", "T2", "T3", "T4", "T5"), 
         col = treatment_colors,
         pch = treatment_pch,
         lty = line_types,
         lwd = line_widths,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  #------------------------------------------------------------------
  # PLOT 2: ADJUSTED DATA BY CONCENTRATION (from plot_concentration_response)
  #------------------------------------------------------------------
  
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- results$means_adjusted[results$means_adjusted$Treatment %in% treatments_to_include, ]
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_adj_var]] - plot_data[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_adj_var]] + plot_data[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up plotting area with adjusted margins
  par(mar = c(5, 5, 4, 2))  # c(bottom, left, top, right)
  
  # Get unique days and set up colors
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Get min/max concentration values
  min_conc <- min(results$fungicide_conc[treatments_to_include])
  max_conc <- max(results$fungicide_conc[treatments_to_include])
  
  # Handle log scale appropriately
  # For x-axis (concentration), we always use log scale
  # For y-axis, check if values are appropriate for log scale
  log_param <- "x"  # Always log scale for x-axis (concentration)
  
  # Reset log parameters and multipliers for second plot
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use smaller multiplier for log scale
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with improved y-axis limits
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = paste0("Adjusted ", y_label), 
       main = "Acetone-Adjusted Data by Concentration",
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  
  # Add x-axis with treatment concentrations
  conc_for_axis <- c(0.054, 0.54, 5.4, 54, 540)
  axis_labels <- c("0.054", "0.54", "5.4", "54", "540")
  axis(1, at = conc_for_axis, labels = axis_labels)
  
  # Add grid
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray80")
  
  # Increase line widths to match the thickness in the first plot
  # Using maximum thickness from the first plot (3.5) as a baseline
  enhanced_line_width <- 3.5
  
  # Plot each day separately
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    
    # Sort by concentration
    day_data <- day_data[order(day_data$Concentration), ]
    if(nrow(day_data) == 0) next
    
    # Plot the line with enhanced thickness
    lines(day_data$Concentration, day_data[[results$mean_adj_var]], 
          col = day_colors[i], lwd = enhanced_line_width)
    
    # Plot points and add significance markers if needed
    for (j in 1:nrow(day_data)) {
      # Check if this point is significant
      is_sig <- FALSE
      if (!is.null(results$significant_points) && length(results$significant_points) > 0) {
        for (sig_point in results$significant_points) {
          if (sig_point$day == day && 
              sig_point$concentration == day_data$Concentration[j]) {
            is_sig <- TRUE
            break
          }
        }
      }
      
      # Draw the point
      points(day_data$Concentration[j], day_data[[results$mean_adj_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
      
      # Add significance marker if needed with improved positioning
      if (is_sig) {
        # Adjust the position of the asterisk based on scale
        if (log_param == "xy") {
          # For log scale, use multiplicative offset
          text_y <- day_data[[results$mean_adj_var]][j] * 1.03 + day_data[[results$se_adj_var]][j] * 1.1
        } else {
          # For linear scale, use additive offset
          text_y <- day_data[[results$mean_adj_var]][j] + 1.3 * day_data[[results$se_adj_var]][j]
        }
        
        text(day_data$Concentration[j], text_y, "*", cex = 1.5, font = 2)  # Bold asterisk
      }
    }
    
    # Add error bars
    for (j in 1:nrow(day_data)) {
      arrows(day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
             day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
             angle = 90, code = 3, length = 0.05, col = day_colors[i], lwd = 1.5)
    }
  }
  
  # Add day labels at the right end of each line where possible
  # Find the rightmost point for each day
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    if(nrow(day_data) == 0) next
    
    # Get the point with maximum concentration
    max_conc_idx <- which.max(day_data$Concentration)
    if(length(max_conc_idx) > 0) {
      # Label just to the right of the last point
      text(day_data$Concentration[max_conc_idx] * 1.15, 
           day_data[[results$mean_adj_var]][max_conc_idx], 
           paste("Day", day), 
           col = day_colors[i], cex = 0.8, adj = 0)
    }
  }
  
  # Add compact legend at the upper left
  legend("topleft", 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         lty = 1, 
         lwd = 2,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  # Reset the plotting parameters
  par(mfrow = c(1, 1))
}


# Function to generate just the data visualization plots (no scorecards)
plot_fungicide_data <- function(results, 
                                y_multiplier_linear = 1.1, 
                                y_multiplier_log = 1.05,
                                line_spacing_factor = 1.0) {
  # Get log_scale setting from results
  log_scale <- results$log_transform
  
  # Get y-label based on whether log transform was applied
  y_label <- if(log_scale) {
    paste0(results$original_var_name, " (log", results$log_base, ")")
  } else {
    results$original_var_name
  }
  
  # Define treatment set - only treatments, no controls
  selected_treatments <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  
  # Create a side-by-side layout
  par(mfrow = c(1, 2))
  
  #------------------------------------------------------------------
  # PLOT 1: ALL TREATMENTS ONLY (ADJUSTED DATA)
  #------------------------------------------------------------------
  
  # Prepare data
  data_to_use <- results$means_adjusted
  y_var <- results$mean_adj_var
  se_var <- results$se_adj_var
  
  # Filter data for selected treatments only
  data_subset <- data_to_use[data_to_use$Treatment %in% selected_treatments, ]
  
  # Get unique days
  days <- sort(unique(data_subset$Day))
  
  # Determine the y-axis limits including error bars
  y_min <- min(data_subset[[y_var]] - data_subset[[se_var]], na.rm = TRUE)
  y_max <- max(data_subset[[y_var]] + data_subset[[se_var]], na.rm = TRUE)
  
  # Set up an empty plot with appropriate margins
  par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
  
  # Handle log scale for y-axis when needed
  log_param <- ""
  # Reset upper limit multiplier
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale with non-positive values. Using linear scale instead.")
    } else {
      log_param <- "y"
      # Use smaller multiplier for log scale to avoid excessive space
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with appropriate axes
  plot(NULL, xlim = range(days), 
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Day", ylab = paste0("Adjusted ", y_label), 
       main = "", # No caption
       xaxt = "n", # we'll add custom x-axis
       log = log_param,
       cex.lab = 1.2, cex.axis = 1.1)
  
  # Add custom x-axis with only the specific days
  axis(1, at = days)
  
  # Add grid for readability
  grid(nx = NA, ny = NULL, lty = 2, col = "gray80")
  
  # Define line styles for treatments only - darker grays and thicker lines
  # Treatment lines vary from T1 (dashed) to T5 (extra thick black)
  line_types <- c(2, 1, 1, 1, 1)  # T1 is dashed, others solid
  
  # Apply line spacing factor to line widths
  base_line_widths <- c(1.3, 1.6, 2.0, 2.5, 3.5)  # Original thick lines
  line_widths <- base_line_widths * line_spacing_factor
  
  # Treatment colors and point symbols
  treatment_colors <- c("gray65", "gray55", "gray35", "gray15", "black")
  treatment_pch <- c(17, 18, 19, 20, 21)  # Different point symbols
  
  # Plot each treatment in order of concentration (low to high)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    treatment_data <- data_to_use[data_to_use$Treatment == treatment, ]
    
    # Sort by day to ensure lines connect points in order
    treatment_data <- treatment_data[order(treatment_data$Day), ]
    
    # Plot the line with custom styling
    lines(treatment_data$Day, treatment_data[[y_var]], 
          col = treatment_colors[i], lwd = line_widths[i], lty = line_types[i])
    
    # Plot the points
    points(treatment_data$Day, treatment_data[[y_var]], 
           col = treatment_colors[i], pch = treatment_pch[i], cex = 1.3)
    
    # Add error bars (smaller and less prominent)
    for (j in 1:nrow(treatment_data)) {
      arrows(treatment_data$Day[j], 
             treatment_data[[y_var]][j] - treatment_data[[se_var]][j],
             treatment_data$Day[j], 
             treatment_data[[y_var]][j] + treatment_data[[se_var]][j],
             angle = 90, code = 3, length = 0.03, 
             col = adjustcolor(treatment_colors[i], alpha.f = 0.7), lwd = 0.8)
    }
  }
  
  # Add treatments as text labels at the end of each line
  # Get the last day's data for each treatment
  last_day <- max(days)
  for (i in 1:length(selected_treatments)) {
    treatment <- selected_treatments[i]
    # Extract the last day's data point
    last_point <- data_to_use[data_to_use$Treatment == treatment & data_to_use$Day == last_day, ]
    if (nrow(last_point) > 0) {
      # Add treatment label slightly to the right of the last point
      # Shorten treatment name to just T1, T2, etc.
      label <- paste0("T", substr(treatment, 10, 10))
      
      text(last_day + 0.3, last_point[[y_var]], label, 
           col = treatment_colors[i], cex = 0.9, adj = 0)
    }
  }
  
  # Add legend at top left
  legend("topleft",
         legend = c("T1", "T2", "T3", "T4", "T5"), 
         col = treatment_colors,
         pch = treatment_pch,
         lty = line_types,
         lwd = line_widths,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  #------------------------------------------------------------------
  # PLOT 2: ADJUSTED DATA BY CONCENTRATION
  #------------------------------------------------------------------
  
  # Filter out control groups - include only Treatment1-5
  treatments_to_include <- c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
  data_filtered <- results$means_adjusted[results$means_adjusted$Treatment %in% treatments_to_include, ]
  
  # Convert treatment to concentration values
  plot_data <- data_filtered
  plot_data$Concentration <- results$fungicide_conc[as.character(plot_data$Treatment)]
  
  # Determine y-axis limits including error bars
  y_min <- min(plot_data[[results$mean_adj_var]] - plot_data[[results$se_adj_var]], na.rm = TRUE)
  y_max <- max(plot_data[[results$mean_adj_var]] + plot_data[[results$se_adj_var]], na.rm = TRUE)
  
  # Set up plotting area with adjusted margins
  par(mar = c(5, 5, 2, 2))  # c(bottom, left, top, right)
  
  # Get min/max concentration values
  min_conc <- min(results$fungicide_conc[treatments_to_include])
  max_conc <- max(results$fungicide_conc[treatments_to_include])
  
  # Get unique days and set up colors
  days <- sort(unique(data_filtered$Day))
  num_days <- length(days)
  
  # Create a gradient of blue colors from light to dark
  day_colors <- colorRampPalette(c("lightblue", "darkblue"))(num_days)
  day_pch <- 15:(15+num_days-1)  # Different point symbols for each day
  
  # Handle log scale appropriately
  # For x-axis (concentration), we always use log scale
  # For y-axis, check if values are appropriate for log scale
  log_param <- "x"  # Always log scale for x-axis (concentration)
  
  # Reset log parameters and multipliers for second plot
  upper_limit_multiplier <- y_multiplier_linear  # Default for linear scale
  
  if(log_scale) {
    # Check if data can use log scale for y-axis (all values must be positive)
    if(y_min <= 0) {
      warning("Cannot use log scale for y-axis with non-positive values. Using linear scale for y-axis.")
    } else {
      log_param <- "xy"  # Log scale for both axes
      # Use smaller multiplier for log scale
      upper_limit_multiplier <- y_multiplier_log
    }
  }
  
  # Create empty plot with improved y-axis limits
  plot(NULL, 
       xlim = c(min_conc * 0.5, max_conc * 1.2),
       ylim = c(max(y_min, 0.001), y_max * upper_limit_multiplier),
       xlab = "Fungicide Concentration (μg/L)", 
       ylab = paste0("Adjusted ", y_label), 
       main = "", # No caption
       log = log_param,  # Log scale for x-axis (and y-axis if log_scale=TRUE)
       xaxt = "n", # Custom x-axis
       cex.lab = 1.2, cex.axis = 1.1)
  
  # Add x-axis with treatment concentrations
  conc_for_axis <- c(0.054, 0.54, 5.4, 54, 540)
  axis_labels <- c("0.054", "0.54", "5.4", "54", "540")
  axis(1, at = conc_for_axis, labels = axis_labels)
  
  # Add grid
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray80")
  
  # Increase line widths to match the thickness in the first plot
  # Using maximum thickness from the first plot as a baseline
  enhanced_line_width <- max(line_widths)
  
  # Plot each day separately
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    
    # Sort by concentration
    day_data <- day_data[order(day_data$Concentration), ]
    if(nrow(day_data) == 0) next
    
    # Plot the line with enhanced thickness based on line spacing factor
    lines(day_data$Concentration, day_data[[results$mean_adj_var]], 
          col = day_colors[i], lwd = enhanced_line_width * line_spacing_factor)
    
    # Plot points
    for (j in 1:nrow(day_data)) {
      # Draw the point
      points(day_data$Concentration[j], day_data[[results$mean_adj_var]][j], 
             col = day_colors[i], pch = day_pch[i], cex = 1.5)
      
      # Add error bars
      arrows(day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] - day_data[[results$se_adj_var]][j],
             day_data$Concentration[j], 
             day_data[[results$mean_adj_var]][j] + day_data[[results$se_adj_var]][j],
             angle = 90, code = 3, length = 0.05, col = day_colors[i], lwd = 1.5)
    }
  }
  
  # Add day labels at the right end of each line where possible
  # Find the rightmost point for each day
  for (i in 1:length(days)) {
    day <- days[i]
    day_data <- plot_data[plot_data$Day == day, ]
    if(nrow(day_data) == 0) next
    
    # Get the point with maximum concentration
    max_conc_idx <- which.max(day_data$Concentration)
    if(length(max_conc_idx) > 0) {
      # Label just to the right of the last point
      text(day_data$Concentration[max_conc_idx] * 1.15, 
           day_data[[results$mean_adj_var]][max_conc_idx], 
           paste("Day", day), 
           col = day_colors[i], cex = 0.8, adj = 0)
    }
  }
  
  # Add legend at top left
  legend("topleft", 
         legend = paste("Day", days), 
         col = day_colors,
         pch = day_pch,
         lty = 1, 
         lwd = 2,
         cex = 0.7, 
         bty = "n")  # no box around legend
  
  # Return color schemes for use in scorecards
  return(list(
    treatment_colors = treatment_colors,
    day_colors = day_colors,
    days = days,
    concentrations = concentrations
  ))
}



# Function to create a complete figure with scorecards and data plots
plot_complete_fungicide_analysis <- function(results,
                                             y_multiplier_linear = 1.1,
                                             y_multiplier_log = 1.05,
                                             asterisk_size = 1.2,
                                             line_spacing_factor = 1.0) {
  # Set up the plotting layout - 2 rows, 2 columns
  # First row: scorecards, Second row: data plots
  layout_matrix <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE)
  layout(layout_matrix)
  
  # First create the data plots and get color information
  old_par <- par(no.readonly = TRUE)  # Save current parameters
  par(mfrow = c(1, 2))  # Temporary 1x2 layout for data plots
  
  # Generate color info from the data plots (without actually plotting)
  color_info <- plot_fungicide_data(results, 
                                    y_multiplier_linear = y_multiplier_linear,
                                    y_multiplier_log = y_multiplier_log,
                                    line_spacing_factor = line_spacing_factor)
  
  # Reset parameters and use the layout
  par(old_par)
  layout(layout_matrix)
  
  # Plot the scorecards in the first row
  plot_significance_scorecards(results, color_info, asterisk_size)
  
  # Plot the data in the second row
  plot_fungicide_data(results, 
                      y_multiplier_linear = y_multiplier_linear,
                      y_multiplier_log = y_multiplier_log,
                      line_spacing_factor = line_spacing_factor)
  
  # Reset the plotting parameters
  par(mfrow = c(1, 1))
}