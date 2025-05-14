# Load the data

data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

names(data)

# Run the full analysis and store results

results <- analyze_fungicide_effect(data, "afdm", log_transform = TRUE)

# Create the plot using the create_plot function

create_plot(
  data = results$means_original,
  y_var = "mean_afdm", 
  se_var = "se_afdm", 
  y_label = "afdm (Autotrophic Index)", 
  mafdmn_title = "Original AFDM Values Over Time"
)


create_plot(
  data = results$means_adjusted,
  y_var = "mean_adjusted_afdm", 
  se_var = "se_adjusted_afdm", 
  y_label = "Adjusted afdm (Autotrophic Index)", 
  mafdmn_title = "Acetone-Adjusted afdm Values Over Time"
)

create_selective_plot(
  data = results$means_original,
  y_var = "mean_afdm", 
  se_var = "se_afdm", 
  y_label = "afdm (Autotrophic Index)", 
  mafdmn_title = "Controls and Highest Treatment", 
  selected_treatments = c("Control1", "Control2", "Treatment5")
)

create_selective_plot(
  data = results$means_adjusted,
  y_var = "mean_adjusted_afdm", 
  se_var = "se_adjusted_afdm", 
  y_label = "Adjusted afdm (Autotrophic Index)", 
  mafdmn_title = "All Fungicide Treatments (Acetone-Adjusted)", 
  selected_treatments = c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")
)

# Select specific treatments to compare
selected_treatments = c("Treatment1", "Treatment2", "Treatment3", "Treatment4", "Treatment5")

# Create plot showing significant differences between selected treatments
create_selective_plot_sig_days(
  data = results$means_adjusted,
  y_var = results$mean_adj_var, 
  se_var = results$se_adj_var, 
  y_label = "Adjusted Response", 
  mafdmn_title = "Treatment Comparison with Significance", 
  selected_treatments = selected_treatments,
  significant_points = results$significant_points,
  selected_day = NULL,  # Set to NULL to show all days, or a specific day number (e.g., 7) to focus on one day
  log_scale = T
)







# Define the fungicide concentrations

fungicide_conc <- c(
  "Control1" = 0,      # No fungicide
  "Control2" = 0,      # No fungicide, only acetone
  "Treatment1" = 0.054, # 0.054 μg/L
  "Treatment2" = 0.54,  # 0.54 μg/L
  "Treatment3" = 5.4,   # 5.4 μg/L
  "Treatment4" = 54,    # 54 μg/L
  "Treatment5" = 540    # 540 μg/L
)

create_concentration_plot(
  data = results$means_original,
  y_var = "mean_afdm", 
  se_var = "se_afdm", 
  y_label = "afdm (Autotrophic Index)", 
  mafdmn_title = "Original afdm Data by Fungicide Concentration",
  fungicide_conc = fungicide_conc,
  significant_points = NULL  # No significance testing for original data
)

create_concentration_plot(
  data = results$means_adjusted,
  y_var = "mean_adjusted_afdm", 
  se_var = "se_adjusted_afdm", 
  y_label = "Adjusted afdm (Autotrophic Index)", 
  mafdmn_title = "Acetone-Adjusted afdm Data by Fungicide Concentration\nAsterisks (*) indicate significant differences from Control2",
  fungicide_conc = fungicide_conc,
  significant_points = results$significant_points
)


