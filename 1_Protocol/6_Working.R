# Load the data

data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

names(data)

# Run the full analysis and store results

results <- analyze_fungicide_effect(data, "AI", log_transform = T)

plot_treatment_comparison(results)

plot_selective_treatment_comparison(results)

plot_concentration_response(results)

plot_jittered_treatment_comparison(results)

plot_jittered_concentration_response(results)

plot_selective_treatment_comparison_alt(results)

