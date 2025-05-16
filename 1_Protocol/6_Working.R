# Load the data

data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

names(data)

# Run the full analysis and store results

results <- analyze_fungicide_effect(data, "Total.P", log_transform = F)

plot_selective_treatment_comparison_alt(results)

plot_fungicide_data(results)

plot_significance_scorecards(results)

plot_fungicide_analysis(results)

names(data)

