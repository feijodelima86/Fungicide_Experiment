library(readr)

# Load your data
data <- read.csv("~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM.csv")

# Define old and new column names
old_names <- c("Sample", "Group", "Date", "Day", "Fungi_Est", "Acet_Est", "AI", "N/P", "C/P ratio", "P/D", "SRP")
new_names <- c("sample_id", "treatment_group", "sample_date", "experiment_day", "fungi_concentration", 
               "acetone_concentration", "autotrophic_index", "np_ratio", "cp_ratio", "pd_ratio", 
               "soluble_reactive_phosphorus")

# Find which columns need to be renamed
cols_to_rename <- which(colnames(data) %in% old_names)

# Rename only the matching columns
colnames(data)[cols_to_rename] <- new_names[match(colnames(data)[cols_to_rename], old_names)]

# Transform the autotrophic index
data$autotrophic_index <- log10(data$autotrophic_index + 1)

# Save the renamed data if needed
write.csv(data, "~/GIT_Projects/Fungicide_Experiment/0_Data/Data_Merged_OLRM_renamed.csv", row.names = FALSE)