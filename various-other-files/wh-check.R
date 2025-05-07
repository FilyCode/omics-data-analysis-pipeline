source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp = "NutriNeuro"

datadir = paste0(wd, exp, "/data_wh-check/")
resultsdir = paste0(wd, exp, "/data_wh-check/")

###############################
## Start of Data Processing  ##
###############################

# Get the whole processed data set from file (already prepared csv file with mzMine)
exp_data <- process_merged_files(datadir)
all_data <- exp_data$data
df_ids <- exp_data$df_ids


non_wh_cols <- colnames(all_data)[!grepl("_WH$", colnames(all_data))]
result_list <- list()

# Loop through each non-WH column and divide by its WH pair
for (col in non_wh_cols) {
  wh_col <- paste0(col, "_WH")
  if (wh_col %in% colnames(all_data)) {
    result_list[[col]] <- all_data[[col]] / all_data[[wh_col]]
  } else {
    warning(paste("WH pair not found for column:", col))
  }
}

# Combine into a new dataframe
result_df <- as.data.frame(result_list)
result_df_combined <- cbind(df_ids, result_df)
result_df_combined$CV <- apply(result_df, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
result_df_combined$mean_ratio <- rowMeans(result_df, na.rm = TRUE)


nutrineuro_long <- melt(result_df, variable.name = "Sample", value.name = "Ratio")

ggplot(nutrineuro_long, aes(x = Sample, y = Ratio)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of NutriNeuro Ratios", x = "Sample", y = "Ratio")

ggplot(result_df_combined, aes(x = rt, y = mean_ratio)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Mean NutriNeuro Ratio over Retention Time",
       x = "Retention Time (rt)",
       y = "Mean Ratio")

tlists = paste0(wd, exp, "/Transition-Lists/")
tfile <- "Screening_TransitionList.csv"
transitionFile <- paste0(tlists, tfile)
targeted_experiment_data_FiS <- extract_feature_list_advanced(all_data, df_ids, transitionFile) %>% select(-c(1:4,"Molecule.Name"))

non_wh_cols <- colnames(targeted_experiment_data_FiS)[!grepl("_WH$", colnames(targeted_experiment_data_FiS))]
result_list <- list()

# Loop through each non-WH column and divide by its WH pair
for (col in non_wh_cols) {
  wh_col <- paste0(col, "_WH")
  if (wh_col %in% colnames(targeted_experiment_data_FiS)) {
    result_list[[col]] <- targeted_experiment_data_FiS[[col]] / targeted_experiment_data_FiS[[wh_col]]
  } else {
    warning(paste("WH pair not found for column:", col))
  }
}

# Combine into a new dataframe
result_df_targeted <- as.data.frame(result_list)

nutrineuro_long_targeted <- melt(result_df_targeted, variable.name = "Sample", value.name = "Ratio")

ggplot(nutrineuro_long_targeted, aes(x = Sample, y = Ratio)) +
  geom_boxplot() +
  scale_y_log10() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot of NutriNeuro Ratios", x = "Sample", y = "Ratio")

