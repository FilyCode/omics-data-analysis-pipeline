source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp = "NutriNeuro"

datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/")

info_file <- "/Metadata.csv"
info_file_dir <- paste0(wd, exp, info_file)
add_info <- read.csv(info_file_dir)
if (length(add_info) < 2) {
  add_info <- read.csv2(info_file_dir)
}

###############################
## Start of Data Processing  ##
###############################

# Get the whole processed data set from file (already prepared csv file with mzMine)
exp_data <- process_merged_files(datadir)
all_data <- exp_data$data
df_ids <- exp_data$df_ids

# Look at finger sweat data
resultsdir_FiS = paste0(wd, exp, "/results/FiS/")
all_data_FiS <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "FiS"]))]
exp_data_FiS <- exp_data
exp_data_FiS$data <- all_data_FiS

# Define the molecules to look at from transition list and extract them from the whole data set
significant_abundant_features_object_FiS <- get_significant_abundant_features(datadir, exp_data_FiS, info_file_dir, resultsdir_FiS, figuredir, exp, tlists)
all_needed_features_FiS <- significant_abundant_features_object_FiS$significant_abundant_features

plot_feature_ranking_overview_to_pdf(all_needed_features_FiS, resultsdir_FiS)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_FiS, resultsdir_FiS, add_info, "Timepoint", 0)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_FiS, resultsdir_FiS, add_info, "Timepoint", 1)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_FiS, resultsdir_FiS, add_info, "Sex", "Male")
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_FiS, resultsdir_FiS, add_info, "Sex", "Female")



# Look at blood plasma data
resultsdir_Plasma = paste0(wd, exp, "/results/Plasma/")
all_data_Plasma <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "Plasma"]))]
exp_data_Plasma <- exp_data
exp_data_Plasma$data <- all_data_Plasma

# Define the molecules to look at from transition list and extract them from the whole data set
significant_abundant_features_object_Plasma <- get_significant_abundant_features(datadir, exp_data_Plasma, info_file_dir, resultsdir_Plasma, figuredir, exp, tlists)
all_needed_features_Plasma <- significant_abundant_features_object_Plasma$significant_abundant_features

plot_feature_ranking_overview_to_pdf(all_needed_features_Plasma, resultsdir_Plasma)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_Plasma, resultsdir_Plasma, add_info, "Timepoint", 0)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_Plasma, resultsdir_Plasma, add_info, "Timepoint", 3)
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_Plasma, resultsdir_Plasma, add_info, "Sex", "Male")
plot_feature_ranking_overview_sample_filtering_to_pdf(all_needed_features_Plasma, resultsdir_Plasma, add_info, "Sex", "Female")


