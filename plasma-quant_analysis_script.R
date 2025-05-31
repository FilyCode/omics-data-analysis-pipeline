source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp = "NutriNeuro/Plasma_quant"

datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/")

tlists = paste0(wd, exp, "/Transition-Lists/")
tfile <- "quant-molecules.csv"

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


  # Define the molecules to look at from transition list and extract them from the whole data set
  transitionFile <- paste0(tlists, tfile)
  targeted_experiment_data_Plasma <- extract_feature_list(all_data, df_ids, transitionFile) %>% select(-QC)
  
  plasma_standards <- targeted_experiment_data_Plasma %>%
    filter(grepl("heavy", Molecule.Name, ignore.case = TRUE)) %>%
    select(matches("^NutriNeuro_Plasma_S"), Molecule.Name)
      
  regression_fits <- create_calibration_curves(plasma_standards, add_info, resultsdir)
  show_light_heavy_ratio(targeted_experiment_data_Plasma, add_info, resultsdir, nr_molecules = 22)
  aminoacid_overview(targeted_experiment_data_Plasma, resultsdir, nr_molecules = 22)
  quantitive_AA <- get_quantitative_aminoacid_concentrations(targeted_experiment_data_Plasma, add_info, resultsdir, nr_molecules, output_file = "quantitative-aminoacid-concentrations.csv")
    
      
  # Compare Tyrosine and Tryptophan to LNAAs
    Group_list <- list("Intervention", "Sex")
    plot_to_pdf_Tyr_to_Trp(resultsdir, quantitive_AA, add_info, Group_list)
    plot_to_pdf_Tyr_Trp_to_LNAA(resultsdir, quantitive_AA, add_info, Group_list)
      
    