source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp_quant = "NutriNeuro/Plasma_quant"

datadir_quant = paste0(wd, exp_quant, "/data/")
resultsdir_quant = paste0(wd, exp_quant, "/results/")

tlists = paste0(wd, exp_quant, "/Transition-Lists/")
tfile <- "quant-molecules.csv"

info_file <- "/Metadata.csv"
info_file_dir <- paste0(wd, exp_quant, info_file)
add_info_quant <- read.csv(info_file_dir)
if (length(add_info) < 2) {
  add_info_quant <- read.csv2(info_file_dir)
}

###############################
## Start of Data Processing  ##
###############################

# Get the whole processed data set from file (already prepared csv file with mzMine)
  exp_data_quant <- process_merged_files(datadir)
  all_data_quant <- exp_data_quant$data
  df_ids_quant <- exp_data_quant$df_ids


  # Define the molecules to look at from transition list and extract them from the whole data set
  transitionFile <- paste0(tlists, tfile)
  targeted_experiment_data_Plasma_quant <- extract_feature_list(all_data_quant, df_ids_quant, transitionFile) %>% select(-QC)
  
  targeted_experiment_data_Plasma_quant <- targeted_experiment_data_Plasma_quant[-grep( "Gregor", colnames(targeted_experiment_data_Plasma_quant))] # remove other sample not from NutriNeuro
  
  plasma_standards <- targeted_experiment_data_Plasma_quant %>%
    filter(grepl("heavy", Molecule.Name, ignore.case = TRUE)) %>%
    select(matches("^NutriNeuro_Plasma_S"), Molecule.Name)
      
  regression_fits <- create_calibration_curves(plasma_standards, add_info_quant, resultsdir_quant)
  show_light_heavy_ratio(targeted_experiment_data_Plasma_quant, add_info_quant, resultsdir_quant, nr_molecules = 22)
  aminoacid_overview(targeted_experiment_data_Plasma_quant, resultsdir_quant, nr_molecules = 22)
  quantitive_AA <- get_quantitative_aminoacid_concentrations(targeted_experiment_data_Plasma_quant, add_info_quant, resultsdir_quant, nr_molecules, output_file = "quantitative-aminoacid-concentrations.csv")
    
  # Compare Tyrosine and Tryptophan to LNAAs
    Group_list <- list("Intervention", "Sex")
    plot_to_pdf_Tyr_Trp_to_LNAA_quant(resultsdir_quant, quantitive_AA, add_info_quant, Group_list)
      
    