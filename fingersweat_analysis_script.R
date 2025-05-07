source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp = "NutriNeuro"

datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/")

tlists = paste0(wd, exp, "/Transition-Lists/")
#tfile <- "biological_normalization_molecule_list_updated-version.csv"
tfile <- "Screening_TransitionList.csv"
#tfile <- "potentialy-regulated-molecules.csv"

info_file <- "/Metadata.csv"
info_file_dir <- paste0(wd, exp, info_file)
add_info <- read.csv(info_file_dir)
if (length(add_info) < 2) {
  add_info <- read.csv2(info_file_dir)
}

#normalization_list  = paste0(wd, "normalization/biological_normalization_molecule_list_updated-version.csv")

###############################
## Start of Data Processing  ##
###############################

# Get the whole processed data set from file (already prepared csv file with mzMine)
  exp_data <- process_merged_files(datadir)
  all_data <- exp_data$data
  df_ids <- exp_data$df_ids



# Look at finger sweat data
  resultsdir_FiS = paste0(wd, exp, "/results/FiS/")
  if (!dir.exists(resultsdir_FiS)) {
    dir.create(resultsdir_FiS, recursive = TRUE)
  }
  
  all_data_FiS <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "FiS"]))]
  exp_data_FiS <- exp_data
  exp_data_FiS$data <- all_data_FiS
  
  # Define the molecules to look at from transition list and extract them from the whole data set
  transitionFile <- paste0(tlists, tfile)
  targeted_experiment_data_FiS <- extract_feature_list_advanced(all_data_FiS, df_ids, transitionFile)
  significant_abundant_features_object_FiS <- get_significant_abundant_features(datadir, exp_data_FiS, info_file_dir, resultsdir_FiS, figuredir, exp, tlists)
  all_needed_features_FiS <- significant_abundant_features_object_FiS$significant_abundant_features
  
  
    # Get some overview of the data
      Group_list <- list("Donor", "Intervention", "Timepoint", "Age", "Sex", "BMI", "rel_fatMass", "abs_fatMass", "abs_fatFreeMass",
                         "abs_muscleMass", "waist_c", "weight_in_kg", "height_in_m", "FFMI", "FMI", "Batch.Nr")
      Groups_to_group <- list("Age", "BMI", "rel_fatMass", "abs_fatMass", "abs_fatFreeMass",
                              "abs_muscleMass", "waist_c", "weight_in_kg", "height_in_m", "FFMI", "FMI")
      plot_to_pdf_objective_data_overview(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, Groups_to_group)
      plot_to_pdf_objective_data_overview(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, Groups_to_group, color_donor_only=TRUE)
      
      
      
    # Compare Tyrosine and Tryptophan to LNAAs
      Group_list <- list("Intervention", "Sex")
      plot_to_pdf_Tyr_to_Trp(resultsdir_FiS, targeted_experiment_data_FiS, add_info, Group_list)
      plot_to_pdf_Tyr_Trp_to_LNAA(resultsdir_FiS, targeted_experiment_data_FiS, add_info, Group_list)
      
    # Compare different settings with limma and plot as volcano
      Group_list <- list("Intervention", "Sex")
      ann_feature_list <- list('Nuclease, ribo-', 'gamma-Glutamylglutamine', 'Indole-3-acetate', 'Xylitol', '1-Propene-1,2,3-tricarboxylic acid',
                               'Betaine', 'Orotic Acid', '6-o-Phosphonohexopyranose', 'Hexahydrohippurate', '(3R)-3,9-dihydroxynonanoic acid', 
                               'Glycerin', 'Citrulline, (+/-)-', 'Suberylglycine', 'Glucoheptonic Acid', 'Leu-Phe-Ala', 'Tyr-Ser-Ala', 
                               '11-Cyanoundecanoic acid', 'L-leucyl-L-glutaminyl-L-proline', 'Hexitol', 'N-(p-nitrophenyl)-beta-alaninamide',
                               'DL-alanine', 'Sebacic Acid', 'Nonaethylene glycol', '	Saccharin', 'Phenylalanylvalyllysine', 'Ile-Pro-Phe',
                               'Phenylalanylproline', 'Gln-Tyr')
      
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.7, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.3, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.3, ann_features = ann_feature_list)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = FALSE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = FALSE, adj_p_val = TRUE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_FiS, all_needed_features_FiS, add_info, Group_list, only_annotated_features = FALSE, adj_p_val = FALSE, nr_ann_features = 10)
      
      
      plot_significant_changed_features(all_needed_features_FiS, 'limma_all-resulttables/only-annotated-features_adj-p-value_confidence-cutoff0.3/', resultsdir_FiS, info_file_dir)
      
    # Check for clusters of the features
      fold_change_df <- get_fold_change_dataframe(all_needed_features_FiS, add_info, Group_list = list('Intervention', 'Sex'))
      cluster_data <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df)
      full_prep_data <- prepare_data_for_plot(all_needed_features_FiS, info_file_dir)
      annotated_clustered_features <- as.data.frame(cluster_data[[3]]) %>% rownames_to_column("id") %>%
        left_join(full_prep_data %>% select(id, Annotation, rt, mz) %>% mutate(id = as.character(id)), by = "id") %>%
        mutate(Annotation = if_else(is.na(Annotation), paste0(rt, "@", mz), Annotation)) %>% distinct(id, .keep_all = TRUE)
      
      Group_list = list('Intervention', 'Sex')
      nr_comb <- nrow(expand.grid(lapply(Group_list, function(group) unique(full_prep_data[[group]])), KEEP.OUT.ATTRS = FALSE))
      
      pdf(paste0(resultsdir_FiS,"cluster_hierarchical_automatic_plots.pdf"), width = 25, height = 12)
      plots <- plot_clusters(cluster_data, full_prep_data, Group_list = list('Intervention', 'Sex')) # Generate plots for the current clustering
      grid.arrange(grobs = plots, ncol = nr_comb, nrow = length(plots)/nr_comb)
      dev.off()

      
      
# Look at EDTA blood plasma data
  resultsdir_Plasma = paste0(wd, exp, "/results/Plasma/")
  if (!dir.exists(resultsdir_Plasma)) {
    dir.create(resultsdir_Plasma, recursive = TRUE)
  }
  
  all_data_Plasma <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "Plasma"]))]
  exp_data_Plasma <- exp_data
  exp_data_Plasma$data <- all_data_Plasma
      
  # Define the molecules to look at from transition list and extract them from the whole data set
  transitionFile <- paste0(tlists, tfile)
  targeted_experiment_data_Plasma <- extract_feature_list_advanced(all_data_Plasma, df_ids, transitionFile)
  significant_abundant_features_object_Plasma <- get_significant_abundant_features(datadir, exp_data_Plasma, info_file_dir, resultsdir_Plasma, figuredir, exp, tlists)
  all_needed_features_Plasma <- significant_abundant_features_object_Plasma$significant_abundant_features
  
  add_info$Column.Batch <- paste0(add_info$Column, '_', add_info$Batch.Nr)
  
  # Batch correction for Column and Batch.Nr
    pure_data_plasma <- as.matrix(all_needed_features_Plasma[, grep("^NutriNeuro_Plasma", colnames(all_needed_features_Plasma))])
    rownames(pure_data_plasma) <- all_needed_features_Plasma$id
    meta_aligned <- add_info[match(gsub("\\.", "_", colnames(pure_data_plasma)), add_info$Sample), ]
    batch <- factor(paste0("Col", meta_aligned$Column, "_Batch", meta_aligned$Batch.Nr))

    sel <- complete.cases(meta_aligned[, c("Column", "Batch.Nr", c("Age", "Sex", "Intervention", "BMI", "rel_fatMass", "FFMI", "FMI"))])
    pure_data_plasma_filt <- as.matrix(pure_data_plasma[, sel])      # samples (columns)
    meta_aligned_filt     <- meta_aligned[sel, ]          # samples (rows)
    batch                 <- factor(paste0("Col", meta_aligned_filt$Column, "_Batch", meta_aligned_filt$Batch.Nr))
    mod                   <- model.matrix(~ Age + Sex + Intervention + BMI + rel_fatMass + FFMI + FMI, data=meta_aligned_filt)

    # Batch correction (preserve other covariates with mod)
    mod <- model.matrix(~ Age + Sex + Intervention + BMI + rel_fatMass + FFMI + FMI, data=meta_aligned_filt)
    corrected_expr_matrix <- ComBat(dat=pure_data_plasma_filt, batch=batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)

    all_needed_features_Plasma <- cbind(all_needed_features_Plasma[1:9], as.data.frame(corrected_expr_matrix))
  
    # Get some overview of the data
      Group_list <- list("Donor", "Intervention", "Timepoint", "Age", "Sex", "BMI", "rel_fatMass", "abs_fatMass", "abs_fatFreeMass",
                         "abs_muscleMass", "weight_in_kg", "height_in_m", "FFMI", "FMI", "Column", "Batch.Nr", "Column.Batch")
      Groups_to_group <- list("Age", "BMI", "rel_fatMass", "abs_fatMass", "abs_fatFreeMass",
                              "abs_muscleMass", "weight_in_kg", "height_in_m", "FFMI", "FMI")
      plot_to_pdf_objective_data_overview(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, Groups_to_group)
      plot_to_pdf_objective_data_overview(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, Groups_to_group, color_donor_only=TRUE)
      
      
      
    # Compare Tyrosine and Tryptophan to LNAAs
      Group_list <- list("Intervention", "Sex")
      plot_to_pdf_Tyr_to_Trp(resultsdir_Plasma, targeted_experiment_data_Plasma, add_info, Group_list)
      plot_to_pdf_Tyr_Trp_to_LNAA(resultsdir_Plasma, targeted_experiment_data_Plasma, add_info, Group_list)
      
    # Compare different settings with limma and plot as volcano
      Group_list <- list("Intervention", "Sex")
      
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.7, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.3, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.3, nr_ann_features = 0)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = FALSE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = FALSE, adj_p_val = TRUE, nr_ann_features = 10)
      plot_to_pdf_limma_test_results(resultsdir_Plasma, all_needed_features_Plasma, add_info, Group_list, only_annotated_features = FALSE, adj_p_val = FALSE, nr_ann_features = 10)
      
      plot_significant_changed_features(all_needed_features_Plasma, 'limma_all-resulttables/only-annotated-features_adj-p-value_confidence-cutoff0.3/', resultsdir_Plasma, info_file_dir)
      
    # Check for clusters of the features
      fold_change_df <- get_fold_change_dataframe(all_needed_features_Plasma, add_info, Group_list = list('Intervention', 'Sex'))
      cluster_data <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df)
      full_prep_data <- prepare_data_for_plot(all_needed_features_Plasma, info_file_dir)
      annotated_clustered_features <- as.data.frame(cluster_data[[3]]) %>% rownames_to_column("id") %>%
        left_join(full_prep_data %>% select(id, Annotation, rt, mz) %>% mutate(id = as.character(id)), by = "id") %>%
        mutate(Annotation = if_else(is.na(Annotation), paste0(rt, "@", mz), Annotation)) %>% distinct(id, .keep_all = TRUE)
      
      Group_list = list('Intervention', 'Sex')
      nr_comb <- nrow(expand.grid(lapply(Group_list, function(group) unique(full_prep_data[[group]])), KEEP.OUT.ATTRS = FALSE))
      
      pdf(paste0(resultsdir_Plasma,"cluster_hierarchical_automatic_plots.pdf"), width = 25, height = 12)
      plots <- plot_clusters(cluster_data, full_prep_data, Group_list = list('Intervention', 'Sex')) # Generate plots for the current clustering
      grid.arrange(grobs = plots, ncol = nr_comb, nrow = length(plots)/nr_comb)
      dev.off()
      
      
      
      