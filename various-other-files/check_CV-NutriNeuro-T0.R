source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

# set all the experiment and file names
wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)
exp = "NutriNeuro"

datadir = paste0(wd, exp, "/data/")
resultsdir = paste0(wd, exp, "/results/Plasma/")

tlists = paste0(wd, exp, "/Transition-Lists/")
tfile <- "CV-check.csv"

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

# all_data <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "FiS"]))]
all_data <- all_data[, c(which(colnames(all_data) %in% add_info$Sample[add_info$Sample.Type == "Plasma"]))]

exp_data$all_data <- all_data

# Define the molecules to look at from transition list and extract them from the whole data set
transitionFile <- paste0(tlists, tfile)
targeted_experiment_data <- extract_feature_list(all_data, df_ids, transitionFile) %>% select(!contains("QC"))
significant_abundant_features_object <- get_significant_abundant_features(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists, remove_outliers = FALSE)
all_needed_features <- significant_abundant_features_object$significant_abundant_features

full_prep_data <- prepare_data_for_plot(targeted_experiment_data, info_file_dir) %>% filter(Timepoint == 0)
full_prep_all_data <- prepare_data_for_plot(all_needed_features, info_file_dir) %>% filter(Timepoint == 0)


# Calculate CV for targeted molecule
  # for general CV over all samples
    if ("replicate" %in% colnames(full_prep_all_data)) {
      cv_data <- full_prep_data %>%
        group_by(Molecule.Name, replicate) %>%
        summarise(
          Area = mean(Area, na.rm = TRUE),
          rt = min(rt)
        ) %>%
        ungroup() %>%
        group_by(Molecule.Name) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          rt = min(rt)
        ) %>%
        ungroup()
    } else {
      cv_data <- full_prep_data %>%
        group_by(Molecule.Name) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          rt = min(rt)
        ) %>%
        ungroup()
    }
    
    # Plot: Histogram of CVs
    plot <- ggplot(cv_data, aes(x = CV)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of CV of targeted Molecules",
           x = "CV (%)", y = "Count") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-histogram-targeted.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. Mean Area (Scatter)
    plot <- ggplot(cv_data, aes(x = Mean_Area, y = CV, label = Molecule.Name)) +
      geom_point(color = "blue", alpha = 0.6) +
      geom_text_repel(
        size = 3, 
        box.padding = 0.3,      # Smaller padding around text
        point.padding = 0.2,    # Less space between text and point
        segment.color = "gray50",
        segment.size = 0.3,     # Thinner line for better readability
        min.segment.length = 0.1, # Ensures more lines appear
        force = 0.3             # Reduces label spreading
      ) + 
      scale_x_log10() +  # Log scale for better visualization
      labs(title = "CV vs. Mean Area of targeted Molecules",
           x = "Mean Area (log scale)", y = "CV (%)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-Mean-scatter-plot-targeted.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. RT (Scatter)
    plot <- ggplot(cv_data, aes(x = rt, y = CV, label = Molecule.Name)) +
      geom_point(color = "blue", alpha = 0.6) +
      geom_text_repel(
        size = 3, 
        box.padding = 0.3,      # Smaller padding around text
        point.padding = 0.2,    # Less space between text and point
        segment.color = "gray50",
        segment.size = 0.3,     # Thinner line for better readability
        min.segment.length = 0.1, # Ensures more lines appear
        force = 0.3             # Reduces label spreading
      ) + 
      labs(title = "CV vs. RT of targeted Molecules",
           x = "RT (min)", y = "CV (%)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-RT-scatter-plot-targeted.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    
  
# Calculate CV for each molecule
  # for general CV over all samples
    if ("replicate" %in% colnames(full_prep_all_data)) {
      all_cv_data <- full_prep_all_data %>%
        group_by(id, replicate) %>%
        summarise(
          Area = mean(Area, na.rm = TRUE),
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup() %>%
        group_by(id) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup()
    } else {
      all_cv_data <- full_prep_all_data %>%
        group_by(id) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup()
    }
    
    
    # Plot: Histogram of CVs
    plot <- ggplot(all_cv_data, aes(x = CV)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of CV of all Molecules",
           x = "CV (%)", y = "Count") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-histogram-all.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. Mean Area (Scatter)
    plot <- ggplot(all_cv_data, aes(x = Mean_Area, y = CV, label = Annotation)) +
      geom_point(color = "blue", alpha = 0.6) +
      scale_x_log10() +  # Log scale for better visualization
      labs(title = "CV vs. Mean Area of all Molecules",
           x = "Mean Area (log scale)", y = "CV (%)") +
      theme_minimal()
    ggsave(paste0(resultsdir, "CV-vs-Mean-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    # Plot: CV vs. RT (Scatter)
    plot <- ggplot(all_cv_data, aes(x = rt, y = CV, label = Annotation)) +
      geom_point(color = "blue", alpha = 0.6) +
      labs(title = "CV vs. RT of all Molecules",
           x = "RT (min)", y = "CV (%)") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-RT-scatter-plot-all.png"), plot = plot, width = 8, height = 6, dpi = 300)

    
    
    
    
# normalize with Caffeine-D9
    targeted_experiment_data_norm <- targeted_experiment_data
    caffeine_d9_values <- as.numeric(targeted_experiment_data_norm[targeted_experiment_data_norm$Molecule.Name == "Caffeine-D9", 
                                                                   5:(ncol(targeted_experiment_data) - 1)])
    targeted_experiment_data_norm[,5:(ncol(targeted_experiment_data)-1)] <- sweep(as.matrix(targeted_experiment_data_norm[5:(ncol(targeted_experiment_data) - 1)]), 2, caffeine_d9_values, "/")
    all_needed_features_norm <- all_needed_features
    all_needed_features_norm[10:ncol(all_needed_features_norm)]  <- sweep(as.matrix(all_needed_features_norm[10:ncol(all_needed_features_norm)]), 2, caffeine_d9_values, "/")
    
    full_prep_data_norm <- prepare_data_for_plot(targeted_experiment_data_norm, info_file_dir)
    full_prep_all_data_norm <- prepare_data_for_plot(all_needed_features_norm, info_file_dir)
    
    
    # Calculate CV for targeted molecule
    # for general CV over all samples
    if ("replicate" %in% colnames(full_prep_data)) {
      cv_data_norm <- full_prep_data_norm %>%
        group_by(Molecule.Name, replicate) %>%
        summarise(
          Area = mean(Area, na.rm = TRUE),
          rt = min(rt)
        ) %>%
        ungroup() %>%
        group_by(Molecule.Name) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          rt = min(rt)
        ) %>%
        ungroup()
    } else {
      cv_data_norm <- full_prep_data_norm %>%
        group_by(Molecule.Name) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          rt = min(rt)
        ) %>%
        ungroup()
    }
    
    cv_data_norm$CV_raw <- cv_data$CV
    cv_data_norm$diff <- cv_data_norm$CV_raw - cv_data_norm$CV
    cv_data_norm$color <- ifelse(cv_data_norm$diff < 0, "red", "blue")
    
    # Plot: Histogram of CVs
    plot <- ggplot(cv_data_norm, aes(x = CV)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of CV of normalized targeted Molecules",
           x = "CV (%)", y = "Count") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-histogram-targeted-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. Mean Area (Scatter)
    plot <- ggplot(cv_data_norm, aes(x = Mean_Area, y = CV, label = Molecule.Name, color = color)) +
      geom_point(alpha = 0.6) +
      geom_text_repel(
        size = 3, 
        box.padding = 0.3,      # Smaller padding around text
        point.padding = 0.2,    # Less space between text and point
        segment.color = "gray50",
        segment.size = 0.3,     # Thinner line for better readability
        min.segment.length = 0.1, # Ensures more lines appear
        force = 0.3             # Reduces label spreading
      ) + 
      scale_x_log10() +  # Log scale for better visualization
      labs(title = "CV vs. Mean Area of normalized targeted Molecules",
           x = "Mean Area (log scale)", y = "CV (%)") +
      scale_color_identity() +  # Use the colors from the 'color' column
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-Mean-scatter-plot-targeted-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. RT (Scatter)
    plot <- ggplot(cv_data_norm, aes(x = rt, y = CV, label = Molecule.Name, color = color)) +
      geom_point(alpha = 0.6) +
      geom_text_repel(
        size = 3, 
        box.padding = 0.3,      # Smaller padding around text
        point.padding = 0.2,    # Less space between text and point
        segment.color = "gray50",
        segment.size = 0.3,     # Thinner line for better readability
        min.segment.length = 0.1, # Ensures more lines appear
        force = 0.3             # Reduces label spreading
      ) + 
      labs(title = "CV vs. RT of normalized targeted Molecules",
           x = "RT (min)", y = "CV (%)") +
      scale_color_identity() +  # Use the colors from the 'color' column
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-RT-scatter-plot-targeted-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    
    
    # Calculate CV for each molecule
    # for general CV over all samples
    if ("replicate" %in% colnames(full_prep_all_data)) {
      all_cv_data_norm <- full_prep_all_data_norm %>%
        group_by(id, replicate) %>%
        summarise(
          Area = mean(Area, na.rm = TRUE),
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup() %>%
        group_by(id) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup()
    } else {
      all_cv_data_norm <- full_prep_all_data_norm %>%
        group_by(id) %>%
        summarise(
          Mean_Area = mean(Area, na.rm = TRUE),
          SD_Area = sd(Area, na.rm = TRUE),
          CV = (SD_Area / Mean_Area) * 100,
          Annotation = min(Annotation),
          rt = min(rt)
        ) %>%
        ungroup()
    }
    
    all_cv_data_norm$CV_raw <- all_cv_data$CV
    all_cv_data_norm$diff <- all_cv_data_norm$CV_raw - all_cv_data_norm$CV
    all_cv_data_norm$color <- ifelse(all_cv_data_norm$diff < 0, "red", "blue")
    
    # Plot: Histogram of CVs
    plot <- ggplot(all_cv_data_norm, aes(x = CV)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = "Distribution of CV of all normalized Molecules",
           x = "CV (%)", y = "Count") +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-histogram-all-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    # Plot: CV vs. Mean Area (Scatter)
    plot <- ggplot(all_cv_data_norm, aes(x = Mean_Area, y = CV, label = Annotation, color = color)) +
      geom_point(alpha = 0.6) +
      scale_x_log10() +  # Log scale for better visualization
      labs(title = "CV vs. Mean Area of all normalized Molecules",
           x = "Mean Area (log scale)", y = "CV (%)") +
      scale_color_identity() +  # Use the colors from the 'color' column
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-Mean-scatter-plot-all-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    # Plot: CV vs. RT (Scatter)
    plot <- ggplot(all_cv_data_norm, aes(x = rt, y = CV, label = Annotation, color = color)) +
      geom_point(alpha = 0.6) +
      labs(title = "CV vs. RT of all Molecules",
           x = "RT (min)", y = "CV (%)") +
      scale_color_identity() +  # Use the colors from the 'color' column
      theme_minimal() +
      theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14), 
        plot.title = element_text(size = 16, face = "bold") 
      )
    ggsave(paste0(resultsdir, "CV-vs-RT-scatter-plot-all-normalized.png"), plot = plot, width = 8, height = 6, dpi = 300)
    
    
    
    
    