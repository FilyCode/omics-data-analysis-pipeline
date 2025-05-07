####################################################################################
####################################################################################
### This is a R file that contains all needed functions for checking normalization  
### strategies and evaluate what works best, to be used in other scripts for
### easier readability and handling
####################################################################################
####################################################################################



######################################
## Assisting Functions for plotting ##
######################################

# Function to calculate RSD (Relative Standard Deviation)
calculate_rsd <- function(data, groups, Group) {
  rsd_data <- data.frame(Feature = character(),
                         RSD = numeric(),
                         Group = character())
  for(group in groups) {
    group_df <- subset(data, data[[Group]] == group, select = -c(1, (length(data)-5):length(data)))
    
    rsd_values <- apply(group_df, 2, function(x) {
      return(sd(x) / mean(x) * 100)
    })
    
    rsd_subdata <- data.frame(Feature = colnames(group_df),
                              RSD = rsd_values,
                              Group = group)
    rsd_data <- rbind(rsd_data, rsd_subdata)
  }
  
  return(rsd_data)
}



# function to normalize data with different molecules
normalized_data_with_single_molecules <- function(norm_feature_list, merged_data, pure_data, sample_names, add_info) {
  normalized_data_list <- list()
  rownames(pure_data) <- sample_names
  
  for (norm_feature in norm_feature_list$id) {
    norm_molecule_name <- norm_feature_list$Molecule.Name[norm_feature_list$id == norm_feature]
    norm_merged_data <- merged_data
    norm_pure_data <- pure_data
    
    # Pre-compute normalization factors
    norm_factors <- norm_feature_list[norm_feature_list$id == norm_feature, 5:(4 + length(sample_names))]
    # Apply normalization to merged data
    norm_merged_data[, 2:(length(norm_merged_data) - (length(add_info) -1) )] <- sweep(norm_merged_data[, 2:(length(norm_merged_data) - (length(add_info) -1))], 1, as.matrix(norm_factors), "/")
    # Apply normalization to pure data
    norm_pure_data <- sweep(norm_pure_data, 1, as.matrix(norm_factors), "/")
    # Remove the normalization feature columns from data
    norm_merged_data <- norm_merged_data[, !colnames(norm_merged_data) %in% norm_feature]
    norm_pure_data <- norm_pure_data[, !colnames(norm_pure_data) %in% norm_feature]
    
    # Exchange NA and Inf, as it makes problems later on
    norm_merged_data[, 2:(length(norm_merged_data) - 6)] <- as.data.frame(apply(norm_merged_data[, 2:(length(norm_merged_data) - 6)], 2, replace_na_and_Inf))
    norm_pure_data <- as.data.frame(apply(norm_pure_data, 2, replace_na_and_Inf))
    rownames(norm_pure_data) <- rownames(pure_data)
    
    # Store the normalized dataframes in the list
    normalized_data_list[[norm_feature]] <- list(merged_data = norm_merged_data, pure_data = norm_pure_data, molecule_name = norm_molecule_name)
  }
  
  return(normalized_data_list)
}



# function to normalize data with different normalization methods
normalized_data_with_various_normalization_methods <- function(norm_feature_list, merged_data, pure_data, norm_methods, sample_names, add_info) {
  normalized_data_list <- list()
  rownames(pure_data) <- sample_names
  
  for (method in norm_methods) {
    norm_factors <- calculate_norm_factors(t(pure_data), method, norm_feature_list$Molecule.Name, norm_feature_list)
    norm_factors <- t(as.data.frame(norm_factors))
    
    norm_method_name <- method
    norm_merged_data <- merged_data
    norm_pure_data <- pure_data
    
    # Apply normalization to merged data
    norm_merged_data[, 2:(length(norm_merged_data) - (length(add_info) -1) )] <- sweep(norm_merged_data[, 2:(length(norm_merged_data) - (length(add_info) -1))], 1, as.matrix(norm_factors), "/")
    # Apply normalization to pure data
    norm_pure_data <- sweep(norm_pure_data, 1, as.matrix(norm_factors), "/")
    
    # Remove the normalization feature columns from data (is needed if only one molecule like Tyrosine was used)
    norm_merged_data <- norm_merged_data[, !sapply(norm_merged_data, function(col) all(col == 1))]
    norm_pure_data <- norm_pure_data[, !sapply(norm_pure_data, function(col) all(col == 1))]
    
    # Exchange NA and Inf, as it makes problems later on
    norm_merged_data[, 2:(length(norm_merged_data) - 6)] <- as.data.frame(apply(norm_merged_data[, 2:(length(norm_merged_data) - 6)], 2, replace_na_and_Inf))
    norm_pure_data <- as.data.frame(apply(norm_pure_data, 2, replace_na_and_Inf))
    rownames(norm_pure_data) <- rownames(pure_data)
    
    # Store the normalized dataframes in the list
    normalized_data_list[[method]] <- list(merged_data = norm_merged_data, pure_data = norm_pure_data, molecule_name = norm_method_name)
  }
  
  return(normalized_data_list)
}



# Function to calculate normalization factors for different methods
calculate_norm_factors <- function(data, method, normalization_molecules = NULL, targeted_norm_molecules = NULL) {
  if ("Confidence" %in% colnames(data)) {
    start_row <- 9
  } else{
    start_row <- 0
  }
  
  norm_factors <- numeric(ncol(data) - start_row)
  
  for (i in (start_row+1):ncol(data)) {
    sample_name <- colnames(data)[i]
    
    if (method == "Tyrosine") {
      if("Annotation" %in% colnames(data)) {
        tyrosine_data <- data[data$Annotation %in% "DL-Tyrosine", i, drop = FALSE]
      } else {
        tyrosine_data <- data[targeted_norm_molecules$Molecule.Name %in% "DL-Tyrosine", i, drop = FALSE]
      }
      
      norm_factors[i - start_row] <- max(tyrosine_data, na.rm = TRUE)
      
    } else if (method == "Total_Sum") {
      norm_factors[i - start_row] <- sum(data[, i], na.rm = TRUE)
      
    } else {
      if (is.null(targeted_norm_molecules)) {
        data$rank <- rank(-data[, i])
        norm_data <- data[data$Annotation %in% normalization_molecules, ]
      } else {
        if (start_row == 9) {
          targeted_norm_molecules$rank <- rank(-targeted_norm_molecules[, i-5])
        } else {
          targeted_norm_molecules$rank <- rank(-targeted_norm_molecules[, i+4])
        }
        norm_data <- targeted_norm_molecules[targeted_norm_molecules$Molecule.Name %in% normalization_molecules, ]
        norm_data$Annotation <- norm_data$Molecule.Name
      }
      
      # Get the area values of all norm. molecules from the sample
      # Use the max value if there are duplicates and and write to max_area column
      norm_data_max <- norm_data %>%
        group_by(Annotation) %>%
        slice_max(order_by = !!sym(sample_name)) %>%
        ungroup() %>%
        select(max_area = !!sym(sample_name))
      
      if (method == "non-weighted") {
        norm_factors[i - start_row] <- sum(norm_data_max$max_area, na.rm = TRUE)
        
      } else if (method == "Log2") {
        log_df <- log2(norm_data_max$max_area)
        mean_ <- mean(log_df, na.rm = TRUE)
        if (mean_ <= 0) {
          log_df <- log_df + (abs(mean_) + 1e2) # push all features up so the mean is positive 
        } 
        norm_factors[i - start_row] <- sum(log_df, na.rm = TRUE)
        
      } else if (method == "Log2_Sqrt") {
        norm_factors[i - start_row] <- sum(sqrt(log2(norm_data_max$max_area)), na.rm = TRUE)
        
      } else if (method == "Sqrt") {
        norm_factors[i - start_row] <- sum(sqrt(norm_data_max$max_area), na.rm = TRUE)
        
      } else if (method == "Rank") {
        norm_data_max <- norm_data_max %>%
          mutate(rank = rank(max_area))
        norm_factors[i - start_row] <- sum(norm_data_max$max_area * norm_data$rank, na.rm = TRUE)
        
      } else if (method == "Rank2") {
        norm_data_max <- norm_data_max %>%
          mutate(rank = rank(max_area))
        norm_factors[i - start_row] <- sum(norm_data_max$max_area * (norm_data$rank**2), na.rm = TRUE)
        
      } else if (method == "Weight_shifting") {
        median <- median(norm_data_max$max_area, na.rm = TRUE)
        subset_above_median <- norm_data_max[norm_data_max$max_area > median,] 
        subset_below_median <- norm_data_max[norm_data_max$max_area < median,] 
        sum <- sum(subset_above_median$max_area / ((subset_above_median$max_area - median) * 0.5), na.rm = TRUE)
        sum <- sum + sum(subset_below_median$max_area * ((median - subset_below_median$max_area) * 0.5), na.rm = TRUE)
        norm_factors[i - start_row] <- sum
        
      } else if (method == "Median") {
        median <- median(norm_data_max$max_area, na.rm = TRUE)
        norm_factors[i - start_row] <- sum(norm_data_max$max_area / median, na.rm = TRUE)
        
      } else if (method == "Med.Std.dev") { #basically a z-score
        median <- median(norm_data_max$max_area, na.rm = TRUE)
        stdDev <- sd(norm_data_max$max_area, na.rm = TRUE)
        norm_factors[i - start_row] <- sum((norm_data_max$max_area - median) / stdDev, na.rm = TRUE)
        
      } else if (method == "Mean") {
        mean <- mean(norm_data_max$max_area, na.rm = TRUE)
        norm_factors[i - start_row] <- sum(norm_data_max$max_area / mean, na.rm = TRUE)
        
      } 
    }
  }
  if ('rank' %in% names(data)) { # if there is still the column rank then remove this, would make problems in next line
    data <- subset(data, select = -c(rank))
  }
  names(norm_factors) <- colnames(data)[(start_row+1):ncol(data)]
  return(norm_factors)
}



# Function to normalize the data
normalize_data_with_normFactor <- function(data, norm_factors) {
  for (i in 5:(ncol(data)-1)) {
    sample_name <- colnames(data)[i]
    data[, i] <- data[, i] / norm_factors[sample_name]
    # print(sample_name)
    # print(data[, i])
    # print(norm_factors[sample_name])
    # print(data[, i] / norm_factors[sample_name])
  }
  return(data)
}



# Function to prepare data for plotting
prepare_data_for_plot <- function(data, add_info_file) {
  add_info <- read.csv(add_info_file)
  if (ncol(add_info) < 2) {
    add_info <- read.csv2(add_info_file)
  }
  
  # Reshape data to have samples as rows and features as columns
  if ("rtime_group" %in% colnames(data) ) {
    data <- pivot_longer(data, 
                                  cols = 10:ncol(data),
                                  names_to = "Sample", 
                                  values_to = "Area")
  } else {
    data <- pivot_longer(data,
                                  cols = 5:(ncol(data) - 1),
                                  names_to = "Sample",
                                  values_to = "Area")
  }
  
  data <- merge(data, add_info, by.x = "Sample", by.y = "Sample", all.x = TRUE)
  data$Area_log2 <- log2(data$Area)

  if("date" %in% colnames(data)) {
    data$date <- as.Date(data$date, format = "%d.%m.%Y")
  } 
  
  if("Date" %in% colnames(data)) {
    data$Date <- as.Date(data$Date, format = "%d.%m.%Y")
  } 
  
  if("time" %in% colnames(data)) {
    data$time <- sapply(data$time, convert_to_minutes)
    if ("Timepoint" %in% colnames(data)) {
      data$time[data$Timepoint == "24h"] <- 1440 # from csv get 00:00 for 24h so change "by hand" here
    }
  }
  
  if("Time" %in% colnames(data)) {
    data$Time <- sapply(data$Time, convert_to_minutes)
    if ("Timepoint" %in% colnames(data)) {
      data$Time[data$Timepoint == "24h"] <- 1440 # from csv get 00:00 for 24h so change "by hand" here
    }
    data$time <- data$Time
  }
  
  return(data)
}



# function to change the time to minutes between samples
convert_to_minutes <- function(time_str) {
  time_components <- strsplit(time_str, ":")[[1]]
  hours <- as.numeric(time_components[1])
  minutes <- as.numeric(time_components[2])
  total_minutes <- hours * 60 + minutes
  return(total_minutes)
}



# Function to replace NAs and Infs 
replace_na_and_Inf <- function(x) {
  if (is.numeric(x)) {
    x[is.na(x)] <- min(x, na.rm = TRUE)
    x[is.infinite(x)] <- max(x, na.rm = TRUE)
  }
  return(as.numeric(x))
}


  
###############################################
## Data Normalization and Plotting Functions ##
###############################################

# Function to Plot seceral molecules for several donors
  plot_molecule_curves <- function(resultsdir, full_prep_data, Group = 'Donor') {
    unique_molecules <- unique(full_prep_data$Molecule.Name)
    unique_donors <- unique(full_prep_data[[Group]])
    
    plot_list <- list()
    
    # Loop over each unique molecule
    for (molecule in unique_molecules) {
      # Filter data for the current molecule
      molecule_data <- subset(full_prep_data, Molecule.Name == molecule)
      
      # Create the ggplot
      p <- ggplot(molecule_data, aes(x = time, y = Area, color = !!sym(Group), group = !!sym(Group))) +
        geom_point() +
        geom_line(linewidth = 1) +  # Add lines for each donor
        labs(title = molecule, 
             x = "Time", y = "Area", color = Group) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5))  # Center the title
      
      # Add the plot to the list
      plot_list[[molecule]] <- p
    }
    
    # Save all plots to a single PDF with arranged grid
    pdf(file = paste0(resultsdir, "curve-molecule-check-for-normalization.pdf"), width = 25, height = 4*(length(unique_molecules) %/% 3) + 1) 
    do.call(grid.arrange, c(plot_list, ncol = 3))  # Arrange plots in a grid (3 per row)
    dev.off()
    
  }




# Function: Plot data over time for different normalization parameters
# Input: resultsdir (directory for result files), experiment_data, unique_molecules, normalization_molecules
# Output: PDF with plots
  plot_to_pdf_normalization_check <- function(resultsdir, experiment_data, unique_molecules, normalization_molecules) {
    # Set up PDF device
    pdf(file = paste0(resultsdir, "normalization-feature_check.pdf"), width = 15, height = 3*(length(unique_molecules) %/% 3) + 1) 
    
    # Loop over each unique molecule
    for (molecule in unique_molecules) {
      plot_list <- list()
      
      # subsets for the current molecule
      filtered_data <- subset(experiment_data, Molecule.Name == molecule)
      
      # raw data plot
      raw_plot <- ggplot(data = filtered_data, 
                         aes(x = time, y = Area, 
                             #color = interaction(paper, date, hand)
                             #color = interaction(Donor, Solvent, Hand)
                             color = interaction(Donor, Drink)
                             )) +
        geom_line() +
        geom_point(size = 3) +
        labs(x = "Time (minutes)", y = "Area", color = "Factors", 
             title = paste(molecule, "over time (Raw Data)")) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      # loop over each normalization variant
      for (norm_mol in normalization_molecules) {
        if (norm_mol %in% unique_molecules) {
          # subsets for normalization variant
          norm_data <- subset(experiment_data, Molecule.Name == norm_mol)
          
          # normalization
          filtered_data_norm <- merge(filtered_data, norm_data, by = "Sample", suffixes = c("", "_norm_mol"))
          filtered_data_norm$Normalized_Area <- filtered_data_norm$Area / filtered_data_norm$Area_norm_mol
          
          # normalized data plot
          normalized_plot <- ggplot(data = filtered_data_norm, 
                                    aes(x = time, y = Normalized_Area, 
                                        #color = interaction(paper, date, hand)
                                        #color = interaction(Donor, Solvent, Hand)
                                        color = interaction(Donor, Drink)
                                        )) +
            geom_line() +
            geom_point(size = 3) +
            labs(x = "Time (minutes)", y = "Normalized Area", color = "Factors", 
                 title = paste(norm_mol, "normalized", molecule, "over time")) +
            theme_minimal()
          
          plot_list[[length(plot_list) + 1]] <- normalized_plot
        }
      }
      
      # combine all plots for the current molecule into one plot
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    # Close PDF device
    dev.off()
  }
  
  
  
  # Function: Plot data over time for caffein for several donors 
  # Input: resultsdir (directory for result files), experiment_data, unique_molecules, normalization_molecules
  # Output: PDF with plots
  plot_to_pdf_normalization_check_per_person_caffein <- function(resultsdir, experiment_data, unique_molecules, normalization_molecules) {
    # Set up PDF device
    pdf(file = paste0(resultsdir, "normalization-feature_check-per-donor.pdf"), width = 15, height = 3*(length(unique_molecules) %/% 3) + 1) 
    
    unique_donors <- unique(experiment_data$Donor)
    molecule = "Caffein"
    
    # Loop over each unique molecule
    for (donor in unique_donors) {
      plot_list <- list()
      
      # subsets for the current molecule
      filtered_data <- subset(experiment_data, Donor == donor & Molecule.Name == molecule)
      
      # raw data plot
      raw_plot <- ggplot(data = filtered_data, 
                         aes(x = time, y = Area, 
                             #color = interaction(paper, date, hand)
                             color = interaction(Donor, Solvent, Hand)
                             #color = interaction(Donor, Drink)
                         )) +
        geom_line() +
        geom_point(size = 3) +
        labs(x = "Time (minutes)", y = "Area", color = "Factors", 
             title = paste(molecule, "over time (Raw Data) in Donor", donor)) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      # loop over each normalization variant
      for (norm_mol in normalization_molecules) {
        if (norm_mol %in% unique_molecules) {
          # subsets for normalization variant
          norm_data <- subset(experiment_data, Donor == donor & Molecule.Name == norm_mol)
          
          # normalization
          filtered_data_norm <- merge(filtered_data, norm_data, by = "Sample", suffixes = c("", "_norm_mol"))
          filtered_data_norm$Normalized_Area <- filtered_data_norm$Area / filtered_data_norm$Area_norm_mol
          
          # normalized data plot
          normalized_plot <- ggplot(data = filtered_data_norm, 
                                    aes(x = time, y = Normalized_Area, 
                                        #color = interaction(paper, date, hand)
                                        color = interaction(Donor, Solvent, Hand)
                                        #color = interaction(Donor, Drink)
                                    )) +
            geom_line() +
            geom_point(size = 3) +
            labs(x = "Time (minutes)", y = "Normalized Area", color = "Factors", 
                 title = paste(norm_mol, "normalized", molecule, "over time in Donor", donor)) +
            theme_minimal()
          
          plot_list[[length(plot_list) + 1]] <- normalized_plot
        }
      }
      
      # combine all plots for the current molecule into one plot
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    # Close PDF device
    dev.off()
  }
  
  
  
  
  # Function: Plot data for different normalization molecules to check the quality of normalization with different methods 
  # Input: resultsdir (directory for result files), norm_feature_list, all_needed_features, add_info, Group
  # Output: PDF with plots
  plot_to_pdf_normalization_check_objective_tests <- function(resultsdir, norm_feature_list, all_needed_features, add_info, Group, norm_methods = NULL) {
    
    # Reshape all_needed_features to have samples as rows and features as columns
    if ("rtime_group" %in% colnames(all_needed_features) ) {
      sample_names <- colnames(all_needed_features[10:length(all_needed_features)])
      reshaped_data <- pivot_longer(all_needed_features, 
                                    cols = -c(id, rt, mz, rtime_group, feature_group, Annotation, Formula, Confidence, SMILES),
                                    names_to = "variable", 
                                    values_to = "value")
    } else {
      sample_names <- colnames(all_needed_features[5:length(all_needed_features)])
      reshaped_data <- pivot_longer(all_needed_features,
                                    cols = -c(id, rt, mz, charge),
                                    names_to = "variable",
                                    values_to = "value")
    }
    
    reshaped_data <- dcast(reshaped_data, variable ~ id, value.var = "value")
    colnames(reshaped_data)[1] <- "Sample"
    rownames(reshaped_data) <- reshaped_data[,1]
    pure_data <- reshaped_data[,-1]
    
    nr_unique_molecules <- nrow(norm_feature_list)
    
    # Merge reshaped data with add_info
    merged_data <- merge(reshaped_data, add_info, by.x = "Sample", by.y = "Sample")
    
    merged_data[, 2:(length(merged_data) - 6)] <- as.data.frame(apply(merged_data[, 2:(length(merged_data) - 6)], 2, replace_na_and_Inf))
    pure_data <- as.data.frame(apply(pure_data, 2, replace_na_and_Inf))
    
    # Calculate sample counts
    sample_counts <- table(merged_data[[Group]])
    sample_labels <- paste0(names(sample_counts), "\n(n = ", sample_counts, ")")
    
    # Precompute the normalized dataframes and store in a list
    if (is.null(norm_methods)) {
      normalized_data_list <- normalized_data_with_single_molecules(norm_feature_list, merged_data, pure_data, sample_names, add_info)
    } else {
      normalized_data_list <- normalized_data_with_various_normalization_methods(norm_feature_list, merged_data, pure_data, norm_methods, sample_names, add_info)
    }
    
    # Set up PDF device
    if (is.null(norm_methods)) {
      pdf(file = paste0(resultsdir, "normalization-feature_check_objective-methods_", Group, ".pdf"), width = 25, height = 4*(nr_unique_molecules %/% 3) + 1) 
    } else {
      pdf(file = paste0(resultsdir, "normalization-method_check_objective-methods_", Group, ".pdf"), width = 25, height = 4*(nr_unique_molecules %/% 3) + 1) 
    }
    
    
    # List of comparison methods
    comparison_methods <- c("PCA", "PLSDA", "RLA", "RSD")
    
    # Loop over each comparison method
    for (method in comparison_methods) {
      plot_list <- list()
      
      # Plot raw data for the current comparison method
      if (method == "PCA") {
        pca_raw <- prcomp(pure_data, center = TRUE, scale. = TRUE)
        plot <- autoplot(pca_raw, data = merged_data, colour = Group, frame = TRUE, frame.type = 'norm') +
          labs(title = "PCA (Raw Data)") +
          theme_minimal()
        plot_list[[length(plot_list) + 1]] <- plot
        
      } else if (method == "PLSDA") {
        plsda_raw <- plsda(pure_data, as.factor(merged_data[[Group]]), ncomp = 2)
        plsda_scores_raw <- data.frame(Comp.1 = numeric(nrow(merged_data)),
                                       Comp.2 = numeric(nrow(merged_data)),
                                       Group = factor(nrow(merged_data)),
                                       stringsAsFactors = FALSE)
        plsda_scores_raw$Comp.1 <- plsda_raw$scores[,1]
        plsda_scores_raw$Comp.2 <- plsda_raw$scores[,2]
        plsda_scores_raw$Group <- merged_data[[Group]]
        
        plot <- ggplot(plsda_scores_raw, aes(x = Comp.1, y = Comp.2, color = !!sym("Group"))) +
          geom_point(size = 3) +
          labs(title = "PLSDA (Raw Data)", x = "Component 1", y = "Component 2") +
          theme_minimal()
        plot_list[[length(plot_list) + 1]] <- plot
        
      } else if (method == "RLA") {
        rla_data <- t(pure_data) # need to transform to get features to rows
        rla_raw <- rowRla(rla_data, f = merged_data[[Group]]) #makes rla row wise
        rla_long <- as.data.frame(rla_raw) %>%
          rownames_to_column(var = "Feature") %>% 
          pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")
        
        plot <- ggplot(rla_long, aes(x = Sample, y = Value)) +
          geom_boxplot(outlier.size = 0.1) +
          labs(title = "RLA (Raw Data)", x = "Samples", y = "Relative Log Abundance") +
          theme_minimal() +
          theme(axis.text.x = element_blank())
          #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        
        # Automatically adjust the size of axis text based on the number of samples
        # num_samples <- length(sample_names)
        # if (num_samples > 200) {
        #   plot <- plot + theme(axis.text.x = element_text(size = 4))
        # } else if (num_samples > 100) {
        #   plot <- plot + theme(axis.text.x = element_text(size = 6))
        # } else if (num_samples > 50) {
        #   plot <- plot + theme(axis.text.x = element_text(size = 8))
        # } else if (num_samples > 20) {
        #   plot <- plot + theme(axis.text.x = element_text(size = 10))
        # } else {
        #   plot <- plot + theme(axis.text.x = element_text(size = 12))
        # }
        
        plot_list[[length(plot_list) + 1]] <- plot
        
      } else if (method == "RSD") {
        rsd_raw <- calculate_rsd(merged_data, unique(merged_data[[Group]]), Group)
        plot <- ggplot(rsd_raw, aes(x = Group, y = RSD)) +
          geom_boxplot() +
          scale_x_discrete(labels = sample_labels) +
          labs(title = "RSD (Raw Data)", x = "Group", y = "Relative Standard Deviation (%)") +
          theme_minimal()
        #plot_list[[length(plot_list) + 1]] <- plot # dont plot separate but add to the normalized ones as direct comparison
      }
      
      # set what to use for normalization loop
      if (is.null(norm_methods)) {
        normalization_type <- norm_feature_list$id
      } else {
        normalization_type <- norm_methods
      }
      
      # Loop over each normalization variant
      for (norm_feature in normalization_type) {
        # Get the normalized data from precomputed list
        norm_data <- normalized_data_list[[norm_feature]]
        norm_merged_data <- norm_data$merged_data
        norm_pure_data <- norm_data$pure_data
        norm_molecule_name <- norm_data$molecule_name
        
        
        if (method == "PCA") {
          # PCA normalized data
          pca_norm <- prcomp(norm_pure_data, center = TRUE, scale. = TRUE)
          plot <- autoplot(pca_norm, data = norm_merged_data, colour = Group, frame = TRUE, frame.type = 'norm') +
            labs(title = paste0("PCA (Normalized with ", norm_molecule_name, ")")) +
            theme_minimal()
          plot_list[[length(plot_list) + 1]] <- plot
          
        } else if (method == "PLSDA") {
          # PLSDA normalized data
          plsda_norm <- plsda(norm_pure_data, as.factor(norm_merged_data[[Group]]), ncomp = 2)
          plsda_scores_norm <- data.frame(Comp.1 = numeric(nrow(norm_merged_data)),
                                         Comp.2 = numeric(nrow(norm_merged_data)),
                                         Group = factor(nrow(norm_merged_data)),
                                         stringsAsFactors = FALSE)
          plsda_scores_norm$Comp.1 <- plsda_norm$scores[,1]
          plsda_scores_norm$Comp.2 <- plsda_norm$scores[,2]
          plsda_scores_norm$Group <- norm_merged_data[[Group]]
          
          plot <- ggplot(plsda_scores_norm, aes(x = Comp.1, y = Comp.2, color = !!sym("Group"))) +
            geom_point(size = 3) +
            labs(title = paste0("PLSDA (Normalized  with ", norm_molecule_name, ")"), x = "Component 1", y = "Component 2") +
            theme_minimal()
          plot_list[[length(plot_list) + 1]] <- plot
          
        } else if (method == "RLA") {
          # RLA normalized data
          rla_norm_data <- t(norm_pure_data) # need to transform to get features to rows
          rla_norm <- rowRla(rla_norm_data, f = norm_merged_data[[Group]]) #makes rla row wise
          rla_norm_long <- as.data.frame(rla_norm) %>%
            rownames_to_column(var = "Feature") %>% 
            pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")
          
          plot <- ggplot(rla_norm_long, aes(x = Sample, y = Value)) +
            geom_boxplot(outlier.size = 0.1) +
            labs(title = paste0("RLA (Normalized with ", norm_molecule_name, ")"), x = "Samples", y = "Relative Log Abundance") +
            theme_minimal() +
            theme(axis.text.x = element_blank())
            #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
          
          # Automatically adjust the size of axis text based on the number of samples
          # num_samples <- length(sample_names)
          # if (num_samples > 200) {
          #   plot <- plot + theme(axis.text.x = element_text(size = 4))
          # } else if (num_samples > 100) {
          #   plot <- plot + theme(axis.text.x = element_text(size = 6))
          # } else if (num_samples > 50) {
          #   plot <- plot + theme(axis.text.x = element_text(size = 8))
          # } else if (num_samples > 20) {
          #   plot <- plot + theme(axis.text.x = element_text(size = 10))
          # } else {
          #   plot <- plot + theme(axis.text.x = element_text(size = 12))
          # }
          
          plot_list[[length(plot_list) + 1]] <- plot
          
        } else if (method == "RSD") {
          # RSD normalized data
          rsd_norm <- calculate_rsd(norm_merged_data, unique(norm_merged_data[[Group]]), Group)
          
          rsd_raw$Typ <- "raw"
          rsd_norm$Typ <- "norm"
          rsd_merged <- rbind(rsd_raw, rsd_norm)
          rsd_merged$Typ <- factor(rsd_merged$Typ, levels = c("raw", "norm"))
          
          # Perform one-sided t-test and collect p-values
          p_values <- sapply(unique(rsd_merged$Group), function(group) {
            raw_rsd <- subset(rsd_raw, Group == group)$RSD
            norm_rsd <- subset(rsd_norm, Group == group)$RSD
            t_test <- t.test(raw_rsd, norm_rsd, alternative = "greater", paired = TRUE)
            return(t_test$p.value)
          })
          
          # Prepare data frame for p-values
          p_values_df <- data.frame(Group = unique(rsd_merged$Group), p_value = p_values)
          
          plot <- ggplot(rsd_merged, aes(x = Group, y = RSD, fill = Typ)) +
            geom_boxplot(outlier.size = 0.5) +
            scale_x_discrete(labels = sample_labels) +
            labs(title = paste0("RSD (Normalized with ", norm_molecule_name, ")"), x = "Group", y = "Relative Standard Deviation (%)") +
            theme_minimal() +
            facet_wrap(~ Typ, scales = "free_x") +  # Separate panels for raw and normalized
            geom_text(data = p_values_df, aes(x = Group, y = max(rsd_merged$RSD) + 1, label = paste("p =", round(p_value, 4))), 
                      inherit.aes = FALSE, vjust = -0.5)
          
          plot_list[[length(plot_list) + 1]] <- plot
        }
      }
      
      # Combine all plots for the current comparison method into one plot
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    
    # Close PDF device
    dev.off()
  }
  
  
  
  
  # Function: Plot data over time for different normalization parameters
  # Input: resultsdir (directory for result files), experiment_data, unique_molecules, normalization_molecules
  # Output: PDF with plots
  plot_to_pdf_normalization_methods_check <- function(targeted_experiment_data, unique_molecules, normalization_molecules, 
                                                      all_needed_features, add_info_file, resultsdir, norm_feature_list = NULL) {
    
    # Normalization methods
    #norm_methods <- c("Tyrosine", "Total_Sum", "non-weighted", "Log2_Sqrt", "Sqrt", "Rank", "Rank2", "Weight_shifting")
    norm_methods <- c("Tyrosine", "Total_Sum", "non-weighted", "Sqrt", "Rank2", "Median", "Med.Std.dev", "Mean")
    
    norm_factors_list <- list()
    for (method in norm_methods) {
      norm_factors_list[[method]] <- calculate_norm_factors(all_needed_features, method, normalization_molecules, norm_feature_list)
    }
    
    # Set up PDF device
    pdf(file = paste0(resultsdir, "normalization-method-comparison.pdf"), width = 15, height = 8) 
    
    # Loop over each unique molecule
    for (molecule in unique_molecules) {
      plot_list <- list()
      
      # Subset for the current molecule
      filtered_data <- subset(targeted_experiment_data, Molecule.Name == molecule)
      prepared_data <- prepare_data_for_plot(filtered_data, add_info_file)
      
      # Raw data plot
      raw_plot <- ggplot(data = prepared_data, 
                         aes(x = time, y = Area_log2, 
                             color = interaction(paper, date, hand),
                             #color = interaction(Donor, Solvent, Hand)
                             )) +
                  geom_line() +
                  geom_point(size = 3) +
                  labs(x = "Time (minutes)", y = "log2(AUC)", color = "Factors", 
                       title = paste(molecule, "over time (Raw Data)")) +
                  theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      
      for (i in 1:length(norm_methods)) {
        filtered_data_norm <- normalize_data_with_normFactor(filtered_data, norm_factors_list[[i]])
        prepared_data <- prepare_data_for_plot(filtered_data_norm, add_info_file)
        
        normalized_plot <- ggplot(data = prepared_data, 
                                  aes(x = time, y = Area_log2, 
                                      color = interaction(paper, date, hand),
                                      #color = interaction(Donor, Solvent, Hand)
                                      )) +
          geom_line() +
          geom_point(size = 3) +
          labs(x = "Time (minutes)", y = "log2(nAUC)", color = "Factors", 
               title = paste(norm_methods[i], "normalized", molecule, "over time")) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- normalized_plot
      }
      
      
      # combine all plots for the current molecule into one plot and save to pdf
        do.call(grid.arrange, plot_list)
    }
    
    dev.off()
  }
  
    
    
  # Function: Plot data over time for different normalization parameters
  # Input: resultsdir (directory for result files), experiment_data, unique_molecules, normalization_molecules
  # Output: PDF with plots
  plot_to_pdf_normalization_methods_check_per_person <- function(targeted_experiment_data, unique_molecules, normalization_molecules, 
                                                      all_needed_features, add_info_file, resultsdir, norm_feature_list = NULL) {
    
    # Normalization methods
    norm_methods <- c("Tyrosine", "Total_Sum", "non-weighted", "Log2_Sqrt", "Sqrt", "Rank", "Rank2", "Weight_shifting")
    
    norm_factors_list <- list()
    for (method in norm_methods) {
      norm_factors_list[[method]] <- calculate_norm_factors(all_needed_features, method, normalization_molecules, norm_feature_list)
    }
    
    # Set up PDF device
    pdf(file = paste0(resultsdir, "normalization-method-comparison.pdf"), width = 15, height = 8) 
    
    # Loop over each unique molecule
    for (molecule in unique_molecules) {
      # Subset for the current molecule
      filtered_data <- subset(targeted_experiment_data, Molecule.Name == molecule)
      prep_data_all_persons <- prepare_data_for_plot(filtered_data, add_info_file)
      
      unique_person <- unique(prep_data_all_persons$Donor)
      
      for (person in unique_person) {
        plot_list <- list()
        
        filtered_data_per_person <- filtered_data %>% select(1:4, contains(person), length(filtered_data))
        prepared_data <- prepare_data_for_plot(filtered_data_per_person, add_info_file)
        
          # Raw data plot
          raw_plot <- ggplot(data = prepared_data, 
                             aes(x = time, y = Area_log2, 
                                 #color = interaction(paper, date, hand),
                                 #color = interaction(Solvent, Hand))) +
                                 color = interaction(Drink))) +
          geom_line() +
          geom_point(size = 3) +
          labs(x = "Time (minutes)", 
               y = "log2(AUC)", 
               #y = "AUC", 
               color = "Factors", 
               title = paste(molecule, "over time (Raw Data) from Donor", person)) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- raw_plot
        
        
        for (i in 1:length(norm_methods)) {
          filtered_data_norm <- normalize_data_with_normFactor(filtered_data_per_person, norm_factors_list[[i]])
          prepared_data <- prepare_data_for_plot(filtered_data_norm, add_info_file)
          
          normalized_plot <- ggplot(data = prepared_data, 
                                    aes(x = time, y = Area_log2, 
                                        #color = interaction(paper, date, hand),
                                        #color = interaction(Solvent, Hand))) +
                                        color = interaction(Drink))) +
            geom_line() +
            geom_point(size = 3) +
            labs(x = "Time (minutes)", 
                 y = "log2(AUC)", 
                 #y = "AUC", 
                 color = "Factors", 
                 title = paste(norm_methods[i], "normalized", molecule, "over time from Donor", person)) +
            theme_minimal()
          
          plot_list[[length(plot_list) + 1]] <- normalized_plot
        }
        
        
        # combine all plots for the current molecule into one plot and save to pdf
        do.call(grid.arrange, plot_list)
      }
    }
      
    dev.off()
  }
  
  
  # this part is to get the needed parts for the function below
    # targeted_experiment_data <- targeted_experiment_data_raw
    # all_needed_features <- raw_data
    # final_results <- final_raw_results
    # ranked_results <- ranked_raw_results %>% 
    #   mutate(Norm.ID.Method = paste(Norm.ID, "+", Norm.Method))  %>%
    #   filter(Norm.ID.Method %in% (filtered_counted_norm_combinations_raw %>%
    #                                 filter(count >= 2) %>%
    #                                 pull(Norm.ID.Method)))
  
  # Function to plot data over time for tested normalization combinations (check top 40 combinations)
  plot_to_pdf_normalization_combinations_check <- function(targeted_experiment_data, unique_molecules, add_info_file, all_needed_features,
                                                           norm_feature_list, resultsdir, final_results, ranked_results) {
    
    top5percent <- round((nrow(final_results) / length(unique(final_results$Donor))) * 0.1)
    ranked_results_topx <- rank_combinations(final_results, top_x =  top5percent)
    
    filtered_counted_norm_combinations <- ranked_results_topx %>%
      mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
      group_by(Norm.ID.Method) %>%
      summarise(count = n()) %>%
      ungroup() %>%
      arrange(desc(count)) %>%
      left_join(ranked_results_topx %>%
                  distinct(Norm.ID, Norm.Method, Norm.Molecule) %>%
                  mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                  select(-Norm.ID, -Norm.Method),
                by = "Norm.ID.Method")
    
    fwrite(filtered_counted_norm_combinations, paste0(resultsdir, "counted_norm_combinations_of_top", top5percent, ".csv"))
     
    # Rank normalization methods by residuals over all groups
    overall_ranked_results <- ranking_norm_methods_by_residuals_over_all_groups(final_results, ranked_results_topx, filtered_counted_norm_combinations)
    overall_ranked_results <- bind_rows(lapply(overall_ranked_results, as.data.frame))
    #overall_ranked_results <- overall_ranked_results[, 1:4]
    overall_ranked_results <- overall_ranked_results %>%
      distinct(norm_id, norm_method, norm_molecules)
    
    # Calculate normalization factors for top 40 normalization methods
    top_40_norm_methods <- overall_ranked_results[1:40, ]
    
    norm_factors_list <- lapply(1:nrow(top_40_norm_methods), function(i) {
      method <- top_40_norm_methods$norm_method[i]
      molecules <- unlist(strsplit(top_40_norm_methods$norm_molecules[i], "\\+"))
      norm_factors <- calculate_norm_factors(all_needed_features, method, molecules, norm_feature_list)
      list(norm_id = top_40_norm_methods$norm_id[i], norm_method = method, norm_factors = norm_factors)
    })
    
    pdf(file = paste0(resultsdir, "normalization-combination_check.pdf"), width = 15, height = 70)
    
    for (molecule in unique_molecules) {
      plot_list <- list()
      
      # Subset for the current molecule
      filtered_data <- subset(targeted_experiment_data_vsn, Molecule.Name == molecule)
      prepared_data <- prepare_data_for_plot(filtered_data, add_info_file)
      
      # Raw data plot
      raw_plot <- ggplot(data = prepared_data, 
                         aes(x = time, y = Area_log2, 
                             color = interaction(paper, date, hand))) +
        geom_line() +
        geom_point(size = 3) +
        labs(x = "Time (minutes)", y = "log2(AUC)", color = "Factors", 
             title = paste(molecule, "over time (Raw Data)")) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      # Plot each normalization method
      for (i in 1:length(norm_factors_list)) {
        filtered_data_norm <- normalize_data_with_normFactor(filtered_data, norm_factors_list[[i]]$norm_factors)
        prepared_data_norm <- prepare_data_for_plot(filtered_data_norm, add_info_file)
        
        normalized_plot <- ggplot(data = prepared_data_norm, 
                                  aes(x = time, y = Area_log2, 
                                      color = interaction(paper, date, hand))) +
          geom_line() +
          geom_point(size = 3) +
          labs(x = "Time (minutes)", y = "log2(nAUC)", color = "Factors", 
               title = paste("Normalized", molecule, "over time\n(Norm ID:", norm_factors_list[[i]]$norm_id, ", Norm Method:", norm_factors_list[[i]]$norm_method)) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- normalized_plot
      }
      
      # Combine all plots for the current molecule into one plot and save to pdf
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    
    dev.off()
  }
  
  
  
  # Bateman Function with 1 compartment
  bateman <- function(t, ka, kel, Vd=60) {
    BA <- 1
    D <- 80 
    (D * BA)/Vd * (ka / (ka - kel))  * (exp(-kel * t) - exp(-ka * t))
  }
  
  # Bateman Function, as a delayed model with 1 compartments
  bateman_delayed_1_compartment <- function(t, ka, kel, tlag) {
    # Concentration after absorption delay
    C <- ifelse(
      t < tlag,
      0, # Before the lag time, concentration is zero
      (ka / (ka - kel)) * (exp(-kel * (t - tlag)) - exp(-ka * (t - tlag)))
    )
    return(C)
  }
  
  
  # Bateman Function, as a transit model with 2 compartments
  bateman_transit_2_compartment <- function(t, ka, k12, kel) {
    # Blood + interstitial (Compartment 1)
    C1 <- ka / (ka - k12) * (exp(-k12 * t) - exp(-ka * t))
    
    # Sweat (Compartment 2)
    C2 <- C1 * (k12 / (k12 - kel)) * (exp(-(kel) * t) - exp(-(k12) * t))
    
    return(C2)
  }
  
  
  # Bateman Function, as a transit model with 2 compartments
  bateman_2_compartment_peripheral_no_reflux <- function(t, ka, kel, k12, kel2, Vt = 15, D = 80, BA = 1) {
    # Ensure numerical stability by avoiding division by zero or near-zero differences
    if (abs(ka - (kel + k12)) < 1e-6 || abs((k12 - kel2)) < 1e-6) {
      stop("Numerical instability detected. Adjust parameters.")
    }
    
    # Calculate drug concentration in peripheral compartment (compartment 2)
    # Equation for C1(t) (using exponential solution for first-order decay)
    C1 <- (D * BA * ka / (ka - (kel + k12))) * (exp(-(kel + k12) * t) - exp(-(ka) * t))  # Assuming initial concentration C0 = 0 at t=0
    
    # C2(t) is driven by the transfer from C1 and elimination from C2
    C2 <- (C1 / Vt) * (k12 / (k12 - kel2)) * (exp(-(kel2) * t) - exp(-(k12) * t)) 
    
    # Set concentration to zero if it drops below a small threshold (to avoid numerical artifacts)
    C2[C2 < 0] <- 0
    
    return(C2)
  }
  
  
  
  # Function to Normalize and Fit Caffeine Data
  plot_to_pdf_various_normalization_fits_for_biological_kinetics <- function(targeted_experiment_data, norm_feature_list, resultsdir, add_info_file, molecule) {
    
    full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
    full_prep_norm_list <- prepare_data_for_plot(norm_feature_list, add_info_file)
    unique_donors <- unique(full_prep_data$Donor)
    norm_molecules <- norm_feature_list$Molecule.Name
    
    #unique_donors <- unique_donors[-1] # remove first one as for this experiment the first one has outliers
    
    # Set up PDF device
    pdf(file = paste0(resultsdir, "normalization-quality-check_caffeine-curve-fitting.pdf"), width = 15, height = 3*(length(norm_molecules) %/% 3) + 1)
    
    # Data frame for results
    results <- data.frame(matrix(ncol = 10, nrow = 0))
    colnames(results) <- c('Norm.Molecule', 'Donor', 'Res.stand.err', 'Res.stand.err.norm',
                           'ka', 'kel', 'Dose', 'ka.norm', 'kel.norm', 'Dose.norm')
    i <- 1
    
    for (donor in unique_donors) {
      plot_list <- list()
      
      # Subset for current donor and for the wanted molecule
      donor_data <- subset(full_prep_data, Donor == donor & Molecule.Name == molecule)
      donor_data$scaled_conc <- donor_data$Area - min(donor_data$Area) # put lowest point to 0
      donor_data$scaled_conc <- donor_data$scaled_conc / max(donor_data$scaled_conc[donor_data$time > 0 & donor_data$time < 100]) # second point should be the highest, normalize this to 1
      
      # add a fixed point for caffeine that after 10 hours the caffeine is fully metabolized
      # add_line_for_fixed_point <- as.data.frame(donor_data[nrow(donor_data),])
      # add_line_for_fixed_point$Sample <- "added_fixed_point"
      # add_line_for_fixed_point$time <- 600
      # add_line_for_fixed_point$scaled_conc <- 0
      # donor_data <- rbind(donor_data, add_line_for_fixed_point)
      # 
      # # set weights so the last point is weighted strongest
      we <- rep(1, nrow(donor_data))
      # we[nrow(donor_data)] <- 10
      

      # Fit raw data
      fit <- nlsLM(scaled_conc ~ bateman(time, ka, kel, D), 
                   data = donor_data,
                   start = list(ka = 0.03, kel = 0.02, D = 80),
                   lower = c(ka = 0.01, kel = 0.01, D = 60),
                   upper = c(ka = 0.05, kel = 0.03, D = 100),
                   weights = we,
                   control = list(maxiter = 100))
      
      fitted_data <- data.frame(time = seq(0, max(donor_data$time)))
      fitted_data$scaled_conc <- predict(fit, list(time = seq(0, max(donor_data$time))))
      summary_model <- summary(fit)
      
      # Raw data plot
      raw_plot <- ggplot(data = donor_data, aes(x = time)) +
        #geom_point(aes(y = scaled_conc, color = interaction(Donor, paper, hand)), size = 3) +
        geom_point(aes(y = scaled_conc, color = interaction(Donor, Drink)), size = 3) +
        #geom_line(aes(y = scaled_conc, color = interaction(Donor, paper, hand))) +
        geom_line(data = fitted_data, aes(x = time, y = scaled_conc), color = "red") +
        annotate("text", x = Inf, y = Inf, label = paste("RSE:", round(summary_model$sigma, 3)), 
                 hjust = 1.1, vjust = 1.1, size = 5) +
        labs(x = "Time (minutes)", y = "scaled Concentration", color = "Factors", 
             title = paste(molecule, "over time (Raw Data) in Donor", donor)) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      # Loop over each normalization variant
      for (norm_mol in normalization_molecules) {
        norm_data <- subset(full_prep_norm_list, Donor == donor & Molecule.Name == norm_mol)
        if (nrow(norm_data) > 0) {
          donor_data_norm <- merge(donor_data, norm_data, by = "Sample", suffixes = c("", "_norm"))
          donor_data_norm$Normalized_Area <- donor_data_norm$Area / donor_data_norm$Area_norm
          donor_data_norm$scaled_conc_norm <- donor_data_norm$Normalized_Area - min(donor_data_norm$Normalized_Area) # put lowest point to 0
          donor_data_norm$scaled_conc_norm <- donor_data_norm$scaled_conc_norm / max(donor_data_norm$scaled_conc_norm[donor_data_norm$time > 0 & donor_data_norm$time < 100]) # second point should be the highest, normalize this to 1
          
          # # add a fixed point for caffeine that after 10 hours the caffeine is fully metabolized
          # add_line_for_fixed_point <- as.data.frame(donor_data_norm[nrow(donor_data_norm),])
          # add_line_for_fixed_point$Sample <- "added_fixed_point"
          # add_line_for_fixed_point$time <- 600
          # add_line_for_fixed_point$scaled_conc_norm <- 0
          # donor_data_norm <- rbind(donor_data_norm, add_line_for_fixed_point)
          # 
          # # set weights so the last point is weighted strongest
          # we <- rep(1, nrow(donor_data_norm))
          # we[nrow(donor_data_norm)] <- 10
          
          # Fit normalized data
          fit_norm <- try(nlsLM(scaled_conc_norm ~ bateman(time, ka, kel, D), data = donor_data_norm,
                                start = list(ka = 0.03, kel = 0.02, D = 80),
                                lower = c(ka = 0.01, kel = 0.01, D = 60),
                                upper = c(ka = 0.05, kel = 0.03, D = 100),
                                weights = we,
                                control = list(maxiter = 100)),
                          silent = TRUE)
          
          if (class(fit_norm) != "try-error") {
            fitted_norm_data <- data.frame(time = seq(0, max(donor_data_norm$time)))
            fitted_norm_data$scaled_conc_norm <- predict(fit_norm, list(time = seq(0, max(donor_data_norm$time))))
            summary_model_norm <- summary(fit_norm)
            
            # Plot normalized data with fitted curve
            normalized_plot <- ggplot(data = donor_data_norm, aes(x = time)) +
              #geom_point(aes(y = scaled_conc_norm, color = interaction(Donor, paper, hand)), size = 3) +
              geom_point(aes(y = scaled_conc_norm, color = interaction(Donor, Drink)), size = 3) +
              #geom_line(aes(y = Normalized_Area, color = interaction(Donor, paper, hand))) +
              geom_line(data = fitted_norm_data, aes(x = time, y = scaled_conc_norm), color = "red") +
              annotate("text", x = Inf, y = Inf, label = paste("RSE:", round(summary_model_norm$sigma, 3)), 
                       hjust = 1.1, vjust = 1.1, size = 5) +
              labs(x = "Time (minutes)", y = "Normalized scaled Conc", color = "Factors", 
                   title = paste(norm_mol, "normalized", molecule, "in Donor", donor)) +
              theme_minimal()
            
            plot_list[[length(plot_list) + 1]] <- normalized_plot
            
            
            results[i,] <- c(norm_mol,
                             donor, 
                             as.numeric(summary_model$sigma), 
                             as.numeric(summary_model_norm$sigma), 
                             as.numeric(coef(fit)['ka']),
                             as.numeric(coef(fit)['kel']),
                             as.numeric(coef(fit)['D']),
                             as.numeric(coef(fit_norm)['ka']),
                             as.numeric(coef(fit_norm)['kel']),
                             as.numeric(coef(fit_norm)['D']))
            i <- i + 1
          }
        }
      }
      
      # Combine all plots for the current donor into one plot
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    
    # Close PDF device
    dev.off()
    
    # Save results to CSV
    write.csv2(results, paste0(resultsdir, "results_normalization-quality-check.csv"))
  } 
  
  
  
  # Function to rank combinations based on fit quality
  rank_combinations <- function(results, top_x = 20) {
    # Extract non-list columns and filter out NA and Inf values in Res.stand.err.norm
    filtered_non_list_cols <- results %>%
      select(-Fit.Data, -Fit.Curve) %>%
      filter(!is.na(Res.stand.err.norm) & is.finite(Res.stand.err.norm))
    
    # Rank the combinations and select top 20 for each donor
    top_hits <- filtered_non_list_cols %>%
      group_by(Donor) %>%
      arrange(Res.stand.err.norm) %>%
      slice_head(n = top_x) %>%
      ungroup()
    
    # Reattach the list columns
    top_hits <- top_hits %>%
      left_join(results %>% select(Donor, Norm.Molecule, Norm.Method, Fit.Data, Fit.Curve),
                by = c("Donor", "Norm.Molecule", "Norm.Method"))
    
    return(top_hits)
  }
  
  
  # Function to rank combinations based on fit quality and get top 1000 for each donor
  rank_top1000_combinations <- function(results) {
    # Extract non-list columns and filter out NA and Inf values in Res.stand.err.norm
    filtered_non_list_cols <- results %>%
      select(-Fit.Data, -Fit.Curve) %>%
      filter(!is.na(Res.stand.err.norm) & is.finite(Res.stand.err.norm))
    
    # Rank the combinations and select top 20 for each donor
    top_1000 <- filtered_non_list_cols %>%
      group_by(Donor) %>%
      arrange(Res.stand.err.norm) %>%
      slice_head(n = 1000) %>%
      ungroup()
    
    # Reattach the list columns
    top_1000 <- top_1000 %>%
      left_join(results %>% select(Donor, Norm.Molecule, Norm.Method, Fit.Data, Fit.Curve),
                by = c("Donor", "Norm.Molecule", "Norm.Method"))
    
    return(top_1000)
  }
  
  
  
  # Function to rank all methods over all groups by the residuals over all groups
  ranking_norm_methods_by_residuals_over_all_groups <- function(final_results, ranked_results, filtered_counted_norm_combinations, cl_given = NULL) {
    # Identify unique normalization methods that appear in the top x for any donor
    unique_normalizations <- ranked_results %>%
      distinct(Norm.ID, Norm.Method, Norm.Molecule) %>%
      mutate(Norm.ID.Method = paste(Norm.ID, "+", Norm.Method)) %>% 
      filter(Norm.ID.Method %in% (filtered_counted_norm_combinations %>%
                                    filter(count >= 4) %>%
                                    pull(Norm.ID.Method)))
    
    # filter out only the necessary Norm.IDs
    final_results_filtered <- final_results %>% 
      mutate(Norm.ID.Method = paste(Norm.ID, "+", Norm.Method)) %>%
      filter(Norm.ID.Method %in% unique_normalizations$Norm.ID.Method) 
    
    print("Filtering done...")
    
    # Preparation for parallelization
    cl <- cl_given
    if(is.null(cl_given)) {
      cl <- makeCluster(10)
      on.exit(stopCluster(cl), add = TRUE)  # Ensure the cluster is stopped if created here
    }
    
    clusterEvalQ(cl, {
      library(dplyr)
      library(tidyr)
    })
    clusterExport(cl, varlist = c('unique_normalizations', 'final_results_filtered'), envir = environment())
    
    print("Starting calculation...")
    
    # Initialize list to store sum of residuals for each normalization method
    norm_residuals_summary <- parLapply(cl, 1:nrow(unique_normalizations), function(i) {
      norm_id <- unique_normalizations$Norm.ID[i]
      norm_method <- unique_normalizations$Norm.Method[i]
      norm_molecules <- unique_normalizations$Norm.Molecule[i]
      
      # Gather data for this normalization method
      norm_data <- final_results_filtered %>%
        filter(Norm.ID == norm_id & Norm.Method == norm_method)
      
      residuals_list <- list()
      
      # Calculate residuals for each donor
      for (j in 1:nrow(norm_data)) {
        donor <- norm_data$Donor[j]
        fit_data_norm <- norm_data$Fit.Data[[j]]
        fit_curve_norm <- norm_data$Fit.Curve[[j]]
        
        residuals <- fit_data_norm$scaled_conc_norm - fit_curve_norm$scaled_conc_norm[fit_curve_norm$time %in% fit_data_norm$time]
        residuals_list[[j]] <- data.frame(time = fit_data_norm$time, residuals = residuals, donor = donor)
      }
      
      residuals_data <- bind_rows(residuals_list)
      
      # Calculate the mean of residuals for this normalization method
      total_residuals_mean <- mean(abs(residuals_data$residuals))
      total_residuals_sd <- sd(residuals_data$residuals)
      list(norm_id = norm_id, norm_method = norm_method, norm_molecules = norm_molecules, 
           total_residuals_mean = total_residuals_mean, total_residuals_sd = total_residuals_sd, residuals_data = residuals_data)
    })
    
    print("Finished calculation...")
    
    # Sort normalization methods by the mean of residuals
    norm_residuals_summary <- norm_residuals_summary[order(sapply(norm_residuals_summary, function(x) x$total_residuals_sd))]

    # Convert the summary to a data frame for easier handling
    #norm_residuals_summary_df <- bind_rows(lapply(norm_residuals_summary, as.data.frame))
    
    return(norm_residuals_summary)
  }
  
  
  
  # Helper function to apply normalization factors
  apply_normalization <- function(data, norm_factors) {
    normalized_data <- data
    normalized_data$Area_norm <- sweep(as.matrix(data$Area), 1, as.matrix(norm_factors), "/")
    return(normalized_data)
  }
  
  
  
  # Function to plot the top 20 combinations including the raw data
  plot_top_20_curves <- function(ranked_results, raw_results, resultsdir, molecule, pqn_results = NULL, vsn_results = NULL) {
    pdf(file = paste0(resultsdir, "top_20_normalization_curves_for_", molecule, "_fit.pdf"), 
        width = 15, height = 3*((nrow(ranked_results) / length(unique(ranked_results$Donor))) %/% 3 + 1))
    
    for (donor in unique(ranked_results$Donor)) {
      donor_results <- ranked_results %>% filter(Donor == donor)
      raw_donor_results <- raw_results %>% filter(Donor == donor)
      pqn_donor_results <- pqn_results %>% filter(Donor == donor)
      vsn_donor_results <- vsn_results %>% filter(Donor == donor)
      plot_list <- list()
      
      # Plot the raw data
      fit_data <- raw_donor_results$Fit.Data[[1]]
      fit_curve <- raw_donor_results$Fit.Curve[[1]]
      
      raw_plot <- ggplot(data = fit_data, aes(x = time)) +
        geom_point(aes(y = scaled_conc), size = 2) +
        geom_line(data = fit_curve, aes(x = time, y = scaled_conc), color = "red") +
        labs(x = "Time (minutes)", y = "Scaled Concentration", 
             title = paste0("Donor: ", donor, " - Raw Data Fit\n(RSE: ", round(raw_donor_results$Res.stand.err[1], 3), ")")) +
        theme_minimal()
      
      plot_list[[length(plot_list) + 1]] <- raw_plot
      
      # Plot the pqn data
      if (!is.null(pqn_results)) {
        # Plot the raw data first
        fit_data_pqn <- pqn_donor_results$Fit.Data[[1]]
        fit_curve_pqn <- pqn_donor_results$Fit.Curve[[1]]
        
        pqn_plot <- ggplot(data = fit_data_pqn, aes(x = time)) +
          geom_point(aes(y = scaled_conc), size = 2) +
          geom_line(data = fit_curve_pqn, aes(x = time, y = scaled_conc), color = "red") +
          labs(x = "Time (minutes)", y = "Scaled Concentration", 
               title = paste0("Donor: ", donor, " - PQN Data Fit\n(RSE: ", round(pqn_donor_results$Res.stand.err[1], 3), ")")) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- pqn_plot
      }
      
      # Plot the vsn data
      if (!is.null(vsn_results)) {
        # Plot the raw data first
        fit_data_vsn <- vsn_donor_results$Fit.Data[[1]]
        fit_curve_vsn <- vsn_donor_results$Fit.Curve[[1]]
        
        vsn_plot <- ggplot(data = fit_data_vsn, aes(x = time)) +
          geom_point(aes(y = scaled_conc), size = 2) +
          geom_line(data = fit_curve_vsn, aes(x = time, y = scaled_conc), color = "red") +
          labs(x = "Time (minutes)", y = "Scaled Concentration", 
               title = paste0("Donor: ", donor, " - VSN Data Fit\n(RSE: ", round(vsn_donor_results$Res.stand.err[1], 3), ")")) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- vsn_plot
      }
      
      
      
      # Plot the top 20 normalized data fits
      for (i in 1:nrow(donor_results)) {
        fit_data_norm <- donor_results$Fit.Data[[i]]
        fit_curve_norm <- donor_results$Fit.Curve[[i]]
        
        plot <- ggplot(data = fit_data_norm, aes(x = time)) +
          geom_point(aes(y = scaled_conc_norm), size = 2) +
          geom_line(data = fit_curve_norm, aes(x = time, y = scaled_conc_norm), color = "red") +
          labs(x = "Time (minutes)", y = "Scaled Concentration", 
               title = paste0("Donor: ", donor, " - Combination: ", donor_results$Norm.ID[i], "\n", donor_results$Norm.Method[i], " (RSE: ", round(donor_results$Res.stand.err.norm[i], 3), ")")) +
          theme_minimal()
        
        plot_list[[length(plot_list) + 1]] <- plot
      }
      
      do.call(grid.arrange, c(plot_list, ncol = 3))
    }
    dev.off()
  }
  
  
  
  # Function to plot the top 20 combinations including the raw data
  plot_top_20_residuals <- function(ranked_results, final_results, final_raw_results, resultsdir, molecule, final_pqn_results = NULL, final_vsn_results = NULL) {
    # Identify unique normalization methods that appear in the top 20 for any donor
    unique_normalizations <- ranked_results %>%
      distinct(Norm.ID, Norm.Method, Norm.Molecule)
    
    #Define how many plots per page
    plots_per_page <- 60
    
    pdf(file = paste0(resultsdir, "top_20_normalization_residuals_for_", molecule, "_fit.pdf"), 
        width = 15, height = 4 * (plots_per_page %/% 3) + 1)
    
    plot_list <- list()
    
    
    # Plot residuals for raw data
    residuals_list <- list()
    raw_data_donors <- unique(final_raw_results$Donor)
    
    for (donor in raw_data_donors) {
      raw_data <- final_raw_results %>% filter(Donor == donor)
      fit_data <- raw_data$Fit.Data[[1]]
      fit_curve <- raw_data$Fit.Curve[[1]]
      
      residuals <- fit_data$scaled_conc - fit_curve$scaled_conc[fit_curve$time %in% fit_data$time]
      residuals_list[[donor]] <- data.frame(time = fit_data$time, residuals = residuals, donor = donor)
    }
    
    residuals_data <- do.call(rbind, residuals_list)
    
    # Plot residuals for raw data
    raw_residual_plot <- ggplot(data = residuals_data, aes(x = time, y = residuals, color = donor, group = donor)) +
      geom_point() +
      geom_line(linewidth = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      labs(x = "Time (minutes)", y = "Residuals", 
           title = "Raw Data Residuals") +
      theme_minimal() +
      annotate("text", x = Inf, y = Inf, label = paste("Sum of Residuals:", round(sum(residuals_data$residuals), 3), 
                                                       "\nSD of Residuals:", round(sd(residuals_data$residuals), 3)), 
               hjust = 1.1, vjust = 1.1, size = 5)
    
    plot_list[[length(plot_list) + 1]] <- raw_residual_plot
    
    
    
    # Plot residuals for pqn data
    if (!is.null(final_pqn_results)) {
      residuals_list <- list()
      pqn_data_donors <- unique(final_pqn_results$Donor)
      
      for (donor in pqn_data_donors) {
        pqn_data <- final_pqn_results %>% filter(Donor == donor)
        fit_data <- pqn_data$Fit.Data[[1]]
        fit_curve <- pqn_data$Fit.Curve[[1]]
        
        residuals <- fit_data$scaled_conc - fit_curve$scaled_conc[fit_curve$time %in% fit_data$time]
        residuals_list[[donor]] <- data.frame(time = fit_data$time, residuals = residuals, donor = donor)
      }
      
      residuals_data <- do.call(rbind, residuals_list)
      
      # Plot residuals for pqn data
      pqn_residual_plot <- ggplot(data = residuals_data, aes(x = time, y = residuals, color = donor, group = donor)) +
        geom_point() +
        geom_line(linewidth = 0.1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(x = "Time (minutes)", y = "Residuals", 
             title = "PQN Data Residuals") +
        theme_minimal() +
        annotate("text", x = Inf, y = Inf, label = paste("Sum of Residuals:", round(sum(residuals_data$residuals), 3), 
                                                         "\nSD of Residuals:", round(sd(residuals_data$residuals), 3)), 
                 hjust = 1.1, vjust = 1.1, size = 5)
      
      plot_list[[length(plot_list) + 1]] <- pqn_residual_plot
    }
    
    
    
    # Plot residuals for vsn data
    if (!is.null(final_vsn_results)) {
      residuals_list <- list()
      vsn_data_donors <- unique(final_vsn_results$Donor)
      
      for (donor in vsn_data_donors) {
        vsn_data <- final_vsn_results %>% filter(Donor == donor)
        fit_data <- vsn_data$Fit.Data[[1]]
        fit_curve <- vsn_data$Fit.Curve[[1]]
        
        residuals <- fit_data$scaled_conc - fit_curve$scaled_conc[fit_curve$time %in% fit_data$time]
        residuals_list[[donor]] <- data.frame(time = fit_data$time, residuals = residuals, donor = donor)
      }
      
      residuals_data <- do.call(rbind, residuals_list)
      
      # Plot residuals for vsn data
      vsn_residual_plot <- ggplot(data = residuals_data, aes(x = time, y = residuals, color = donor, group = donor)) +
        geom_point() +
        geom_line(linewidth = 0.1) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        labs(x = "Time (minutes)", y = "Residuals", 
             title = "VSN Data Residuals") +
        theme_minimal() +
        annotate("text", x = Inf, y = Inf, label = paste("Sum of Residuals:", round(sum(residuals_data$residuals), 3), 
                                                         "\nSD of Residuals:", round(sd(residuals_data$residuals), 3)), 
                 hjust = 1.1, vjust = 1.1, size = 5)
      
      plot_list[[length(plot_list) + 1]] <- vsn_residual_plot
    }
    
    
    
    # Initialize list to store sum of residuals for each normalization method
    norm_residuals_summary <- list()
    
    # Loop through each unique normalization method
    for (i in 1:nrow(unique_normalizations)) {
      norm_id <- unique_normalizations$Norm.ID[i]
      norm_method <- unique_normalizations$Norm.Method[i]
      norm_molecules <- unique_normalizations$Norm.Molecule[i]
        
      # Gather data for this normalization method
      norm_data <- final_results %>%
        filter(Norm.ID == norm_id & Norm.Method == norm_method)
      
      residuals_list <- list()

      # Calculate residuals for each donor
      for (j in 1:nrow(norm_data)) {
        donor <- norm_data$Donor[j]
        fit_data_norm <- norm_data$Fit.Data[[j]]
        fit_curve_norm <- norm_data$Fit.Curve[[j]]
        
        residuals <- fit_data_norm$scaled_conc_norm - fit_curve_norm$scaled_conc_norm[fit_curve_norm$time %in% fit_data_norm$time]
        residuals_list[[j]] <- data.frame(time = fit_data_norm$time, residuals = residuals, donor = donor)
      }
      
      residuals_data <- do.call(rbind, residuals_list)
      
      # Calculate the mean of residuals for this normalization method
      total_residuals <- mean(abs(residuals_data$residuals))
      norm_residuals_summary[[i]] <- list(norm_id = norm_id, norm_method = norm_method, norm_molecules = norm_molecules, total_residuals = total_residuals, residuals_data = residuals_data)
    }
    
    # Sort normalization methods by the mean of residuals
    #norm_residuals_summary <- norm_residuals_summary[order(sapply(norm_residuals_summary, function(x) x$total_residuals))]
    norm_residuals_summary <- norm_residuals_summary[order(sapply(norm_residuals_summary, function(x) sd(x$residuals_data$residuals)))]
    
    # Plot residuals for each sorted normalization method
    for (norm_residual in norm_residuals_summary) {
      residual_plot <- ggplot(data = norm_residual$residuals_data, aes(x = time, y = residuals, color = donor, group = donor)) +
        geom_point() +
        geom_line(linewidth = 0.5) +
        labs(x = "Time (minutes)", y = "Residuals", 
             title = paste("Normalization:", norm_residual$norm_id, "\n", norm_residual$norm_method)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        theme_minimal() +
        if (length(unique(norm_residual$residuals_data$donor)) < length(raw_data_donors)) {
        annotate("text", x = Inf, y = Inf, label = paste("Mean of Residuals:", round(norm_residual$total_residuals, 3), 
                                                         "\nSD of Residuals:", round(sd(norm_residual$residuals_data$residuals), 3)), 
                 hjust = 1.1, vjust = 1.1, size = 5, color = "red")
        } else {
          annotate("text", x = Inf, y = Inf, label = paste("Mean of Residuals:", round(norm_residual$total_residuals, 3), 
                                                           "\nSD of Residuals:", round(sd(norm_residual$residuals_data$residuals), 3)), 
                   hjust = 1.1, vjust = 1.1, size = 5)
        }
      
      plot_list[[length(plot_list) + 1]] <- residual_plot
    }
    
    # Calculate number of pages needed
    num_pages <- ceiling(length(plot_list) / plots_per_page)
    
    # Create a for loop to arrange and save each page
    for (page in 1:num_pages) {
      # Determine which plots to include in this page
      start_index <- (page - 1) * plots_per_page + 1
      end_index <- min(page * plots_per_page, length(plot_list))
      plots_to_arrange <- plot_list[start_index:end_index]
      
      # Arrange the plots using grid.arrange
      grid.arrange(grobs = plots_to_arrange, ncol = 3)
    }
    
    dev.off()
  }
  
  
  
  # Function to plot histograms and scatter/density plots for given data
  plot_histograms_and_scatter_density <- function(final_results, final_raw_results, resultsdir, molecule) {
    
    donors <- unique(final_results$Donor)
    
    # Create output PDF
    pdf(file = paste0(resultsdir, "overview_of_curve_parameters_after_normalization_combinations.pdf"), height = 15, width = 15)
    
    # Define boundaries 
    lower <- c(ka = 0.009, kel = 0.009, Vd = 24)
    upper <- c(ka = 0.051, kel = 0.031, Vd = 46)
    
    # Plot raw data for all donors
    # Histograms
    p1 <- ggplot(final_raw_results, aes(x = ka, fill = factor(Donor))) + 
      geom_histogram(binwidth = 0.001, alpha = 0.7) + 
      labs(title = paste("Histogram of ka.raw for", molecule, "curves\nAll Donors"), x = "ka", y = "Frequency") + 
      xlim(lower["ka"], upper["ka"]) +
      theme_minimal()
    
    p2 <- ggplot(final_raw_results, aes(x = kel, fill = factor(Donor))) + 
      geom_histogram(binwidth = 0.001, alpha = 0.7) + 
      labs(title = paste("Histogram of kel.raw for", molecule, "curves\nAll Donors"), x = "kel", y = "Frequency") + 
      xlim(lower["kel"], upper["kel"]) +
      theme_minimal()
    
    p3 <- ggplot(final_raw_results, aes(x = Vd, fill = factor(Donor))) + 
      geom_histogram(binwidth = 1, alpha = 0.7) + 
      labs(title = paste("Histogram of Vd.raw for", molecule, "curves\nAll Donors"), x = "Vd", y = "Frequency") + 
      xlim(lower["Vd"], upper["Vd"]) +
      theme_minimal()
    
    p4 <- p1 + scale_y_continuous(trans='log10')
    
    p5 <- p2 + scale_y_continuous(trans='log10')
    
    p6 <- p3 + scale_y_continuous(trans='log10')
    
    # Scatter plots with density
    p7 <- ggplot(final_raw_results, aes(x = ka, y = kel)) + 
      geom_bin2d(bins = 100) + 
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of ka.raw vs kel.raw for", molecule, "curves\nAll Donors"), x = "ka", y = "kel") + 
      xlim(lower["ka"], upper["ka"]) +
      ylim(lower["kel"], upper["kel"]) +
      theme_minimal()
    
    p8 <- ggplot(final_raw_results, aes(x = ka, y = Vd)) + 
      geom_bin2d(bins = 100) +        
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of ka.raw vs Vd.raw for", molecule, "curves\nAll Donors"), x = "ka", y = "Vd") + 
      xlim(lower["ka"], upper["ka"]) +
      ylim(lower["Vd"], upper["Vd"]) +
      theme_minimal()
    
    p9 <- ggplot(final_raw_results, aes(x = kel, y = Vd)) + 
      geom_bin2d(bins = 100) +        
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of kel.norm vs Vd.norm for", molecule, "curves\nAll Donors"), x = "kel", y = "Vd") + 
      xlim(lower["kel"], upper["kel"]) +
      ylim(lower["Vd"], upper["Vd"]) +
      theme_minimal()
    
    # Arrange plots in a grid
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
    
    
    
    # Plot normalization results for all donors
    # Histograms
    p1 <- ggplot(final_results, aes(x = ka.norm, fill = factor(Donor))) +
      geom_histogram(binwidth = 0.001, alpha = 0.7) +
      labs(title = paste("Histogram of ka.norm for", molecule, "curves\nAll Donors"), x = "ka.norm", y = "Frequency") +
      xlim(lower["ka"], upper["ka"]) +
      theme_minimal()

    p2 <- ggplot(final_results, aes(x = kel.norm, fill = factor(Donor))) +
      geom_histogram(binwidth = 0.001, alpha = 0.7) +
      labs(title = paste("Histogram of kel.norm for", molecule, "curves\nAll Donors"), x = "kel.norm", y = "Frequency") +
      xlim(lower["kel"], upper["kel"]) +
      theme_minimal()

    p3 <- ggplot(final_results, aes(x = Vd.norm, fill = factor(Donor))) +
      geom_histogram(binwidth = 1, alpha = 0.7) +
      labs(title = paste("Histogram of Vd.norm for", molecule, "curves\nAll Donors"), x = "Vd.norm", y = "Frequency") +
      xlim(lower["Vd"], upper["Vd"]) +
      theme_minimal()
    
    p4 <- p1 + scale_y_continuous(trans='log10')
    
    p5 <- p2 + scale_y_continuous(trans='log10')
    
    p6 <- p3 + scale_y_continuous(trans='log10')

    # Scatter plots with density
    p7 <- ggplot(final_results, aes(x = ka.norm, y = kel.norm)) +
      geom_bin2d(bins = 100) +
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of ka.norm vs kel.norm for", molecule, "curves\nAll Donors"), x = "ka.norm", y = "kel.norm") +
      xlim(lower["ka"], upper["ka"]) +
      ylim(lower["kel"], upper["kel"]) +
      theme_minimal()

    p8 <- ggplot(final_results, aes(x = ka.norm, y = Vd.norm)) +
      geom_bin2d(bins = 100) +
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of ka.norm vs Vd.norm for", molecule, "curves\nAll Donors"), x = "ka.norm", y = "Vd.norm") +
      xlim(lower["ka"], upper["ka"]) +
      ylim(lower["Vd"], upper["Vd"]) +
      theme_minimal()

    p9 <- ggplot(final_results, aes(x = kel.norm, y = Vd.norm)) +
      geom_bin2d(bins = 100) +
      scale_fill_gradientn(colours = rainbow(4)) +
      labs(title = paste("Scatter plot of kel.norm vs Vd.norm for", molecule, "curves\nAll Donors"), x = "kel.norm", y = "Vd.norm") +
      xlim(lower["kel"], upper["kel"]) +
      ylim(lower["Vd"], upper["Vd"]) +
      theme_minimal()

    # Arrange plots in a grid
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
    
    
    
    # Plot normalization results for each donor
    for (donor in donors) {
      donor_data <- final_results %>% filter(Donor == donor)

      # Histograms
      p1 <- ggplot(donor_data, aes(x = ka.norm)) +
        geom_histogram(binwidth = 0.001, fill = 'blue', alpha = 0.7) +
        labs(title = paste("Histogram of ka.norm for", molecule, "curves\nDonor:", donor), x = "ka.norm", y = "Frequency") +
        xlim(lower["ka"], upper["ka"]) +
        theme_minimal()

      p2 <- ggplot(donor_data, aes(x = kel.norm)) +
        geom_histogram(binwidth = 0.001, fill = 'green', alpha = 0.7) +
        labs(title = paste("Histogram of kel.norm for", molecule, "curves\nDonor:", donor), x = "kel.norm", y = "Frequency") +
        xlim(lower["kel"], upper["kel"]) +
        theme_minimal()

      p3 <- ggplot(donor_data, aes(x = Vd.norm)) +
        geom_histogram(binwidth = 1, fill = 'red', alpha = 0.7) +
        labs(title = paste("Histogram of Vd.norm for", molecule, "curves\nDonor:", donor), x = "Vd.norm", y = "Frequency") +
        xlim(lower["Vd"], upper["Vd"]) +
        theme_minimal()
      
      p4 <- p1 + scale_y_continuous(trans='log10')
      
      p5 <- p2 + scale_y_continuous(trans='log10')
      
      p6 <- p3 + scale_y_continuous(trans='log10')

      # Scatter plots with density
      p7 <- ggplot(donor_data, aes(x = ka.norm, y = kel.norm)) +
        geom_bin2d(bins = 100) +
        scale_fill_gradientn(colours = rainbow(4)) +
        labs(title = paste("Scatter plot of ka.norm vs kel.norm for", molecule, "curves\nDonor:", donor), x = "ka.norm", y = "kel.norm") +
        xlim(lower["ka"], upper["ka"]) +
        ylim(lower["kel"], upper["kel"]) +
        theme_minimal()

      p8 <- ggplot(donor_data, aes(x = ka.norm, y = Vd.norm)) +
        geom_bin2d(bins = 100) +
        scale_fill_gradientn(colours = rainbow(4)) +
        labs(title = paste("Scatter plot of ka.norm vs Vd.norm for", molecule, "curves\nDonor:", donor), x = "ka.norm", y = "Vd.norm") +
        xlim(lower["ka"], upper["ka"]) +
        ylim(lower["Vd"], upper["Vd"]) +
        theme_minimal()

      p9 <- ggplot(donor_data, aes(x = kel.norm, y = Vd.norm)) +
        geom_bin2d(bins = 100) +
        scale_fill_gradientn(colours = rainbow(4)) +
        labs(title = paste("Scatter plot of kel.norm vs Vd.norm for", molecule, "curves\nDonor:", donor), x = "kel.norm", y = "Vd.norm") +
        xlim(lower["kel"], upper["kel"]) +
        ylim(lower["Vd"], upper["Vd"]) +
        theme_minimal()

      # Arrange plots in a grid
      grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol = 3)
    }
    
    # Close PDF
    dev.off()
  }
  
  
  
  # Function to plot a histogram of different pre-normalization strategies
  plot_histogram_of_residuals_of_different_normalization_strategies <- 
    function(top1000_raw, top1000_pqn, top1000_vsn, resultsdir, molecule) {
    
      # Add a new column to each dataframe indicating the group
      top1000_raw <- top1000_raw %>% mutate(Group = "raw")
      top1000_pqn <- top1000_pqn %>% mutate(Group = "pqn")
      top1000_vsn <- top1000_vsn %>% mutate(Group = "vsn")
      
      # Row bind the dataframes
      all_normalizations <- bind_rows(top1000_raw, top1000_pqn, top1000_vsn)
     
      # Get unique donors
      donors <- unique(all_normalizations$Donor)
      
      
      pdf(file = paste0(resultsdir, "histogram_of_residuals_of_different_normalization_strategies_for_", molecule, "_fit.pdf"), 
          width = 15, height = 4 * (length(donors) %/% 3) + 1)
      
      plot_list <- list()
      
      combined_plot <- ggplot(all_normalizations, aes(x = Res.stand.err.norm, fill = factor(Group))) +
        geom_histogram(alpha = 0.5, position = 'identity') + 
        labs(title = paste("Histogram of Res.stand.err.norm for", molecule, "curves for different normalization strategies\nAll Donors"), 
             x = "Res.stand.err.norm", y = "Frequency") +
        theme_minimal()
      
      plot_list[[1]] <- combined_plot
      
      
      # Plot histogram for each donor separately
      for (i in seq_along(donors)) {
        donor <- donors[i]
        donor_data <- all_normalizations %>% filter(Donor == donor)
        
        donor_plot <- ggplot(donor_data, aes(x = Res.stand.err.norm, fill = factor(Group))) +
          geom_histogram(alpha = 0.5, position = 'identity') + 
          labs(title = paste("Histogram of Res.stand.err.norm for", molecule, "curves for different normalization strategies\nDonor:", donor), 
               x = "Res.stand.err.norm", y = "Frequency") +
          theme_minimal()
        
        plot_list[[i + 1]] <- donor_plot
      }
      
      do.call(grid.arrange, c(plot_list, ncol = 3))
      
      
      plot_list <- list()
      
      combined_plot <- ggplot(all_normalizations, aes(x = Res.stand.err.norm, fill = factor(Group))) +
        geom_histogram(alpha = 0.7) + 
        labs(title = paste("Histogram of Res.stand.err.norm for", molecule, "curves for different normalization strategies\nAll Donors"), 
             x = "Res.stand.err.norm", y = "Frequency") +
        theme_minimal()
      
      plot_list[[1]] <- combined_plot
      
      
      # Plot histogram for each donor separately
      for (i in seq_along(donors)) {
        donor <- donors[i]
        donor_data <- all_normalizations %>% filter(Donor == donor)
        
        donor_plot <- ggplot(donor_data, aes(x = Res.stand.err.norm, fill = factor(Group))) +
          geom_histogram(alpha = 0.7) + 
          labs(title = paste("Histogram of Res.stand.err.norm for", molecule, "curves for different normalization strategies\nDonor:", donor), 
               x = "Res.stand.err.norm", y = "Frequency") +
          theme_minimal()
        
        plot_list[[i + 1]] <- donor_plot
      }
      
      do.call(grid.arrange, c(plot_list, ncol = 3))
      
      dev.off()
  }
  
  
  
  # Plot boxplot of different normalization strategies 
  # take only the top hits that worked for all donors for each strategie and compare them
  plot_boxplot_of_methods_with_best_residuals <- function(final_raw_results, final_pqn_results, final_vsn_results, resultsdir, molecule) {
    
    # Rank the normalization methods and get top20 of each donor back
    ranked_raw_results <- rank_combinations(final_raw_results, top_x = 10000)
    ranked_pqn_results <- rank_combinations(final_pqn_results, top_x = 10000)
    ranked_vsn_results <- rank_combinations(final_vsn_results, top_x = 10000)
    
    filtered_counted_norm_combinations_raw <- ranked_raw_results %>%
                                              mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                              group_by(Norm.ID.Method) %>%
                                              summarise(count = n()) %>%
                                              ungroup() %>%
                                              #filter(count >= max(count) - 2) %>%
                                              arrange(desc(count)) %>%
                                              left_join(ranked_raw_results %>%
                                                          distinct(Norm.ID, Norm.Method, Norm.Molecule) %>%
                                                          mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                                          select(-Norm.ID, -Norm.Method),
                                                        by = "Norm.ID.Method")
    
    filtered_counted_norm_molecules_raw <- ranked_raw_results %>%
                                          group_by(Norm.ID) %>%
                                          summarise(count = n()) %>%
                                          ungroup() %>%
                                          #filter(count >= max(count) - 2) %>%
                                          arrange(desc(count)) %>%
                                          left_join(ranked_raw_results %>%
                                                      distinct(Norm.ID, Norm.Molecule),
                                                    by = "Norm.ID")
    
    filtered_counted_norm_combinations_pqn <- ranked_pqn_results %>%
                                              mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                              group_by(Norm.ID.Method) %>%
                                              summarise(count = n()) %>%
                                              ungroup() %>%
                                              #filter(count >= max(count) - 2) %>%
                                              arrange(desc(count)) %>%
                                              left_join(ranked_pqn_results %>%
                                                          distinct(Norm.ID, Norm.Method, Norm.Molecule) %>%
                                                          mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                                          select(-Norm.ID, -Norm.Method),
                                                        by = "Norm.ID.Method")
    
    filtered_counted_norm_molecules_pqn <- ranked_pqn_results %>%
                                          group_by(Norm.ID) %>%
                                          summarise(count = n()) %>%
                                          ungroup() %>%
                                          #filter(count >= max(count) - 2) %>%
                                          arrange(desc(count)) %>%
                                          left_join(ranked_pqn_results %>%
                                                      distinct(Norm.ID, Norm.Molecule),
                                                    by = "Norm.ID")
    
    filtered_counted_norm_combinations_vsn <- ranked_vsn_results %>%
                                              mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                              group_by(Norm.ID.Method) %>%
                                              summarise(count = n()) %>%
                                              ungroup() %>%
                                              #filter(count >= max(count) - 2) %>%
                                              arrange(desc(count)) %>%
                                              left_join(ranked_vsn_results %>%
                                                          distinct(Norm.ID, Norm.Method, Norm.Molecule) %>%
                                                          mutate(Norm.ID.Method = paste(Norm.ID, Norm.Method, sep = " + ")) %>%
                                                          select(-Norm.ID, -Norm.Method),
                                                        by = "Norm.ID.Method")
    
    filtered_counted_norm_molecules_vsn <- ranked_vsn_results %>%
                                          group_by(Norm.ID) %>%
                                          summarise(count = n()) %>%
                                          ungroup() %>%
                                          #filter(count >= max(count) - 2) %>%
                                          arrange(desc(count)) %>%
                                          left_join(ranked_vsn_results %>%
                                                      distinct(Norm.ID, Norm.Molecule),
                                                    by = "Norm.ID")
    
    fwrite(filtered_counted_norm_combinations_raw, paste0(resultsdir, "counted_norm_combinations_of_top10000_raw.csv"))
    fwrite(filtered_counted_norm_combinations_pqn, paste0(resultsdir, "counted_norm_combinations_of_top10000_pqn.csv"))
    fwrite(filtered_counted_norm_combinations_vsn, paste0(resultsdir, "counted_norm_combinations_of_top10000_vsn.csv"))
    
    count_filtered_counted_norm_combinations_raw <- filtered_counted_norm_combinations_raw %>%
      group_by(count) %>%
      summarise(count_count = n()) %>%
      ungroup() %>%
      arrange(desc(count)) %>%
      mutate(count = as.character(count)) %>%
      bind_rows(data.frame(count = "sum", count_count = sum(.$count_count)))
    
    count_filtered_counted_norm_combinations_pqn <- filtered_counted_norm_combinations_pqn %>%
      group_by(count) %>%
      summarise(count_count = n()) %>%
      ungroup() %>%
      arrange(desc(count))%>%
      mutate(count = as.character(count)) %>%
      bind_rows(data.frame(count = "sum", count_count = sum(.$count_count)))
    
    count_filtered_counted_norm_combinations_vsn <- filtered_counted_norm_combinations_vsn %>%
      group_by(count) %>%
      summarise(count_count = n()) %>%
      ungroup() %>%
      arrange(desc(count)) %>%
      mutate(count = as.character(count)) %>%
      bind_rows(data.frame(count = "sum", count_count = sum(.$count_count)))
    
    fwrite(count_filtered_counted_norm_combinations_raw, paste0(resultsdir, "count_of_counted_norm_combinations_of_top10000_raw.csv"))
    fwrite(count_filtered_counted_norm_combinations_pqn, paste0(resultsdir, "count_of_counted_norm_combinations_of_top10000_pqn.csv"))
    fwrite(count_filtered_counted_norm_combinations_vsn, paste0(resultsdir, "count_of_counted_norm_combinations_of_top10000_vsn.csv"))
    
    molecule_counts_count_bigger6_norm_comb_raw  <- filtered_counted_norm_combinations_raw %>%
      filter(count > 6) %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_bigger6_norm_comb_pqn  <- filtered_counted_norm_combinations_pqn %>%
      filter(count > 6) %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_bigger6_norm_comb_vsn  <- filtered_counted_norm_combinations_vsn %>%
      filter(count > 6) %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_bigger9_norm_comb_pqn  <- filtered_counted_norm_combinations_pqn %>%
      filter(count > 9) %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_bigger9_norm_comb_vsn  <- filtered_counted_norm_combinations_vsn %>%
      filter(count > 9) %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    fwrite(molecule_counts_count_bigger6_norm_comb_raw, paste0(resultsdir, "count_of_molecules_of_top10000_bigger6donors_raw.csv"))
    fwrite(molecule_counts_count_bigger6_norm_comb_pqn, paste0(resultsdir, "count_of_molecules_of_top10000_bigger6donors_pqn.csv"))
    fwrite(molecule_counts_count_bigger6_norm_comb_vsn, paste0(resultsdir, "count_of_molecules_of_top10000_bigger6donors_vsn.csv"))
    fwrite(molecule_counts_count_bigger9_norm_comb_pqn, paste0(resultsdir, "count_of_molecules_of_top10000_bigger9donors_pqn.csv"))
    fwrite(molecule_counts_count_bigger9_norm_comb_vsn, paste0(resultsdir, "count_of_molecules_of_top10000_bigger9donors_vsn.csv"))
    
    
    weighting_counts_count_bigger6_norm_comb_raw  <- filtered_counted_norm_combinations_raw %>%
      filter(count > 6) %>%
      mutate(Norm.Method = str_extract(Norm.ID.Method, "(?<=\\+ ).*")) %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_bigger6_norm_comb_pqn  <- filtered_counted_norm_combinations_pqn %>%
      filter(count > 6) %>%
      mutate(Norm.Method = str_extract(Norm.ID.Method, "(?<=\\+ ).*")) %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_bigger6_norm_comb_vsn  <- filtered_counted_norm_combinations_vsn %>%
      filter(count > 6) %>%
      mutate(Norm.Method = str_extract(Norm.ID.Method, "(?<=\\+ ).*")) %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_bigger9_norm_comb_pqn  <- filtered_counted_norm_combinations_pqn %>%
      filter(count > 9) %>%
      mutate(Norm.Method = str_extract(Norm.ID.Method, "(?<=\\+ ).*")) %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_bigger9_norm_comb_vsn  <- filtered_counted_norm_combinations_vsn %>%
      filter(count > 9) %>%
      mutate(Norm.Method = str_extract(Norm.ID.Method, "(?<=\\+ ).*")) %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    fwrite(weighting_counts_count_bigger6_norm_comb_raw, paste0(resultsdir, "count_of_weighting_of_top10000_bigger6donors_raw.csv"))
    fwrite(weighting_counts_count_bigger6_norm_comb_pqn, paste0(resultsdir, "count_of_weighting_of_top10000_bigger6donors_pqn.csv"))
    fwrite(weighting_counts_count_bigger6_norm_comb_vsn, paste0(resultsdir, "count_of_weighting_of_top10000_bigger6donors_vsn.csv"))
    fwrite(weighting_counts_count_bigger9_norm_comb_pqn, paste0(resultsdir, "count_of_weighting_of_top10000_bigger9donors_pqn.csv"))
    fwrite(weighting_counts_count_bigger9_norm_comb_vsn, paste0(resultsdir, "count_of_weighting_of_top10000_bigger9donors_vsn.csv"))
    
    
    
    
    # Get a ranking of the best methods by checking the residues from all Donors with this method
    cl <- makeCluster(60)
    norm_residuals_summary_raw <- ranking_norm_methods_by_residuals_over_all_groups(final_raw_results, ranked_raw_results, filtered_counted_norm_combinations_raw, cl)
    norm_residuals_summary_pqn <- ranking_norm_methods_by_residuals_over_all_groups(final_pqn_results, ranked_pqn_results, filtered_counted_norm_combinations_pqn, cl)
    norm_residuals_summary_vsn <- ranking_norm_methods_by_residuals_over_all_groups(final_vsn_results, ranked_vsn_results, filtered_counted_norm_combinations_vsn, cl)
    stopCluster(cl)
    
    # Function to filter the lists based on the number of unique donors
    filter_lists <- function(data_list) {
      filtered_list <- lapply(data_list, function(sub_list) {
        unique_donors <- unique(sub_list$residuals_data$donor)
        if (length(unique_donors) >= 12) {
          return(sub_list)
        } else {
          return(NULL)
        }
      })
      filtered_list <- Filter(Negate(is.null), filtered_list)
      return(filtered_list)
    }
    
    # Combine the data from the three lists into a single data frame
    combine_data <- function(data_list, group_name) {
      data <- lapply(data_list, function(sub_list) {
        data.frame(total_residuals_sd = sub_list$total_residuals_sd, total_residuals_mean = sub_list$total_residuals_mean, Group = group_name,
                   norm_id = sub_list$norm_id, norm_method = sub_list$norm_method, norm_molecule = sub_list$norm_molecule)
      })
      do.call(bind_rows, data)
    }
    
    # Filter out the ones that dont have data from all Donors, as they dont work for all Donors
    norm_residuals_summary_raw_filtered <- filter_lists(norm_residuals_summary_raw)
    norm_residuals_summary_pqn_filtered <- filter_lists(norm_residuals_summary_pqn)
    norm_residuals_summary_vsn_filtered <- filter_lists(norm_residuals_summary_vsn)
    
    # Combine the data and give them a Group name
    raw_data <- combine_data(norm_residuals_summary_raw_filtered, "Raw") %>% 
      mutate(Norm.ID.Method = paste(norm_id, norm_method, sep = " + "))  %>% 
      left_join(filtered_counted_norm_combinations_raw, by = "Norm.ID.Method")
    pqn_data <- combine_data(norm_residuals_summary_pqn_filtered, "PQN") %>% 
      mutate(Norm.ID.Method = paste(norm_id, norm_method, sep = " + "))  %>% 
      left_join(filtered_counted_norm_combinations_pqn, by = "Norm.ID.Method")
    vsn_data <- combine_data(norm_residuals_summary_vsn_filtered, "VSN") %>% 
      mutate(Norm.ID.Method = paste(norm_id, norm_method, sep = " + "))  %>% 
      left_join(filtered_counted_norm_combinations_vsn, by = "Norm.ID.Method")
    
    # Combine all into one data frame
    all_data <- bind_rows(raw_data, pqn_data, vsn_data)
    
    # Filter Top100
    raw_data_top100 <- raw_data[1:100,]
    pqn_data_top100 <- pqn_data[1:100,]
    vsn_data_top100 <- vsn_data[1:100,]
    
    # Combine all into one data frame
    all_data_top100 <- bind_rows(raw_data_top100, pqn_data_top100, vsn_data_top100)
    
    # Filter Top10
    raw_data_top10 <- raw_data[1:10,]
    pqn_data_top10 <- pqn_data[1:10,]
    vsn_data_top10 <- vsn_data[1:10,]
    
    # Combine all into one data frame
    all_data_top10 <- bind_rows(raw_data_top10, pqn_data_top10, vsn_data_top10)
    
    
    # Filter Top20
    raw_data_top20 <- raw_data[1:20,]
    pqn_data_top20 <- pqn_data[1:20,]
    vsn_data_top20 <- vsn_data[1:20,]
    
    # Combine all into one data frame
    all_data_top20 <- bind_rows(raw_data_top20, pqn_data_top20, vsn_data_top20)
    
    
    molecule_counts_count_top10_of_top10000_norm_comb_raw  <- raw_data_top10 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_top20_of_top10000_norm_comb_raw  <- raw_data_top20 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_top10_of_top10000_norm_comb_pqn  <- pqn_data_top10 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_top20_of_top10000_norm_comb_pqn  <- pqn_data_top20 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_top10_of_top10000_norm_comb_vsn  <- vsn_data_top10 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    molecule_counts_count_top20_of_top10000_norm_comb_vsn  <- vsn_data_top20 %>%
      mutate(norm_molecule = strsplit(norm_molecule, "\\+")) %>%
      unnest(norm_molecule) %>%
      group_by(norm_molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    fwrite(molecule_counts_count_top10_of_top10000_norm_comb_raw, paste0(resultsdir, "count_of_molecules_of_top10-of-top10000_raw.csv"))
    fwrite(molecule_counts_count_top20_of_top10000_norm_comb_raw, paste0(resultsdir, "count_of_molecules_of_top20-of-top10000_raw.csv"))
    fwrite(molecule_counts_count_top10_of_top10000_norm_comb_pqn, paste0(resultsdir, "count_of_molecules_of_top10-of-top10000_pqn.csv"))
    fwrite(molecule_counts_count_top20_of_top10000_norm_comb_pqn, paste0(resultsdir, "count_of_molecules_of_top20-of-top10000_pqn.csv"))
    fwrite(molecule_counts_count_top10_of_top10000_norm_comb_vsn, paste0(resultsdir, "count_of_molecules_of_top10-of-top10000_vsn.csv"))
    fwrite(molecule_counts_count_top20_of_top10000_norm_comb_vsn, paste0(resultsdir, "count_of_molecules_of_top20-of-top10000_vsn.csv"))
    
    
    weighting_counts_count_top10_of_top10000_norm_comb_raw  <- raw_data_top10 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_top20_of_top10000_norm_comb_raw  <- raw_data_top20 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_top10_of_top10000_norm_comb_pqn  <- pqn_data_top10 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_top20_of_top10000_norm_comb_pqn  <- pqn_data_top20 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_top10_of_top10000_norm_comb_vsn  <- vsn_data_top10 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts_count_top20_of_top10000_norm_comb_vsn  <- vsn_data_top20 %>%
      group_by(norm_method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    fwrite(weighting_counts_count_top10_of_top10000_norm_comb_raw, paste0(resultsdir, "count_of_weighting_of_top10-of-top10000_raw.csv"))
    fwrite(weighting_counts_count_top20_of_top10000_norm_comb_raw, paste0(resultsdir, "count_of_weighting_of_top20-of-top10000_raw.csv"))
    fwrite(weighting_counts_count_top10_of_top10000_norm_comb_pqn, paste0(resultsdir, "count_of_weighting_of_top10-of-top10000_pqn.csv"))
    fwrite(weighting_counts_count_top20_of_top10000_norm_comb_pqn, paste0(resultsdir, "count_of_weighting_of_top20-of-top10000_pqn.csv"))
    fwrite(weighting_counts_count_top10_of_top10000_norm_comb_vsn, paste0(resultsdir, "count_of_weighting_of_top10-of-top10000_vsn.csv"))
    fwrite(weighting_counts_count_top20_of_top10000_norm_comb_vsn, paste0(resultsdir, "count_of_weighting_of_top20-of-top10000_vsn.csv"))
    
    
    
    # Plot boxplot to pdf
    pdf(file = paste0(resultsdir, "boxplot_of_top10000_biggerequal-count6_normalization_residuals_for_different_normalization_strategies_for_", molecule, "_fit_colored.pdf"), 
        width = 12, height = 32)
    
    plot_list <- list()
    
    
    # Calculate the number of points for each group
    counts <- as.data.frame(table(all_data$Group))
    colnames(counts) <- c("Group", "Count")
    
    # Plot SD of Residuals
    plot <- ggplot(all_data, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group, log10 scaled",
           x = "Group",
           y = "SD of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group, log10 scaled",
           x = "Group",
           y = "Mean of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot SD of Residuals
    plot <- ggplot(all_data, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group",
           x = "Group",
           y = "SD of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group",
           x = "Group",
           y = "Mean of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    
    
    
    # Calculate the number of points for each group
    counts <- as.data.frame(table(all_data_top100$Group))
    colnames(counts) <- c("Group", "Count")
    
    # Plot SD of Residuals
    plot <- ggplot(all_data_top100, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group of Top100, log10 scaled",
           x = "Group",
           y = "SD of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data_top100, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group of Top100, log10 scaled",
           x = "Group",
           y = "Mean of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot SD of Residuals
    plot <- ggplot(all_data_top100, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group of Top100",
           x = "Group",
           y = "SD of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data_top100, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group of Top100",
           x = "Group",
           y = "Mean of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    
    
    # Calculate the number of points for each group
    counts <- as.data.frame(table(all_data_top10$Group))
    colnames(counts) <- c("Group", "Count")
    
    # Plot SD of Residuals
    plot <- ggplot(all_data_top10, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group of Top10, log10 scaled",
           x = "Group",
           y = "SD of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data_top10, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group of Top10, log10 scaled",
           x = "Group",
           y = "Mean of Total Residuals") +
      scale_y_log10() +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot SD of Residuals
    plot <- ggplot(all_data_top10, aes(x = Group, y = total_residuals_sd, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of SD of Total Residuals by Group of Top10",
           x = "Group",
           y = "SD of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Plot mean of Residuals
    plot <- ggplot(all_data_top10, aes(x = Group, y = total_residuals_mean, fill = Group)) +
      geom_boxplot(alpha = 0.5) +
      geom_jitter(aes(color = count), position = position_jitter(width = 0.2), size = 1.5, alpha = 0.7) +
      labs(title = "Boxplot of Mean of Total Residuals by Group of Top10",
           x = "Group",
           y = "Mean of Total Residuals") +
      theme_minimal() +
      scale_color_gradient(low = "blue", high = "red")
    
    plot <- plot + geom_text(data = counts, aes(x = Group, y = 0.05, label = paste("n =", Count)), vjust = -0.5)
    show(plot)
    plot_list[[length(plot_list)+1]] <- plot
    
    # Arrange the plots
    do.call(grid.arrange, c(plot_list, ncol = 2))
    
    dev.off()
  }
  
  
  # Function to split a list into chunks
  split_chunks_ <- function(lst, n) {
    split_size <- ceiling(length(lst) / n)
    split(lst, rep(1:n, each = split_size, length.out = length(lst)))
  }
  
  # Function to split a list into chunks
  split_chunks <- function(lst, n) {
    split(lst, cut(seq_along(lst), n, labels = FALSE))
  }
  
  
  
  # function to scale the concentration
  scale_concentration <- function(donor_data, norm = FALSE) {
    if (norm == FALSE) {
      donor_data$scaled_conc <- donor_data$Area - min(donor_data$Area)
      donor_data$scaled_conc <- donor_data$scaled_conc / max(donor_data$scaled_conc[donor_data$time > 0 & donor_data$time < 100])
      return(donor_data$scaled_conc)
    } else {
      donor_data$scaled_conc_norm <- donor_data$Area_norm - min(donor_data$Area_norm)
      donor_data$scaled_conc_norm <- donor_data$scaled_conc_norm / max(donor_data$scaled_conc_norm[donor_data$time > 0 & donor_data$time < 100])
      return(donor_data$scaled_conc)
    }
  }
  
  
  
  
  # Start and initialize a cluster
  initialize_cluster <- function(nCore = detectCores(), varlist) {
    # Check to dont start cluster with more cores then are available
    nCore <- min(nCore - 1, detectCores() - 1)
    
    cl <- makeCluster(nCore - 1)   # Create a cluster with the desired number of cores
    
    clusterEvalQ(cl, {
      library(dplyr)
      library(minpack.lm)
      library(nlstools)
      library(tidyr)
      library(ggplot2)
      library(BayesianTools)
      library(data.table)
      library(limma)
    })
    
    clusterExport(cl, varlist)
    
    return(cl)
  }
  
  
  # Stop the cluster
  finalize_cluster <- function(cl) {
    stopCluster(cl)
  }
  
  
  
  # if there is cis and trans uricanic acid in the data then sum it up and give back the dataframe
  sum_up_urocanic_acid <- function(targeted_exp_data) {
    # If Urocanic Acid is used in list than need to sum up cis and trans forms
    if (all(c('trans-Urocanic acid', 'cis-Urocanic acid') %in% targeted_exp_data$Molecule.Name)) {
      # Sum up cis and trans
      summed_row <- targeted_exp_data %>%
        filter(Molecule.Name %in% c('trans-Urocanic acid', 'cis-Urocanic acid')) %>%
        summarise(across(5:(ncol(targeted_exp_data)-1), sum)) %>%
        mutate(Molecule.Name = 'Urocanic acid sum')
      
      # Extract the metadata from trans-Urocanic acid
      trans_metadata <- targeted_exp_data %>%
        filter(Molecule.Name == 'trans-Urocanic acid') %>%
        select(id, rt, mz, charge) 
      
      # Combine the metadata and the summed values
      summed_row <- bind_cols(trans_metadata, summed_row)
      
      # Remove the original rows and add the new row
      targeted_exp_data <- targeted_exp_data %>%
        filter(!Molecule.Name %in% c('trans-Urocanic acid', 'cis-Urocanic acid')) %>%
        bind_rows(summed_row)
      
      row.names(targeted_exp_data) <- targeted_exp_data$id
    }
    return(targeted_exp_data)
  }
  
  
  # a function to filter the list of targeted normalization molecules for certain choosen ones
  filter_list_of_norm_molecules <- function(norm_feature_list, norm_molecules) {
    
    norm_feature_list <- norm_feature_list %>% filter(Molecule.Name %in% norm_molecules)
    
    return(norm_feature_list)
  }
  
  
  
  # Function to evaluate all combinations of normalization molecules and normalization methods to be able to choose the best
  # do this by fitting to a known biological kinetik curve like caffeine
  calculate_all_different_normalization_combinations_check_fits_for_biological_kinetics <- 
    function(targeted_experiment_data, norm_feature_list, resultsdir, add_info_file, molecule, Group) {
      
    #norm_methods <- c("non-weighted", "Sqrt", "Rank2", "Median", "Med.Std.dev", "Mean")
    norm_methods <- c("non-weighted", "Rank2", "Median", "Med.Std.dev")
      
    norm_molecules <- c("Citrulline", "Glutamic acid", "Hypoxanthine", "Xanthine", "Creatinine", "DL-Tyrosine", "Iso-Leucine1", "DL-Phenylalanine", "DL-Tryptophan", "Adenosine", "Uric Acid" )
    #norm_molecules <- c("Citrulline", "Glutamic acid", "Hypoxanthine", "Xanthine", "Creatinine", "DL-Tyrosine", "Iso-Leucine1", "DL-Phenylalanine", "DL-Tryptophan", "Adenosine", "Uric Acid" )
    
    norm_feature_list <- sum_up_urocanic_acid(norm_feature_list)
    norm_feature_list_pqn <- sum_up_urocanic_acid(norm_feature_list_pqn)
    norm_feature_list_vsn <- sum_up_urocanic_acid(norm_feature_list_vsn)
    
    # norm_feature_list <- filter_list_of_norm_molecules(norm_feature_list, norm_molecules)
    # norm_feature_list_pqn <- filter_list_of_norm_molecules(norm_feature_list_pqn, norm_molecules)
    # norm_feature_list_vsn <- filter_list_of_norm_molecules(norm_feature_list_vsn, norm_molecules)
    
    
    
    combinations <- lapply(1:length(norm_feature_list$Molecule.Name), function(i) {
      combn(norm_feature_list$Molecule.Name, i, simplify = FALSE)
    }) %>% unlist(recursive = FALSE)
    
    combinations <- lapply(1:length(norm_feature_list_pqn$Molecule.Name), function(i) {
      combn(norm_feature_list_pqn$Molecule.Name, i, simplify = FALSE)
    }) %>% unlist(recursive = FALSE)
    
    combinations <- lapply(1:length(norm_feature_list_vsn$Molecule.Name), function(i) {
      combn(norm_feature_list_vsn$Molecule.Name, i, simplify = FALSE)
    }) %>% unlist(recursive = FALSE)
    
    full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
    full_prep_data_pqn <- prepare_data_for_plot(targeted_experiment_data_pqn, add_info_file)
    full_prep_data_vsn <- prepare_data_for_plot(targeted_experiment_data_vsn, add_info_file)
    
    
    # want to look at H2O and not IPA for now
    full_prep_data <- full_prep_data[full_prep_data$Solvent == "H2O",]
    full_prep_data_pqn <- full_prep_data_pqn[full_prep_data_pqn$Solvent == "H2O",]
    full_prep_data_vsn <- full_prep_data_vsn[full_prep_data_vsn$Solvent == "H2O",]
    
    # want to look only at certain paper
    full_prep_data <- full_prep_data[full_prep_data$paper == "Chrom",]
    full_prep_data_pqn <- full_prep_data_pqn[full_prep_data_pqn$paper == "Chrom",]
    full_prep_data_vsn <- full_prep_data_vsn[full_prep_data_vsn$paper == "Chrom",]
    
    # remove outlier
    #outlier <- c("230602_SimonBAProjekt_041_P2_IPA_t0_li", "230602_SimonBAProjekt_054_P2_IPA_t3_re", "230606_SimonBAProjekt_021_P5_H2O_t0_li", "230606_SimonBAProjekt_027_P5_H2O_t6_li", "230606_SimonBAProjekt_043_P6_IPA_t2_li")
    #outlier <- c("240415_PhilKoff_Chrom0411_t0", "240415_PhilKoff_Chrom0411_t4")
    outlier <- c("200922_finger_sweat_capsule_cg_11h", "200922_finger_sweat_capsule_cg_7h", "200922_finger_sweat_capsule_cg_8h", 
                 "200922_finger_sweat_capsule_cg_9h", "200922_finger_sweat_capsule_hehe_12h", "200922_finger_sweat_capsule_mago_24h", 
                 "CC_D10_26_h", "CC_D12_60_min", "CC_D12__60_min", "CC_D13_25_h", "CC_D13__25_h", "CC_D16_24h",                     
                 "CC_D17_12h", "CC_D19_11h", "CC_D19_7h", "CC_D19_8h", "CC_D19_9h", "CC_D2_3_h", "CC_D2_gm_3_h", 
                 "CC_D3_8_h", "CC_D3_ab_8_h", "CC_D9_0_min", "CC_D9_4_h")
    full_prep_data <- full_prep_data[!full_prep_data$Sample %in% outlier, ]
    full_prep_data <- full_prep_data[!is.na(full_prep_data$Donor),]
    norm_feature_list <- norm_feature_list %>%
      filter(Molecule.Name %in% norm_feature_list_pqn$Molecule.Name)
    norm_feature_list <- norm_feature_list %>%
      select(-one_of(outlier))
    
    unique_donors <- unique(full_prep_data[[Group]])
  
    combinations_all <- combinations # to save the combinations
    #combinations <- combinations_all # to save the combinations
    
    #combinations <- c(combinations[1:333], combinations[100001:100333], combinations[10000000:10000333]) #reduce the list for testing and checking how long it takes
    #combinations <- c(combinations[1:3], combinations[100001:100003], combinations[10000000:10000003]) #reduce the list for testing and checking how long it takes
    
    
    
    process_donor <- function(donor) { 
      # calculate the pure raw data
        raw_results <- data.frame(
                                  Donor = character(),          # Initialize Donor as character
                                  Res.stand.err = numeric(),    # Initialize Res.stand.err as numeric
                                  ka = numeric(),               # Initialize ka as numeric
                                  kel = numeric(),              # Initialize kel as numeric
                                  Vd = numeric(),               # Initialize Vd as numeric
                                  Fit.Data = I(list()),         # Initialize Fit.Data as list
                                  Fit.Curve = I(list())         # Initialize Fit.Curve as list
                                 )
        
        donor_data <- full_prep_data[full_prep_data[[Group]] == donor & full_prep_data$Molecule.Name == molecule, ]
        donor_data$scaled_conc <- donor_data$Area - min(donor_data$Area, na.rm = TRUE)
        donor_data$scaled_conc <- donor_data$scaled_conc / max(donor_data$scaled_conc[donor_data$time > 0 & donor_data$time < 100])
        
        pure_donor_data_for_norm <- subset(targeted_experiment_data, Molecule.Name == molecule)
        pure_donor_data_for_norm <- pure_donor_data_for_norm[colnames(pure_donor_data_for_norm) %in% donor_data$Sample]
        
        we <- rep(1, nrow(donor_data))
        
        # Fit raw data first
        fit <- try(nlsLM(scaled_conc ~ bateman(time, ka, kel, Vd), data = donor_data,
                         start = list(ka = 0.03, kel = 0.02, Vd = 30),
                         lower = c(ka = 0.01, kel = 0.001, Vd = 20),
                         upper = c(ka = 0.05, kel = 0.05, Vd = 80),
                         weights = we,
                         control = list(maxiter = 100)), silent = TRUE)
        
        if (class(fit) != "try-error") {
          summary_model <- summary(fit)
  
          fit_curve <- data.frame(time = seq(0, max(donor_data$time)))
          fit_curve$scaled_conc <- predict(fit, list(time = fit_curve$time))
          
          raw_results <- bind_rows(raw_results, data.frame(
                                                        Donor = donor, 
                                                        Res.stand.err = as.numeric(summary_model$sigma), 
                                                        ka = as.numeric(coef(fit)['ka']),
                                                        kel = as.numeric(coef(fit)['kel']),
                                                        Vd = as.numeric(coef(fit)['Vd']),
                                                        Fit.Data = I(list(donor_data)),
                                                        Fit.Curve = I(list(fit_curve))
                                                      ))
        }
      
      # calcualte the pqn data
        pqn_results <- data.frame(
          Donor = character(),          # Initialize Donor as character
          Res.stand.err = numeric(),    # Initialize Res.stand.err as numeric
          ka = numeric(),               # Initialize ka as numeric
          kel = numeric(),              # Initialize kel as numeric
          Vd = numeric(),               # Initialize Vd as numeric
          Fit.Data = I(list()),         # Initialize Fit.Data as list
          Fit.Curve = I(list())         # Initialize Fit.Curve as list
        )
        
        donor_data_pqn <- full_prep_data_pqn[full_prep_data_pqn[[Group]] == donor & full_prep_data_pqn$Molecule.Name == molecule, ]
        donor_data_pqn$scaled_conc <- donor_data_pqn$Area - min(donor_data_pqn$Area, na.rm = TRUE)
        donor_data_pqn$scaled_conc <- donor_data_pqn$scaled_conc / max(donor_data_pqn$scaled_conc[donor_data_pqn$time > 0 & donor_data_pqn$time <= 120], na.rm = TRUE)
        
        pure_donor_data_for_norm_pqn <- subset(targeted_experiment_data_pqn, Molecule.Name == molecule)
        pure_donor_data_for_norm_pqn <- pure_donor_data_for_norm_pqn[colnames(pure_donor_data_for_norm_pqn) %in% donor_data_pqn$Sample]
        
        we <- rep(1, nrow(donor_data_pqn))
        
        # Fit raw data first
        fit_pqn <- try(nlsLM(scaled_conc ~ bateman(time, ka, kel, Vd), data = donor_data_pqn,
                         start = list(ka = 0.03, kel = 0.02, Vd = 30),
                         lower = c(ka = 0.01, kel = 0.001, Vd = 20),
                         upper = c(ka = 0.05, kel = 0.05, Vd = 80),
                         weights = we,
                         control = list(maxiter = 100)), silent = TRUE)
        
        if (class(fit_pqn) != "try-error") {
          summary_model_pqn <- summary(fit_pqn)
          
          fit_curve_pqn <- data.frame(time = seq(0, max(donor_data_pqn$time)))
          fit_curve_pqn$scaled_conc <- predict(fit_pqn, list(time = fit_curve_pqn$time))
          
          pqn_results <- bind_rows(pqn_results, data.frame(
            Donor = donor, 
            Res.stand.err = as.numeric(summary_model_pqn$sigma), 
            ka = as.numeric(coef(fit_pqn)['ka']),
            kel = as.numeric(coef(fit_pqn)['kel']),
            Vd = as.numeric(coef(fit_pqn)['Vd']),
            Fit.Data = I(list(donor_data_pqn)),
            Fit.Curve = I(list(fit_curve_pqn))
          ))
        }
        
      # calculate the vsn data
        vsn_results <- data.frame(
          Donor = character(),          # Initialize Donor as character
          Res.stand.err = numeric(),    # Initialize Res.stand.err as numeric
          ka = numeric(),               # Initialize ka as numeric
          kel = numeric(),              # Initialize kel as numeric
          Vd = numeric(),               # Initialize Vd as numeric
          Fit.Data = I(list()),         # Initialize Fit.Data as list
          Fit.Curve = I(list())         # Initialize Fit.Curve as list
        )
        
        donor_data_vsn <- full_prep_data_vsn[full_prep_data_vsn[[Group]] == donor & full_prep_data_vsn$Molecule.Name == molecule, ]
        donor_data_vsn$scaled_conc <- donor_data_vsn$Area - min(donor_data_vsn$Area, na.rm = TRUE)
        donor_data_vsn$scaled_conc <- donor_data_vsn$scaled_conc / max(donor_data_vsn$scaled_conc[donor_data_vsn$time > 0 & donor_data_vsn$time <= 120], na.rm = TRUE)
        
        pure_donor_data_for_norm_vsn <- subset(targeted_experiment_data_vsn, Molecule.Name == molecule)
        pure_donor_data_for_norm_vsn <- pure_donor_data_for_norm_vsn[colnames(pure_donor_data_for_norm_vsn) %in% donor_data_vsn$Sample]
        
        we <- rep(1, nrow(donor_data_vsn))
        
        # Fit raw data first
        fit_vsn <- try(nlsLM(scaled_conc ~ bateman(time, ka, kel, Vd), data = donor_data_vsn,
                             start = list(ka = 0.03, kel = 0.02, Vd = 30),
                             lower = c(ka = 0.01, kel = 0.001, Vd = 20),
                             upper = c(ka = 0.05, kel = 0.05, Vd = 80),
                             weights = we,
                             control = list(maxiter = 100)), silent = TRUE)
        
        if (class(fit_vsn) != "try-error") {
          summary_model_vsn <- summary(fit_vsn)
          
          fit_curve_vsn <- data.frame(time = seq(0, max(donor_data_vsn$time)))
          fit_curve_vsn$scaled_conc <- predict(fit_vsn, list(time = fit_curve_vsn$time))
          
          vsn_results <- bind_rows(vsn_results, data.frame(
            Donor = donor, 
            Res.stand.err = as.numeric(summary_model_vsn$sigma), 
            ka = as.numeric(coef(fit_vsn)['ka']),
            kel = as.numeric(coef(fit_vsn)['kel']),
            Vd = as.numeric(coef(fit_vsn)['Vd']),
            Fit.Data = I(list(donor_data_vsn)),
            Fit.Curve = I(list(fit_curve_vsn))
          ))
        }
      
        
      
      # here can change what data (raw, pwn, vsn) shall be used for normalization
        # pure_donor_data_for_norm <- pure_donor_data_for_norm_pqn
        # norm_feature_list <- norm_feature_list_pqn
        # donor_data <- donor_data_pqn
        # summary_model <- summary_model_pqn
        # fit <- fit_pqn
      
      calculate_norm_combinations_chunks <- function(comb_chunk, start_index) {
        results <- data.frame(
                              Donor = character(),          # Initialize Donor as character
                              Norm.ID = integer(),          # Initialize Norm.ID as integer
                              Norm.Molecule = character(),  # Initialize Norm.Molecule as character
                              Norm.Method = character(),    # Initialize Norm.Method as character
                              Res.stand.err = numeric(),    # Initialize Res.stand.err as numeric
                              Res.stand.err.norm = numeric(), # Initialize Res.stand.err.norm as numeric
                              ka = numeric(),               # Initialize ka as numeric
                              kel = numeric(),              # Initialize kel as numeric
                              Vd = numeric(),               # Initialize Vd as numeric
                              ka.norm = numeric(),          # Initialize ka.norm as numeric
                              kel.norm = numeric(),         # Initialize kel.norm as numeric
                              Vd.norm = numeric(),          # Initialize Vd.norm as numeric
                              Fit.Data = I(list()),         # Initialize Fit.Data as list
                              Fit.Curve = I(list())         # Initialize Fit.Curve as list
                             )

        for (i in seq_along(comb_chunk)) {
          comb <- comb_chunk[[i]]
          norm_id <- start_index + i - 1

          for (method in norm_methods) {
            if (method == "Med.Std.dev" & length(comb) < 3) next  # need this as we get problem with Std.dev otherwise

            norm_factors <- calculate_norm_factors(pure_donor_data_for_norm, method, comb, norm_feature_list)
            donor_data_norm <- apply_normalization(donor_data, norm_factors)

            donor_data_norm$scaled_conc_norm <- donor_data_norm$Area_norm - min(donor_data_norm$Area_norm, na.rm = TRUE)
            donor_data_norm$scaled_conc_norm <- donor_data_norm$scaled_conc_norm / max(donor_data_norm$scaled_conc_norm[donor_data_norm$time > 0 & donor_data_norm$time < 100])

            if (is.na(mean(donor_data_norm$scaled_conc_norm)) | median(donor_data_norm$scaled_conc_norm) < 0.05) next # sometimes molecule possibly not found, normalization doesnt work accordincly could still get good RSE score, want to remove these
            if (any(is.infinite(donor_data_norm$scaled_conc_norm))) next # sometimes norm. molecule possibly not found at one time point div. with 0 gets Inf, remove

            fit_norm <- try(nlsLM(scaled_conc_norm ~ bateman(time, ka, kel, Vd), data = donor_data_norm,
                                  start = list(ka = 0.03, kel = 0.02, Vd = 30),
                                  lower = c(ka = 0.01, kel = 0.001, Vd = 20),
                                  upper = c(ka = 0.05, kel = 0.05, Vd = 80),
                                  weights = we,
                                  control = list(maxiter = 100)), silent = TRUE)

            if (class(fit_norm) != "try-error") {
              summary_model_norm <- summary(fit_norm)

              fit_curve_norm <- data.frame(time = seq(0, max(donor_data_norm$time)))
              fit_curve_norm$scaled_conc_norm <- predict(fit_norm, list(time = fit_curve_norm$time))

              results  <- bind_rows(results , data.frame(
                                                          Donor = donor,
                                                          Norm.ID = norm_id,
                                                          Norm.Molecule = paste(comb, collapse = "+"),
                                                          Norm.Method = method,
                                                          Res.stand.err = as.numeric(summary_model$sigma),
                                                          Res.stand.err.norm = as.numeric(summary_model_norm$sigma),
                                                          ka = as.numeric(coef(fit)['ka']),
                                                          kel = as.numeric(coef(fit)['kel']),
                                                          Vd = as.numeric(coef(fit)['Vd']),
                                                          ka.norm = as.numeric(coef(fit_norm)['ka']),
                                                          kel.norm = as.numeric(coef(fit_norm)['kel']),
                                                          Vd.norm = as.numeric(coef(fit_norm)['Vd']),
                                                          Fit.Data = I(list(donor_data_norm)),
                                                          Fit.Curve = I(list(fit_curve_norm))
                                                        ))

            }
          }
        }

        return(results)
      }

      # Parallel processing of donors and combinations
      clusterExport(cl, list('norm_methods', "process_donor", "full_prep_data", 'targeted_experiment_data',
                             'norm_feature_list', 'molecule', 'combinations', 'donor_data',
                             'pure_donor_data_for_norm', 'donor', 'we', 'fit', 'summary_model',
                             'calculate_norm_combinations_chunks'), envir = environment())

      # Split combinations into chunks of max 1 million each as large chunks make problmes with cluster
      split_combinations <- function(combinations, max_size = 1e6) {
        split(combinations, ceiling(seq_along(combinations) / max_size))
      }

      results_listed <- list()
      i = 1

      for(combinations_split in split_combinations(combinations, max_size = 1e6)) {
        combinations_chunks <- split_chunks(combinations_split, detectCores() - 1)

        clusterExport(cl, list('combinations_chunks'), envir = environment())

        #Parallel processing of donors
        results_list <- parLapply(cl, seq_along(combinations_chunks), function(i)
          calculate_norm_combinations_chunks(combinations_chunks[[i]], (i - 1) * length(combinations_chunks[[i]]) + 1))

        results_listed[[i]] <- bind_rows(results_list)
        i <- i + 1
      }

      results <- bind_rows(results_listed)

      results_filtered <- results %>%
        filter(Res.stand.err.norm < 1.5 * Res.stand.err)
      
      # Collect and return raw_results and results for this donor
      donor_results <- list(raw_results = raw_results, pqn_results = pqn_results, vsn_results = vsn_results, results = results, results_filtered = results_filtered)
      print(paste0("Donor: ", donor, " done!"))
      return(donor_results)
    }
    
    
    # Start the cluster for parallelization of next function
    cl <- initialize_cluster(60, varlist = c('prepare_data_for_plot', 'calculate_norm_factors', 'apply_normalization', 
                                             'bateman', 'scale_concentration')) 
    
    start_time <- Sys.time()
    # run the code and check all combinations in the cluster in parallel
    all_results_list <- lapply(unique_donors, process_donor)
    print("FINISHED!")  
    
    end_time <- Sys.time()
    end_time - start_time
    
    # stop cluster after parallelizations is finished
    finalize_cluster(cl) 
    
    # Combine results from all donors into a single data frame
    
    final_raw_results <- bind_rows(lapply(all_results_list, function(x) x$raw_results))
    final_raw_pqn_results <- bind_rows(lapply(all_results_list, function(x) x$pqn_results))
    final_raw_vsn_results <- bind_rows(lapply(all_results_list, function(x) x$vsn_results))
    final_results <- bind_rows(lapply(all_results_list, function(x) x$results))
    end_time1 <- Sys.time()
    end_time1 - start_time
    
    # filter the results (remove rows that have values at the boarder of kel or Vd and ones where kel is bigger then kel and keep only one row for each combination)
    final_results_filtered <- final_results %>% #select(-Fit.Data, -Fit.Curve) %>%
                                                #filter(!(kel.norm %in% c(0.01, 0.03))) %>% 
                                                #filter(!(Vd.norm %in% c(20, 45)))  %>% 
                                                filter(Res.stand.err.norm < Res.stand.err *0.5) #%>%
                                                #distinct(Donor, Norm.ID, .keep_all = TRUE)
    
    combinations_filtered <- as.list(unique(final_results_filtered$Norm.Molecule))
    combinations_filtered <- lapply(combinations_filtered, function(x) {
      strsplit(x, "\\+")[[1]]
    })
    combinations <- combinations_filtered
    
    # count how often each molecule occurs in the different good combinations
    molecule_counts  <- final_results_filtered %>%
      mutate(Norm.Molecule = strsplit(Norm.Molecule, "\\+")) %>%
      unnest(Norm.Molecule) %>%
      group_by(Norm.Molecule) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    weighting_counts  <- final_results_filtered %>%
      group_by(Norm.Method) %>% 
      summarise(Count = n(), .groups = 'drop') %>%  
      arrange(desc(Count))  # Sort by Count in descending order
    
    write.csv(molecule_counts, file = paste0(resultsdir, "molecule_counts_normalization_combinations_filtered.csv"))
    write.csv(weighting_counts, file = paste0(resultsdir, "weighting_counts_normalization_combinations_filtered.csv"))
    write.csv(molecule_counts, file = paste0(resultsdir, "molecule_counts_normalization_combinations_filtered_Res0.5.csv"))
    write.csv(weighting_counts, file = paste0(resultsdir, "weighting_counts_normalization_combinations_filtered_Res0.5.csv"))
    
    
    # make outputs (write to file and make plots)
    ranked_results <- rank_combinations(final_results)
    fwrite(ranked_results[-c(13:14)], file = paste0(resultsdir, "top_20_normalization_results.csv"))
    saveRDS(ranked_results, file = paste0(resultsdir, "top_20_normalization_results_with_curve_data.rds")) # can be read in with readRDS
    end_time2 <- Sys.time()
    end_time2 - start_time
    
    plot_top_20_curves(ranked_results, final_raw_results, resultsdir, molecule, final_raw_pqn_results, final_raw_vsn_results)
    plot_top_20_residuals(ranked_results, final_results, final_raw_results, resultsdir, molecule, final_raw_pqn_results, final_raw_vsn_results)
    plot_histograms_and_scatter_density(final_results, final_raw_vsn_results, resultsdir, molecule)
    
    fwrite(final_results[-c(13:14)], file = paste0(resultsdir, "all_results_of_all_normalization_combinations.csv"))
    saveRDS(final_results, file = paste0(resultsdir, "all_results_of_all_normalization_combinations_with_curve_data.rds")) # can be read in with readRDS
    fwrite(final_raw_results[-c(6:7)], file = paste0(resultsdir, "all_raw_results.csv"))
    saveRDS(final_raw_results, file = paste0(resultsdir, "all_raw_results_with_curve_data.rds")) # can be read in with readRDS
    fwrite(final_raw_pqn_results[-c(6:7)], file = paste0(resultsdir, "all_pqn_results.csv"))
    saveRDS(final_raw_pqn_results, file = paste0(resultsdir, "all_raw_pqn_results_with_curve_data.rds")) # can be read in with readRDS
    fwrite(final_raw_vsn_results[-c(6:7)], file = paste0(resultsdir, "all_raw_vsn_results.csv"))
    saveRDS(final_raw_vsn_results, file = paste0(resultsdir, "all_raw_vsn_results_with_curve_data.rds")) # can be read in with readRDS
    end_time3 <- Sys.time()
    end_time3 - start_time
    
    # save.image(paste0(resultsdir, "normalization_combinations_test_finished.RData"))
    # end_time4 <- Sys.time()
    # end_time4 - start_time
    
    print("All Done")
  }
  
  
  
  # Function to cluster data with kmeans and calculate silhouette_scores 
  get_silhouette_scores_with_kmeans_clustering <- function(pure_data_for_clustering, n_clusters = 3) {
    # Perform kmeans clustering
    kmeans_result <- kmeans(pure_data_for_clustering, centers = n_clusters) # Perform k-means clustering
    cluster_assignments <- kmeans_result$cluster # Add cluster assignments to the data
    
    silhouette_scores <- silhouette(cluster_assignments, dist(pure_data_for_clustering)) # Compute silhouette scores
    avg_silhouette_per_cluster <- aggregate(silhouette_scores[, 3],   # Calculate average silhouette score for each cluster
                                                by = list(Cluster = silhouette_scores[, 1]), 
                                                FUN = mean)
    sd_silhouette_per_cluster <- aggregate(silhouette_scores[, 3],   # Calculate average silhouette score for each cluster
                                               by = list(Cluster = silhouette_scores[, 1]), 
                                               FUN = sd)
    
    return(list(avg_silhouette_per_cluster, sd_silhouette_per_cluster, cluster_assignments))
  }
  
  
  
  # Function to cluster data with hierarchical clustering and calculate silhouette scores
  get_silhouette_scores_with_hierarchical_clustering <- function(pure_data_for_clustering, n_clusters = NULL, method = "ward.D2") {
    distance_matrix <- dist(pure_data_for_clustering) # Compute the distance matrix
    hclust_result <- hclust(distance_matrix, method = method) # Perform hierarchical clustering
    
    if (is.null(n_clusters)) {
      # Find optimal number of clusters automatically
      best_n_clusters <- 2  # Start from 2 clusters (1 cluster isn't useful for silhouette analysis)
      best_avg_silhouette <- -1 # start from worst possible score
      silhouette_scores_list <- list()
      max_clusters <- 15
      
      # Loop through different cluster numbers to find the best one
      for (k in 2:max_clusters) {
        cluster_assignments <- cutree(hclust_result, k = k)
        silhouette_scores <- silhouette(cluster_assignments, distance_matrix)
        avg_silhouette <- mean(silhouette_scores[, 3])  # Average silhouette score for all points
        
        silhouette_scores_list[[as.character(k)]] <- avg_silhouette
        
        if (avg_silhouette > best_avg_silhouette) {
          best_avg_silhouette <- avg_silhouette
          best_n_clusters <- k
        }
      }
      
      n_clusters = best_n_clusters
    }
    
    cluster_assignments <- cutree(hclust_result, k = n_clusters) # Cut the tree into n_clusters
    silhouette_scores <- silhouette(cluster_assignments, distance_matrix) # Compute silhouette scores
    avg_silhouette_per_cluster <- aggregate(silhouette_scores[, 3], # Calculate average silhouette score for each cluster
                                            by = list(Cluster = silhouette_scores[, 1]), 
                                            FUN = mean)
    sd_silhouette_per_cluster <- aggregate(silhouette_scores[, 3], # Calculate standard deviation of silhouette scores for each cluster
                                           by = list(Cluster = silhouette_scores[, 1]), 
                                           FUN = sd)
    
    return(list(avg_silhouette_per_cluster, sd_silhouette_per_cluster, cluster_assignments))
  }
  
  
  
  # Function to generate plots for a given cluster assignment
  plot_clusters <- function(cluster_data, data, Group_list = None) {
    avg_silhouette_per_cluster <- cluster_data[[1]]
    sd_silhouette_per_cluster <- cluster_data[[2]]
    cluster_assignments <- cluster_data[[3]]
    n_clusters <- max(cluster_assignments)
    plots <- list()
    plot_index <- 1  # Index for storing plots
    
    # Generate unique group combinations only if Group_list is provided
    if (!is.null(Group_list) && length(Group_list) > 0) {
      unique_groups <- lapply(Group_list, function(group) unique(data[[group]]))
      names(unique_groups) <- Group_list
      group_combinations <- expand.grid(unique_groups, KEEP.OUT.ATTRS = FALSE)
    } else {
      group_combinations <- data.frame(dummy = 1)  # Create a single-row data frame
    }
    
    for (i in 1:n_clusters) {
      if ('time' %in% colnames(data)) {
        cluster_plot_data <- data %>%
          filter(id %in% names(cluster_assignments)[cluster_assignments == i]) %>%
          group_by(id, time, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%  # Group by id and time (and other relevant grouping vars)
          summarise(Area = mean(Area, na.rm = TRUE), .groups = "drop")  # Compute mean of Area
      } else {
        cluster_plot_data <- data %>%
          filter(id %in% names(cluster_assignments)[cluster_assignments == i]) %>%
          group_by(id, Timepoint, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%  # Group by id and time (and other relevant grouping vars)
          summarise(Area = mean(Area, na.rm = TRUE), .groups = "drop")  # Compute mean of Area
      }
      
      # Loop through each group combination
      for (j in 1:nrow(group_combinations)) {
        if (!is.null(Group_list) && length(Group_list) > 0) {
          # Convert group_combinations[j, ] to a named list
          group_values <- as.list(group_combinations[j, ])
          
          # Subset data based on the current combination of group values
          subset_data <- cluster_plot_data[
            Reduce(`&`, Map(function(col, val) cluster_plot_data[[col]] == val, names(group_values), group_values)), 
          ]
        } else {
          subset_data <- cluster_plot_data  # If no grouping variables, use full cluster data
        }
        
        # Skip empty subsets
        if (nrow(subset_data) == 0) next 
        
        group_labels <- sapply(names(group_combinations), function(col) {
          levels(group_combinations[[col]])[group_combinations[j, col]]
        })
        
        
        # Create the plot
        if ('time' %in% colnames(data)) {
          p <- ggplot(subset_data, aes(
            x = time, 
            y = Area, 
            group = id
          )) +
            geom_line() +
            labs(
              title = paste0("Cluster ", i, " | ", 
                             paste(group_labels, collapse = ", "), 
                             "\nAvg.Sil.Score=", round(avg_silhouette_per_cluster[i, 2], 3),
                             "\nStd.Dev.Sil.Score=", round(sd_silhouette_per_cluster[i, 2], 3)), 
              x = "Time", 
              y = "Area"
            ) +
            theme_minimal()
        } else {
          p <- ggplot(subset_data, aes(
            x = Timepoint, 
            y = Area, 
            group = id
          )) +
            geom_line() +
            labs(
              title = paste0("Cluster ", i, " | ", 
                             paste(names(group_labels), group_labels, collapse = ", "), 
                             "\nAvg.Sil.Score=", round(avg_silhouette_per_cluster[i, 2], 3),
                             "\nStd.Dev.Sil.Score=", round(sd_silhouette_per_cluster[i, 2], 3)), 
              x = "Timepoint", 
              y = "Area"
            ) +
            theme_minimal()
        }
        
        # Store the plot
        plots[[plot_index]] <- p
        plot_index <- plot_index + 1
      }
    }
    
    return(plots)
  }
  

  
  # prepare data for limma test
  prepare_data_for_limma <- function(all_needed_features, add_info, Group_list = NULL) {
    # Reshape all_needed_features to have samples as rows and features as columns
    if ("rtime_group" %in% colnames(all_needed_features)) {
      sample_names <- colnames(all_needed_features[10:length(all_needed_features)])
      reshaped_data <- pivot_longer(all_needed_features, 
                                    cols = -c(id, rt, mz, rtime_group, feature_group, Annotation, Formula, Confidence, SMILES),
                                    names_to = "variable", 
                                    values_to = "value")
    } else {
      sample_names <- colnames(all_needed_features[5:(length(all_needed_features)-1)])
      reshaped_data <- pivot_longer(all_needed_features,
                                    cols = -c(id, rt, mz, charge, Molecule.Name),
                                    names_to = "variable",
                                    values_to = "value")
    }
    
    reshaped_data <- reshape2::dcast(reshaped_data, variable ~ id, value.var = "value")
    colnames(reshaped_data)[1] <- "Sample"
    rownames(reshaped_data) <- reshaped_data[,1]
    pure_data <- reshaped_data[,-1]
    
    # Merge reshaped data with add_info
    merged_data <- merge(reshaped_data, add_info, by.x = "Sample", by.y = "Sample")
    
    merged_data[, 2:(length(merged_data) - (ncol(add_info) -1))] <- as.data.frame(apply(merged_data[, 2:(length(merged_data) - (ncol(add_info) -1))], 2, replace_na_and_Inf))
    pure_data <- as.data.frame(apply(pure_data, 2, replace_na_and_Inf))
    pure_data <- log2(pure_data)
    row.names(pure_data) <- merged_data$Sample 
    
    # Merge pure_data with metadata
    full_data <- cbind(merged_data[, c("Sample", "Donor", "Timepoint", unlist(Group_list))], pure_data)
    
    return(full_data)
  }
  
  
  # Get Limma Fit for checking a significant change over the time
  get_limma_fit_for_signifint_changes_over_time <- function(full_data) {
    # Define the design matrix
    design <- model.matrix(~ 0 + factor(full_data$Timepoint), data = full_data)
    colnames(design) <- paste0("Timepoint", unique(full_data$Timepoint))
    
    # Define contrasts for each consecutive timepoint pair (e.g., Timepoint 1 vs. Timepoint 2, Timepoint 2 vs. Timepoint 3, etc.)
    timepoints <- unique(full_data$Timepoint)
    contrast_list <- list()
    
    # Generate contrasts for each consecutive timepoint pair
    for (i in 1:(length(timepoints) - 1)) {
      cond1 <- paste0("Timepoint", timepoints[i])
      cond2 <- paste0("Timepoint", timepoints[i+1])
      
      contrast_list[[paste0(cond2, "_vs_", cond1)]] <- paste0(cond2, " - ", cond1)
    }
    
    # Create the contrast matrix from the list
    contrast_matrix <- makeContrasts(contrasts = unlist(contrast_list), levels = design)
    
    # Prepare the expression matrix (excluding non-numeric columns)
    nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(full_data)))
    expression_matrix <- t(full_data[, -c(1:nonnumeric_cols)])
    
    # Define the blocking factor for repeated measures (subjects)
    blocking_factor <- factor(full_data$Donor)
    
    # Estimate the correlation within subjects using duplicateCorrelation
    corfit <- duplicateCorrelation(expression_matrix, design, block = blocking_factor)
    
    # Fit the linear model
    if (!is.na(corfit$consensus)) {
      fit <- lmFit(expression_matrix, design, block = blocking_factor, correlation = corfit$consensus)
    } else {
      fit <- lmFit(expression_matrix, design)
    }
    
    # Apply the contrasts
    fit <- contrasts.fit(fit, contrast_matrix)
    
    return(topTable(eBayes(fit), number = Inf))
  }
  
  
    
  # Function to get fold changes for time series data
  get_fold_change_dataframe <- function(data, add_info, Group_list) {
    
    full_data <- prepare_data_for_limma(data, add_info, Group_list)
    
    # Ensure the data is sorted by Donor, Timepoint and all Group_list variables, ensuring correct grouping
    sort_cols <- c("Donor", unlist(Group_list), "Timepoint") # Ensure Donor is included
    sort_cols <- sort_cols[sort_cols %in% colnames(full_data)] # Keep only existing columns
    full_data <- full_data[do.call(order, full_data[sort_cols]), ]
    
    # Get the list of unique time points in ascending order
    time_points <- sort(unique(full_data$Timepoint))
    
    # Get feature columns (excluding metadata like Sample, Donor, Timepoint, etc.)
    feature_cols <- setdiff(colnames(full_data), c("Sample", "Donor", "Timepoint", unlist(Group_list)))
    
    # Initialize an empty list to store fold changes
    fold_change_list <- list()
    
    # Loop over consecutive time points
    for (i in seq_along(time_points)[-1]) {
      time_prev <- time_points[i - 1]
      time_curr <- time_points[i]
      
      # Subset data for the two time points
      data_prev <- full_data[full_data$Timepoint == time_prev, ]
      data_curr <- full_data[full_data$Timepoint == time_curr, ]
      
      # Merge using all grouping variables
      merge_cols <- c("Donor", unlist(Group_list)) # Ensure grouping is maintained
      merge_cols <- merge_cols[merge_cols %in% colnames(full_data)] # Keep only valid columns
      
      merged_data <- merge(
        data_prev[, c(merge_cols, feature_cols)], 
        data_curr[, c(merge_cols, feature_cols)], 
        by = merge_cols, suffixes = c("_prev", "_curr")
      )
      
      # Ensure all values are positive for log2 transformation
      min_value <- min(merged_data[-seq_along(merge_cols)])
      if (min_value <= 0) {
        merged_data[-seq_along(merge_cols)] <- merged_data[-seq_along(merge_cols)] + abs(min_value) + 0.1
      }
      
      # Compute log2 fold change for each feature
      fold_changes <- log2(merged_data[, paste0(feature_cols, "_curr")] / merged_data[, paste0(feature_cols, "_prev")])
      
      # Generate dynamic group identifiers
      group_id <- apply(merged_data[, c(unlist(Group_list)), drop = FALSE], 1, paste, collapse = "_")
      
      # Store results with stratification
      for (grp in unique(group_id)) {
        fold_change_list[[paste0(grp, "_T", time_prev, "_to_T", time_curr)]] <- 
          colMeans(fold_changes[group_id == grp, , drop = FALSE], na.rm = TRUE)
      }
    }
    
    # Combine into a final dataframe
    fold_change_df <- do.call(cbind, fold_change_list)
    rownames(fold_change_df) <- feature_cols
    
    return(fold_change_df)
  }
  
    
    
  # Function to make plots for various clustering of the data
  get_timeseries_plot_for_various_clustering <- function(clustered_data) {
    full_prep_cluster_data <- prepare_data_for_plot(clustered_data, add_info_file)
    fold_change_df <- get_fold_change_dataframe(clustered_data, add_info, Group_list = list('Intervention'))
    
    # # Define the PDF output file
    # pdf(paste0(resultsdir,"cluster_hierarchical_plots.pdf"), width = 25, height = 12)
    # 
    # # Loop through different cluster numbers
    # for (n_clusters in 2:15) {
    #   #cluster_data <- get_silhouette_scores_with_kmeans_clustering(fold_change_df, n_clusters = n_clusters)
    #   cluster_data <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df, n_clusters = n_clusters)
    #   
    #   # Generate plots for the current clustering
    #   plots <- plot_clusters(cluster_data, full_prep_cluster_data)
    #   
    #   # New page for each clustering
    #   grid.arrange(grobs = plots, ncol = 5, nrow = 3)
    # }
    # 
    # # Close the PDF device
    # dev.off()
    
    
    # Define the PDF output file
    pdf(paste0(resultsdir,"cluster_hierarchical_automatic_plots.pdf"), width = 25, height = 12)
    
    cluster_data <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df) # automatically find best cluster number
    full_prep_cluster_data <- prepare_data_for_plot(clustered_data, add_info_file)
    plots <- plot_clusters(cluster_data, full_prep_cluster_data) # Generate plots for the current clustering
    grid.arrange(grobs = plots, ncol = 5, nrow = 3)
    
    # Close the PDF device
    dev.off()

  } 
  
  
  
  # Function to compute standard deviation of each column from a list of dataframes
  # and exclude certain cells with a identical list with values for filtering
  compute_sd_after_filtering <- function(corr_list, area_list, threshold = 1000) {
    sds <- lapply(seq_along(corr_list), function(i) {
      df_corr <- corr_list[[i]]
      df_area <- area_list[[i]]
      
      # Ensure dimensions match
      if (!all(dim(df_corr) == dim(df_area))) {
        stop("Dimension mismatch between correction and area dataframes at index ", i)
      }
      
      # Set values in df_corr to NA where df_area < threshold
      df_corr[df_area < threshold] <- NA
      
      # Compute standard deviation for each column
      apply(df_corr, 2, sd, na.rm = TRUE)
    })
    
    # Compute the average of standard deviations across all dataframes
    mean_sd <- mean(unlist(sds), na.rm = TRUE)
    return(mean_sd)
  }
  
  
  
  # Function to get pqn normalization factors for each sample/column
  # needs a data matrix with retention time (columns are samples and rows are features)
  get_pqn_factors <- function(df, stable = FALSE) {
    if ("rtime_group" %in% colnames(df) ) {
      if (stable) {
        df <- df %>% filter(rt > 0.8) %>% filter(rt < 4)
      }
      data_matrix <- as.matrix(df[-c(1:9)])
    } else {
      if (stable) {
        df <- df %>% filter(rt > 0.8) %>% filter(rt < 4)
      }
      data_matrix <- as.matrix(df[-c(1:4, ncol(df))])
    }
    
    data_matrix[data_matrix < 0] <- 0 # Ensure non-negative values in the data (VSN requires non-negative data)
    reference_spectrum <- apply(data_matrix, 1, median, na.rm = TRUE) # Calculate the reference spectrum (median of all samples)
    pqn_factors <- apply(data_matrix, 2, function(sample) { # Compute normalization factors for each sample
      median(sample / reference_spectrum, na.rm = TRUE)
    })
    pqn_normalization_factors <- t(as.data.frame(pqn_factors))
    return(pqn_normalization_factors)
  }
  
  
  # Function to get vsn normalization factors for each sample/column
  # needs a data matrix with retention time (columns are samples and rows are features)
  get_vsn_factors <- function(df, stable = FALSE) {
    if ("rtime_group" %in% colnames(df) ) {
      if (stable) {
        df <- df %>% filter(rt > 0.8) %>% filter(rt < 4)
      }
      data_matrix <- as.matrix(df[-c(1:9)])
    } else {
      if (stable) {
        df <- df %>% filter(rt > 0.8) %>% filter(rt < 4)
      }
      data_matrix <- as.matrix(df[-c(1:4, ncol(df))])
    }
    
    data_matrix[data_matrix < 0] <- 0 # Ensure non-negative values in the data (VSN requires non-negative data)
    vsn_fit <- vsn2(data_matrix)
    vsn_normalized_data <- predict(vsn_fit, data_matrix)  # Normalize the data
    vsn_normalization_factors <- exp(colMeans(vsn_normalized_data))   # compute normalization factors (e.g., mean or scaling for each sample)
    vsn_normalization_factors <- t(vsn_normalization_factors)
    rownames(vsn_normalization_factors) <- "vsn_normalization_factors"
    return(vsn_normalization_factors)
  }
  
  
  
  # Make PCA plot
  get_PCA_plot <- function(pure_data, merged_data, Group, title = "PCA Data", resultsdir, filename) {
    pca <- prcomp(pure_data, center = TRUE, scale. = TRUE)
    plot <- autoplot(pca, data = merged_data, colour = Group, frame = TRUE, frame.type = 'norm') +
      labs(title = title) +
      theme_minimal()
    save_plot(paste0(resultsdir, filename), plot)
  }
  
  
  # Function to VSN normalize data and perform all necessary transformations for bayesOpt
  vsn_normalize_and_transform <- function(needed_data, all_data, df_ids, transitionFile, transitionFile_cluster, normalization_list, add_info_file, donor, molecule) {
    # Extract metadata and pure data
    metadata <- needed_data[, 1:9]
    pure_data <- needed_data[, -c(1:9)]
    all_needed_features_norm <- needed_data
    
    # Perform VSN normalization
    vsn_fit <- vsn2(as.matrix(pure_data))
    all_needed_features_norm[, -c(1:9)] <- predict(vsn_fit, as.matrix(pure_data))
    all_needed_features_norm[, -c(1:9)] <- all_needed_features_norm[, -c(1:9)] + (abs(min(all_needed_features_norm[, -c(1:9)], na.rm = TRUE)) + 0.1) # shift everything into positive
    
    vsn_fit <- vsn2(as.matrix(all_data))
    all_data_norm <- all_data
    all_data_norm <- predict(vsn_fit, as.matrix(all_data))
    all_data_norm <- all_data_norm + (abs(min(all_data_norm, na.rm = TRUE)) + 0.1) # shift everything into positive
    
    # Extract metadata and pure data
    pure_data_norm <- all_needed_features_norm[, -c(1:9)]
    
    # Extract targeted and clustered experiment data
    targeted_experiment_data <- extract_feature_list(all_data_norm, df_ids, transitionFile)
    clustered_data <- extract_feature_list(all_data_norm, df_ids, transitionFile_cluster)
    colnames(clustered_data) <- normalize_names(colnames(clustered_data))
    norm_feature_list <- extract_feature_list(all_data_norm, df_ids, normalization_list)
    norm_feature_list <- norm_feature_list %>% distinct(id, .keep_all = TRUE) # remove duplicates
    colnames(norm_feature_list) <- normalize_names(colnames(norm_feature_list))
    
    # Remove outlier from dataframes
    targeted_experiment_data <- targeted_experiment_data[, c(colnames(targeted_experiment_data[c(1:4)]), intersect(colnames(targeted_experiment_data[-c(1:4, ncol(targeted_experiment_data))]), 
                                                                                                                   colnames(needed_data)), colnames(targeted_experiment_data[ncol(targeted_experiment_data)]))]
    clustered_data <- clustered_data[, c(colnames(clustered_data[c(1:4)]), intersect(colnames(clustered_data[-c(1:4, ncol(clustered_data))]), 
                                                                                     colnames(needed_data)), colnames(clustered_data[ncol(clustered_data)]))]
    norm_feature_list <- norm_feature_list[, c(colnames(norm_feature_list[c(1:4)]), intersect(colnames(norm_feature_list[-c(1:4, ncol(norm_feature_list))]), 
                                                                                              colnames(needed_data)), colnames(norm_feature_list[ncol(norm_feature_list)]))]
    
    # Prepare pure data for normalization
    pure_data_for_norm <- targeted_experiment_data %>% select(-c(1:4, Molecule.Name))
    
    full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
    full_prep_data <- full_prep_data[!is.na(full_prep_data$Donor),]
    
    # Prepare data for plotting
    full_prep_data <- prepare_data_for_plot(targeted_experiment_data, add_info_file)
    full_prep_data <- full_prep_data[!is.na(full_prep_data$Donor), ] %>%
      filter(!negControl) %>%
      mutate(
        time = case_when(
          Molecule.Name %in% c("Tryptanthrin", "3,7-Dihydro-1-butyl-7-(5,6-dihydroxyhexyl)-3-methyl-1H-purine-2,6-dione") ~ time - 120,
          TRUE ~ time - 180
        )
      ) %>%
      filter(time >= 0) %>%
      mutate(Area = ifelse(time == 0, 0, Area)) # here set to 0 instead of 300 as done for the raw values, as this already scaled and not the raw Abundances
    
    # Filter donor data
    donor_data <- full_prep_data %>% filter(!!sym(Group) == donor, Molecule.Name == molecule)

    # Prepare pure donor data for normalization
    pure_donor_data_for_norm <- as.data.frame(t(donor_data$Area))
    colnames(pure_donor_data_for_norm) <- donor_data$Sample
    adjusted_colnames <- sub("^X", "", colnames(pure_donor_data_for_norm))
    pure_donor_data_for_norm <- pure_donor_data_for_norm[adjusted_colnames %in% donor_data$Sample]
    
    # Process norm_feature_list for donor samples
    adjusted_colnames <- sub("^X", "", colnames(norm_feature_list))
    keep_columns <- c(1:4, which(adjusted_colnames %in% donor_data$Sample)[match(colnames(pure_donor_data_for_norm), adjusted_colnames[which(adjusted_colnames %in% donor_data$Sample)])], ncol(norm_feature_list))
    norm_feature_list_donor <- norm_feature_list[, keep_columns]
    
    return(list(donor_data=donor_data, pure_donor_data_for_norm=pure_donor_data_for_norm, norm_feature_list_donor=norm_feature_list_donor,
                clustered_data=clustered_data, norm_feature_list=norm_feature_list, pure_data_for_norm=pure_data_for_norm, all_needed_features_norm=all_needed_features_norm))
  }
  
  
  

  # Function to evaluate all combinations of normalization molecules and normalization methods to be able to choose the best
  # do this by fitting to a known biological kinetik curve like caffeine
  calculate_normalization_combinations_check_fits_for_biological_kinetics_with_optimization_algorithm <- 
    function(all_needed_features, all_data, targeted_experiment_data, norm_feature_list, full_prep_data, add_info, clustered_data_raw, resultsdir, 
             Group, df_ids, transitionFile, transitionFile_cluster, normalization_list, bayesoptim = TRUE, run=1, testing = FALSE) {
  
      # Prepare the data
      #norm_feature_list <- sum_up_urocanic_acid(norm_feature_list)
      #targeted_experiment_data <- targeted_experiment_data[targeted_experiment_data$Molecule.Name == molecule,]
      #full_prep_data <- full_prep_data[full_prep_data$Molecule.Name == molecule,]
      unique_donors <- unique(full_prep_data[[Group]])
      unique_molecules <- unique(full_prep_data$Molecule.Name)
      
      pqn_vsn_list <- as.data.frame(rbind(get_pqn_factors(all_needed_features, TRUE), get_vsn_factors(all_needed_features, TRUE)))
      
      full_data_raw <- prepare_data_for_limma(all_needed_features, add_info)
      limma_raw <- get_limma_fit_for_signifint_changes_over_time(full_data_raw)
      significant_features_raw <- length(limma_raw$adj.P.Val[limma_raw$adj.P.Val < 0.05])
      
      # Cluster raw data
      fold_change_df_raw <- get_fold_change_dataframe(clustered_data_raw, add_info, Group_list = list('Intervention'))
      cluster_data_raw <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df_raw)
      
      
      # Parameters for Bateman function fitting
      # # parameters for 1 compartment
      # fit_params <- list(
      #   start = list(ka = 0.03, kel = 0.02, Vd = 30),
      #   lower = c(ka = 0.01, kel = 0.001, Vd = 20),
      #   upper = c(ka = 0.1, kel = 0.1, Vd = 80)
      # )
      
      # fit_params <- list(
      #   start = list(ka = 0.03, kel = 0.01),
      #   lower = c(ka = 0.001, kel = 0.001),
      #   upper = c(ka = 0.1, kel = 0.1)
      # )
      
      # paramerers for 2-compartment model perpheral no reflux
      fit_params <- list(
        start = list(ka = 1, kel = 0.005, k12 = 0.002, kel2 = 0.005),
        lower = c(ka = 0.01, kel = 0.000001, k12 = 0.0002, kel2 = 0.000001),
        upper = c(ka = 10, kel = 0.01, k12 = 0.02, kel2 = 0.1)
      )
      
      
      # Function to evaluate combinations of normalization methods and molecules
      evaluate_combination <- function(method, weights, pqn_vsn_none) {
        
        # Select molecules with weight > 0
        comb <- norm_feature_list$Molecule.Name[weights > 0.5]
        #if (length(comb) == 0) return(1e6) # Penalize if no molecules are selected
        #if (sum(weights > 0.5) < 3) total_residual_error <- ((3 - sum(weights > 0.5)) * 1e4) # Soft penality for low molecule count
        #if (sum(weights > 0.5) > 9) total_residual_error <- ((sum(weights > 0.5) - 9) * 1e4) # Soft penality for high molecule count

        total_residual_error <- 0
        total_residual_error_raw <- 0
        plot_list <- list()
        volumn_corr_list_raw <- list()
        volumn_corr_list_norm <- list()
        area_values_list <- list()
        norm_worse_raw_count <- 0
        curve_params_raw <- data.frame()
        curve_params_norm <- data.frame()
        i <- 0
        
        # Make post-PQN/VSN (if selected)
        if (!is.na(pqn_vsn_none) && (pqn_vsn_none == 'post-PQN' | pqn_vsn_none == 'post-VSN')) {
          # Calculate normalization factors for all samples (needed for post-PQN/VSN)
          pure_data_for_norm <- targeted_experiment_data[-c(1:4)] %>% select(-Molecule.Name)
          all_norm_factors <- calculate_norm_factors(pure_data_for_norm, method, comb, norm_feature_list)
          
          # Normalize all samples
          all_needed_features_norm <- all_needed_features
          all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / all_norm_factors
          
          # Get post-PQN/VSN factors
          post_pqn_vsn_list <- as.data.frame(rbind(get_pqn_factors(all_needed_features_norm, TRUE), get_vsn_factors(all_needed_features_norm, TRUE)))
        }

        # Iterate over all molecules
        for (molecule in unique_molecules) {

          # Iterate over all donors
          for (donor in unique_donors) {
            
            # Subset data for the current donor
            donor_data <- full_prep_data[full_prep_data[[Group]] == donor,]
            donor_data <- subset(donor_data, Molecule.Name == molecule)
            donor_data$Area_unmodified <- donor_data$Area
            pure_donor_data_for_norm <- as.data.frame(t(donor_data$Area))
            colnames(pure_donor_data_for_norm) <- donor_data$Sample
            
            if (sum(!is.na(donor_data$Area)) < 3) {
              next
            }

            # Remove the 'X' prefix (if any) from the column names in pure_donor_data_for_norm
            adjusted_colnames <- sub("^X", "", colnames(pure_donor_data_for_norm))
            # Filter columns that match donor_data$Sample
            pure_donor_data_for_norm <- pure_donor_data_for_norm[adjusted_colnames %in% donor_data$Sample]
            
            # Remove the 'X' prefix (if any) from the column names in pqn_vsn_list
            adjusted_colnames <- sub("^X", "", colnames(pqn_vsn_list))
            # Filter columns that match donor_data$Sample
            pqn_vsn_list_donor <- pqn_vsn_list[adjusted_colnames %in% donor_data$Sample]
            
            # Remove the 'X' prefix (if any) from the column names in pqn_vsn_list
            adjusted_colnames <- sub("^X", "", colnames(norm_feature_list))
            # Identify the columns to keep based on the condition
            # Keep the first 4 and last columns and reorder the columns so they match pure_donor_data_for_norm
            keep_columns <- c(1:4, which(adjusted_colnames %in% donor_data$Sample)[match(colnames(pure_donor_data_for_norm), adjusted_colnames[which(adjusted_colnames %in% donor_data$Sample)])], ncol(norm_feature_list))
            # Subset the norm_feature_list to keep the selected columns
            norm_feature_list_donor <- norm_feature_list[, keep_columns]
            
            
            if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'PQN') {
              donor_data$Area <-  unlist(donor_data$Area) / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'pqn_factors',])
              pure_donor_data_for_norm <- pure_donor_data_for_norm / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'pqn_factors',])
              norm_feature_list_donor[-c(1:4,ncol(norm_feature_list_donor))] <- norm_feature_list_donor[-c(1:4,ncol(norm_feature_list_donor))] / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'pqn_factors',])
              
            } else if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'VSN') {
              donor_data$Area <-  unlist(donor_data$Area) / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'vsn_normalization_factors',])
              pure_donor_data_for_norm <- pure_donor_data_for_norm / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'vsn_normalization_factors',])
              norm_feature_list_donor[-c(1:4,ncol(norm_feature_list_donor))] <- norm_feature_list_donor[-c(1:4,ncol(norm_feature_list_donor))] / unlist(pqn_vsn_list_donor[row.names(pqn_vsn_list_donor) == 'vsn_normalization_factors',])
            
            } else if (!is.na(pqn_vsn_none) && pqn_vsn_none == 'fullVSN') {
              list_vsn_data <- vsn_normalize_and_transform(all_needed_features, all_data, df_ids, transitionFile, transitionFile_cluster, normalization_list, add_info_file, donor, molecule)
              donor_data_tmp <- list_vsn_data$donor_data
              donor_data$Area <- donor_data_tmp$Area
              pure_donor_data_for_norm <- list_vsn_data$pure_donor_data_for_norm
              norm_feature_list_donor <- list_vsn_data$norm_feature_list_donor
            }
            
            # Calculate normalization factors for this donor
            norm_factors <- calculate_norm_factors(pure_donor_data_for_norm, method, comb, norm_feature_list_donor)
            
            # Apply normalization
            donor_data$Area_norm <- donor_data$Area / norm_factors
            
            # Make post-PQN/VSN (if selected)
            if (!is.na(pqn_vsn_none) && (pqn_vsn_none == 'post-PQN' | pqn_vsn_none == 'post-VSN')) {
              # Remove the 'X' prefix (if any) from the column names in pqn_vsn_list
              adjusted_colnames <- sub("^X", "", colnames(post_pqn_vsn_list))
              # Filter columns that match donor_data$Sample
              post_pqn_vsn_list_donor <- post_pqn_vsn_list[adjusted_colnames %in% donor_data$Sample]
              
              if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-PQN') {
                donor_data$Area_norm <-  unlist(donor_data$Area_norm) / unlist(post_pqn_vsn_list_donor[row.names(post_pqn_vsn_list_donor) == 'pqn_factors',])
                pure_donor_data_for_norm <- pure_donor_data_for_norm / unlist(post_pqn_vsn_list_donor[row.names(post_pqn_vsn_list_donor) == 'pqn_factors',])
                
              } else if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-VSN') {
                donor_data$Area_norm <-  unlist(donor_data$Area_norm) / unlist(post_pqn_vsn_list_donor[row.names(post_pqn_vsn_list_donor) == 'vsn_normalization_factors',])
                pure_donor_data_for_norm <- pure_donor_data_for_norm / unlist(post_pqn_vsn_list_donor[row.names(post_pqn_vsn_list_donor) == 'vsn_normalization_factors',])
              
              }
            } else if (!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-fullVSN') {
              pure_data_for_norm <- targeted_experiment_data[-c(1:4)] %>% select(-Molecule.Name)
              all_norm_factors <- calculate_norm_factors(pure_data_for_norm, method, comb, norm_feature_list)
              
              # Normalize all samples
              all_needed_features_norm <- all_needed_features
              all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / all_norm_factors
              all_data_norm <- all_data / all_norm_factors  
              
              list_vsn_data <- vsn_normalize_and_transform(all_needed_features_norm, all_data_norm, df_ids, transitionFile, transitionFile_cluster, normalization_list, add_info_file, donor, molecule)
              donor_data_tmp <- list_vsn_data$donor_data
              donor_data$Area_norm <- donor_data_tmp$Area
              pure_donor_data_for_norm <- list_vsn_data$pure_donor_data_for_norm
            }
              
            # Apply scaling to raw values
            donor_data$scaled_conc <- donor_data$Area_unmodified - min(donor_data$Area_unmodified, na.rm = TRUE)
            donor_data$scaled_conc <- donor_data$scaled_conc / 
              max(donor_data$scaled_conc[donor_data$time > 0 & donor_data$time < 180], na.rm = TRUE)
            
            
            # Apply scaling to norm. values
            donor_data$scaled_conc_norm <- donor_data$Area_norm - min(donor_data$Area_norm, na.rm = TRUE)
            donor_data$scaled_conc_norm <- donor_data$scaled_conc_norm / 
              max(donor_data$scaled_conc_norm[donor_data$time > 0 & donor_data$time < 180], na.rm = TRUE)
            
            # Add a late time point and set to 0
            donor_data <- donor_data %>%
              # Create new rows with time = 720 and Area = 0 for each unique Donor and Molecule.Name
              bind_rows(
                donor_data %>%
                  distinct(Donor, Molecule.Name) %>%
                  mutate(Timepoint = 1200+1200, time = 720, scaled_conc = 0, scaled_conc_norm = 0) 
              )
            
            # Fit Bateman function
            fit <- try(nlsLM(scaled_conc_norm ~ bateman_2_compartment_peripheral_no_reflux(time, ka, kel, k12, kel2),
                             data = donor_data,
                             start = fit_params$start,
                             lower = fit_params$lower,
                             upper = fit_params$upper,
                             weights = c(10, rep(1, nrow(donor_data) - 2), 10), # add an additional weight to the first and last point
                             control = list(maxiter = 1000)),
                       silent = TRUE)
            
            fit_raw <- try(nlsLM(scaled_conc ~ bateman_2_compartment_peripheral_no_reflux(time, ka, kel, k12, kel2),
                             data = donor_data,
                             start = fit_params$start,
                             lower = fit_params$lower,
                             upper = fit_params$upper,
                             weights = c(10, rep(1, nrow(donor_data) - 2), 10), # add an additional weight to the first and last point
                             control = list(maxiter = 1000)),
                       silent = TRUE)
            
            # Get fit parameters for all curves
            if (class(fit_raw) != "try-error" & class(fit) != "try-error") {
              params <- as.data.frame(t(coef(fit)))  # Convert to data frame with named columns
              params$Donor <- donor  
              params$molecule <- molecule
              params_raw <- as.data.frame(t(coef(fit_raw)))  # Convert to data frame with named columns
              params_raw$Donor <- donor  
              params_raw$molecule <- molecule  
              
              # Bind to existing dataframe
              if (length("curve_params_raw") == 0) {
                curve_params_raw <- params_raw
                curve_params_norm <- params
              } else {
                curve_params_raw <- rbind(curve_params_raw, params_raw)  # Append new row
                curve_params_norm <- rbind(curve_params_norm, params)  # Append new row
              }
            }
            
            
            # Check for failed fits and penalize accordingly
            if (class(fit_raw) == "try-error" & class(fit) == "try-error") { # both cant be fit
              next
              
            } else if (class(fit_raw) != "try-error" & class(fit) == "try-error") { # only normalized curve cant be fit, penalty
              print("Fitting failed for donor:")
              print(donor)
              print(head(donor_data[c('scaled_conc_norm', 'time')]))  # Debug the input data
              total_residual_error <- max(total_residual_error * 1.1, total_residual_error + 100)
              
              residual_error_raw <- sum(abs(donor_data$scaled_conc - predict(fit_raw, list(time = donor_data$time))) /
                                          (donor_data$scaled_conc  + 1e-4), na.rm = TRUE) # add small epsilon value to prevent division by 0 and inflate error from to small values
              
              if (!is.na(residual_error_raw)) {
                total_residual_error_raw <- total_residual_error_raw + (residual_error_raw * 10) # Multiply to amplify the error
              }
              
              next
              
            } else if (class(fit_raw) == "try-error" & class(fit) != "try-error") { # raw curve cant be fit but normalized can, give reward
                total_residual_error <- min(total_residual_error * 0.9, total_residual_error - 10)
                next
                
            } else { # both can be fit, calculate error
              # Calculate residual error
              residual_error <- sum(abs(donor_data$scaled_conc_norm - predict(fit, list(time = donor_data$time))) /
                                      (donor_data$scaled_conc_norm  + 1e-4), na.rm = TRUE) # add small epsilon value to prevent division by 0 and inflate error from to small values
              residual_error_raw <- sum(abs(donor_data$scaled_conc - predict(fit_raw, list(time = donor_data$time))) /
                                      (donor_data$scaled_conc  + 1e-4), na.rm = TRUE) # add small epsilon value to prevent division by 0 and inflate error from to small values
              
              if (!is.na(residual_error)) {
                total_residual_error <- total_residual_error + (residual_error * 10) # Multiply to amplify the error
              } 
              
              if (!is.na(residual_error_raw)) {
                total_residual_error_raw <- total_residual_error_raw + (residual_error_raw * 10) # Multiply to amplify the error
              }
              
              if (!is.na(residual_error) & !is.na(residual_error_raw) & (residual_error > residual_error_raw * 1.1)) {
                norm_worse_raw_count <- norm_worse_raw_count + 1 # Count how often the norm. curve is worse than the raw one
              }
            }
            
            # Get correction factor between fitted curve and measured point, which can be used to estimate the sweat volumn on the sample
              volumn_corr_factors_raw <- as.data.frame(donor_data$scaled_conc / predict(fit_raw, list(time = donor_data$time)))
              colnames(volumn_corr_factors_raw) <- molecule
              rownames(volumn_corr_factors_raw) <- donor_data$time
              
              volumn_corr_factors_norm <- as.data.frame(donor_data$scaled_conc_norm / predict(fit, list(time = donor_data$time)))
              colnames(volumn_corr_factors_norm) <- molecule
              rownames(volumn_corr_factors_norm) <- donor_data$time
              
              area_values_for_points <- as.data.frame(donor_data$Area_unmodified)
              colnames(area_values_for_points) <- molecule
              rownames(area_values_for_points) <- donor_data$time
              
              
              if (length(volumn_corr_list_raw[[donor]]) == 0) {
                volumn_corr_list_raw[[donor]] <- volumn_corr_factors_raw
                volumn_corr_list_norm[[donor]] <- volumn_corr_factors_norm
                area_values_list[[donor]] <- area_values_for_points
              } else {
                volumn_corr_list_raw[[donor]] <- merge(volumn_corr_list_raw[[donor]], volumn_corr_factors_raw, by='row.names', all=TRUE) %>% column_to_rownames(var="Row.names")
                volumn_corr_list_norm[[donor]] <- merge(volumn_corr_list_norm[[donor]], volumn_corr_factors_norm, by='row.names', all=TRUE) %>% column_to_rownames(var="Row.names")
                area_values_list[[donor]] <- merge(area_values_list[[donor]], area_values_for_points, by='row.names', all=TRUE) %>% column_to_rownames(var="Row.names")
              }
            
              
            # Prepare fit for plotting 
              time_seq <- seq(min(donor_data$time), max(donor_data$time), length.out = 100)
              fit_values <- predict(fit, list(time = time_seq))
              fit_data <- data.frame(time = time_seq, fit_values = fit_values) # Combine time_seq and fit_values into a data frame
  
              fit_values_raw <- predict(fit_raw, list(time = time_seq))
              fit_data_raw <- data.frame(time = time_seq, fit_values_raw = fit_values_raw) # Combine time_seq and fit_values into a data frame
  
              # Combine both raw and normalized data into a single data frame for plotting
              donor_data_combined <- donor_data %>%
                mutate(type = "Normalized", value = scaled_conc_norm) %>%
                bind_rows(
                  donor_data %>%
                    mutate(type = "Raw", value = scaled_conc)  # Assuming 'scaled_conc_raw' is your raw concentration column
                )

            # Plot the combined data
            if (testing) {
              p <- ggplot() +
                # Plot normalized and raw data points
                geom_point(data = donor_data_combined, aes(x = time, y = value, color = type), size = 2) +
                # Overlay the Bateman function fit (normalized data fit)
                geom_line(data = fit_data, aes(x = time, y = fit_values, color = "Fit (Normalized)"), linewidth = 1.5) +
                # Overlay the Bateman function fit (raw data fit)
                geom_line(data = fit_data_raw, aes(x = time, y = fit_values_raw, color = "Fit (Raw)"), linewidth = 1.5) +
  
                # Add titles and labels
                labs(
                  title = paste("Donor:", donor, "\nMolecule:", molecule),
                  x = "Time (min)",
                  y = "Concentration",
                  color = "Data Type"
                ) +
                theme_minimal()


            # Add the plot to the list
              plot_list[[paste(donor, molecule, sep = "_")]] <- p
            
            # For debugging purpose
              print(donor)
              print(fit)
              print(fit_raw)
  
  
              print(residual_error)
              print(total_residual_error)
              print(residual_error_raw)
              print(total_residual_error_raw)
            
              i <- i + 1
            }
          }
        }
        
        
        if (testing) {
          print(length(plot_list))
          print(i)
          
          pdf_file <- file.path(resultsdir, "overview_of_metabolic_profiles_test5.pdf")
          pdf(file = pdf_file, width = 100, height = 40)
          grid.arrange(grobs = plot_list, ncol = 9)
          dev.off()
          
          get_PCA_plot(curve_params_raw[c(1:4)], curve_params_raw, "Donor", title = "PCA for Curve Parameters (Raw Data)", resultsdir, "PCA_curve-parameters-raw.png")
          get_PCA_plot(curve_params_norm[c(1:4)], curve_params_norm, "Donor", title = "PCA for Curve Parameters (Norm Data)", resultsdir, "PCA_curve-parameters-norm-test1.png")
          write.csv(curve_params_raw, paste0(resultsdir, 'curve_params_raw.csv'))
          write.csv(curve_params_norm, paste0(resultsdir, 'curve_params_norm-test1.csv'))
        }
        
        # check how well the features (for wine, potato and the others cluster -> 3 clusters)
          # Apply PQN/VSN if selected
          clustered_data_norm <- clustered_data_raw
          pure_data_for_norm <- targeted_experiment_data[-c(1:4)] %>% select(-Molecule.Name)
          norm_feature_list_norm <- norm_feature_list
          all_needed_features_norm <- all_needed_features
          
          if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'PQN') {
            clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] <-  clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'pqn_factors',])
            pure_data_for_norm <- pure_data_for_norm / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'pqn_factors',])
            norm_feature_list_norm[-c(1:4,ncol(norm_feature_list_norm))] <- norm_feature_list_norm[-c(1:4,ncol(norm_feature_list_norm))] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'pqn_factors',])
            all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'pqn_factors',])

          } else if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'VSN') {
            clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] <-  clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'vsn_normalization_factors',])
            pure_data_for_norm <- pure_data_for_norm / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'vsn_normalization_factors',])
            norm_feature_list_norm[-c(1:4,ncol(norm_feature_list_norm))] <- norm_feature_list_norm[-c(1:4,ncol(norm_feature_list_norm))] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'vsn_normalization_factors',])
            all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / unlist(pqn_vsn_list[row.names(pqn_vsn_list) == 'vsn_normalization_factors',])
            
          } else if (!is.na(pqn_vsn_none) && pqn_vsn_none == 'fullVSN') {
            list_vsn_data <- vsn_normalize_and_transform(all_needed_features, all_data, df_ids, transitionFile, transitionFile_cluster, normalization_list, add_info_file, donor, molecule)
            clustered_data_norm <- list_vsn_data$clustered_data
            pure_data_for_norm <- list_vsn_data$pure_data_for_norm
            norm_feature_list_norm <- list_vsn_data$norm_feature_list
            all_needed_features_norm <- list_vsn_data$all_needed_features_norm
          }
          
          # Calculate normalization factors for all samples
          all_norm_factors <- calculate_norm_factors(pure_data_for_norm, method, comb, norm_feature_list_norm)
          
          # Apply normalization
          clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] <- clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] / all_norm_factors
          all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / all_norm_factors
          
          if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-PQN') {
            clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] <- clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] / unlist(post_pqn_vsn_list[row.names(post_pqn_vsn_list) == 'pqn_factors',])
            all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / unlist(post_pqn_vsn_list[row.names(post_pqn_vsn_list) == 'pqn_factors',])
            
          } else if(!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-VSN') {
            clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] <-  clustered_data_norm[-c(1:4,ncol(clustered_data_norm))] / unlist(post_pqn_vsn_list[row.names(post_pqn_vsn_list) == 'vsn_normalization_factors',])
            all_needed_features_norm[-c(1:9)] <- all_needed_features_norm[-c(1:9)] / unlist(post_pqn_vsn_list[row.names(post_pqn_vsn_list) == 'vsn_normalization_factors',])
            
          } else if (!is.na(pqn_vsn_none) && pqn_vsn_none == 'post-fullVSN') {
            all_data_norm <- all_data / all_norm_factors
            
            list_vsn_data <- vsn_normalize_and_transform(all_needed_features_norm, all_data_norm, df_ids, transitionFile, transitionFile_cluster, normalization_list, add_info_file, donor, molecule)
            clustered_data_norm <- list_vsn_data$clustered_data
            pure_data_for_norm <- list_vsn_data$pure_data_for_norm
            norm_feature_list_norm <- list_vsn_data$norm_feature_list
            all_needed_features_norm <- list_vsn_data$all_needed_features_norm
          }
          
          # Get raw cluster data
          avg_silhouette_per_cluster_raw <- cluster_data_raw[[1]]
          sd_silhouette_per_cluster_raw <- cluster_data_raw[[2]]
          
          if (testing) {
            pdf(paste0(resultsdir,"cluster_hierarchical_raw_data_automatic_plots.pdf"), width = 25, height = 12)
            full_prep_cluster_data_raw <- prepare_data_for_plot(clustered_data_raw, add_info_file)
            plots <- plot_clusters(cluster_data_raw, full_prep_cluster_data_raw) # Generate plots for the current clustering
            grid.arrange(grobs = plots, ncol = 5, nrow = 3)
            dev.off()
          }
          
          # Cluster normalized data
          fold_change_df_norm <- get_fold_change_dataframe(clustered_data_norm, add_info, Group_list = list('Intervention'))
          cluster_data_norm <- get_silhouette_scores_with_hierarchical_clustering(fold_change_df_norm)
          avg_silhouette_per_cluster_norm <- cluster_data_norm[[1]]
          sd_silhouette_per_cluster_norm <- cluster_data_norm[[2]]
          
          if (testing) {
            pdf(paste0(resultsdir,"cluster_hierarchical_norm_test5_data_automatic_plots.pdf"), width = 25, height = 12)
            full_prep_cluster_data_norm <- prepare_data_for_plot(clustered_data_norm, add_info_file)
            plots <- plot_clusters(cluster_data_norm, full_prep_cluster_data_norm) # Generate plots for the current clustering
            grid.arrange(grobs = plots, ncol = 5, nrow = 3)
            dev.off()
          }
          
          # Problem is that the cluster numbers are assigned randomly so a direct comparison is hard, 
          # tried with distance matrix but reordering can sometimes assign same cluster twice
          # so now just make a total comparison
          # bigger avg. silhouette is good (1=perfect clustering, 0 overlapping, -1 poor cluster)
          # if norm. avg_silhouette bigger than raw it is good, so want negative value
          avg_diff <- mean(avg_silhouette_per_cluster_raw[,2], na.rm = TRUE) - mean(avg_silhouette_per_cluster_norm[,2], na.rm = TRUE) 
          # smaller sd. silhouette is good, the smaller the better the fit into cluster
          # if norm. avg_silhouette smaller than raw it is good, so want negative value
          sd_diff <- mean(sd_silhouette_per_cluster_norm[,2], na.rm = TRUE) - mean(sd_silhouette_per_cluster_raw[,2], na.rm = TRUE)
          
          total_residual_error_curve_only <- total_residual_error
          
          # if avg/sd_diff are negative then clustering got better so reduce error, otherwise increase it, depending on the change
          total_residual_error <- total_residual_error * (1 + (avg_diff * 0.5)) * (1 + (sd_diff * 2))
          
          # If all curves are better, give a slight extra reward
          if (norm_worse_raw_count == 0) {
            total_residual_error <- total_residual_error * 0.9
          }
          
          # If clustering gets worse then add number of worse curves as slight punishment on top
          if (avg_diff > 0 & sd_diff > 0) { # Clusters are not well defined and not tight
            if (norm_worse_raw_count > 5) { # only add additional score for worse curves if cluster is also bad
              total_residual_error <- total_residual_error * (1 + 0.004 * norm_worse_raw_count)
            }
          }
          
        # Check how well the correction factors for the same sample over different molecules fit together
          mean_sd_raw <- compute_sd_after_filtering(volumn_corr_list_raw, area_values_list, threshold = 3000)
          mean_sd_norm <- compute_sd_after_filtering(volumn_corr_list_norm, area_values_list, threshold = 3000)
          
          # If we get better than sd_norm smaller than sd_raw, so we are negative
          corr_factor_diff <- mean_sd_norm - mean_sd_raw
          
          total_residual_error <- total_residual_error * (1 + (corr_factor_diff * 0.5))
        
          
        
        # Make LIMMA test and see if after normalization more features are significant
          full_data_norm <- prepare_data_for_limma(all_needed_features_norm, add_info)
          limma_norm <- get_limma_fit_for_signifint_changes_over_time(full_data_norm)
          significant_features_norm <- length(limma_norm$adj.P.Val[limma_norm$adj.P.Val < 0.05])
          if (significant_features_norm > 0) {
            significant_features_diff <- significant_features_raw / significant_features_norm # want it to get smaller, so norm has more sign. features 
          } else {
            significant_features_diff <- significant_features_raw / 1 # if no significant features are found than we set to 0, otherwise diff by 0
          }
          
          total_residual_error <- total_residual_error * significant_features_diff
          
          
        if (bayesoptim) {
          # Return the negative residual error (Bayesian optimization maximizes it values so we want to get them closest to 0 as possible)
          return(-total_residual_error)
        } else {
          return(list(Score = -total_residual_error, pqn_vsn_used = toString(pqn_vsn_none), method_name = toString(method), comb = toString(norm_feature_list$Molecule.Name[weights > 0.5]), 
                      norm_curve_error = total_residual_error_curve_only, raw_curve_error = total_residual_error_raw, avg_silhouette_diff = avg_diff, sd_silhouette_diff = sd_diff, 
                      corr_factor_diff = corr_factor_diff, significant_features_diff = significant_features_diff, significant_features_raw = significant_features_raw, significant_features_norm = significant_features_norm))
        }
        
      }
      
      
      # Define normalization methods
      #norm_methods <- c("non-weighted", "Rank2", "Median", "Sqrt", "Log2")
      norm_methods <- c("Median")
      
      #pqn_vsn_methods <- c(NA, "PQN", "post-PQN", "VSN", "post-VSN", "fullVSN", "post-fullVSN")
      #pqn_vsn_methods <- c(NA, "PQN", "post-PQN", "VSN", "post-VSN")
      pqn_vsn_methods <- c("post-VSN")
      
      # Define save interval and output file
      last_save_time <- Sys.time()
      next_save_interval <- runif(1, min = 580, max = 620)  # Random time in seconds
      interim_results_file <- paste0(resultsdir,"interim_results.csv")
      
      # Register a parallel backend
      cl <- initialize_cluster(61, varlist = c("norm_feature_list", "targeted_experiment_data", "column_to_rownames", "add_info", "add_info_file", "resultsdir",
                                               "compute_sd_after_filtering", "get_fold_change_dataframe", "get_silhouette_scores_with_hierarchical_clustering",
                                               "replace_na_and_Inf", "silhouette", "bateman", "calculate_norm_factors", "bateman_2_compartment_peripheral_no_reflux",
                                               "get_pqn_factors", "get_vsn_factors", "vsn2", "predict", "vsn_normalize_and_transform", "normalize_names", "str_detect",
                                               "df_ids", "transitionFile", "transitionFile_cluster", "normalization_list", "extract_feature_list", "formula2mz", "topTable",
                                               "prepare_data_for_plot", "convert_to_minutes", "all_data", "get_limma_fit_for_signifint_changes_over_time", "prepare_data_for_limma")) 
      
      clusterExport(cl, list("full_prep_data", "unique_donors", "norm_methods", "Group", "clustered_data_raw", "unique_molecules", "bayesoptim", 
                             "fit_params", "evaluate_combination", "pqn_vsn_list", "cluster_data_raw", "significant_features_raw",
                             "next_save_interval", "last_save_time", "interim_results_file"), envir = environment())
      
      registerDoParallel(cl)
      
      print("Cluster prepared. Calculation is starting.")
      start_time <- Sys.time()  # Record start time
      
      if (bayesoptim) {
        # Set possible search space for Bayesian Optimization
        search_space <- list(
          # For method, the range is between 1 and the length of norm_methods (integer encoding)
          method = c(1,length(norm_methods)+0.99),  # Map methods to integers
          
          # To decide if pqn or vsn or non at all get added too
          pqn_vsn = c(0,length(norm_methods)+0.99),
          
          # For each molecule, the weight can range from 0 to 1 (both inclusive)
          m1 = c(0,1),
          m2 = c(0,1),
          m3 = c(0,1),
          m4 = c(0,1),
          m5 = c(0,1),
          m6 = c(0,1),
          m7 = c(0,1),
          m8 = c(0,1),
          m9 = c(0,1),
          m10 = c(0,1),
          m11 = c(0,1),
          m12 = c(0,1),
          m13 = c(0,1),
          m14 = c(0,1)
        )
        
        opt_res <- bayesOpt(
          FUN = function(method, pqn_vsn, m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14) {
            weights <- c(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14)
            
            # Decode method from integer to string
            method_decoded <- norm_methods[floor(method)]
            
            pqn_vsn_none <- pqn_vsn_methods[floor(pqn_vsn)]
            
            # Evaluate the combination
            score <- evaluate_combination(method_decoded, weights, pqn_vsn_none)
            return(list(Score = score, pqn_vsn_used = toString(pqn_vsn_none), method_name = toString(method_decoded), comb = toString(norm_feature_list$Molecule.Name[weights > 0.5])))
          },
          bounds = search_space,
          initPoints = 60,    # Initial random evaluations
          iters.n = 1680,     # Number of iterations (14 epochs = 1680)
          iters.k = 120,      # Number of times to sample FUN at each optimization step
          acq = "ei",         # Acquisition function
          eps = -0.1,         # Parameter to focus (a little) more on exploration to avoid getting stuck in a local minimum 
          convThresh = 5e+07, # Stricter convergence, can help prevent oversampling regions that aren't improving
          acqThresh = 0.6,    # Allow local optima with a reasonable utility to be considered, avoiding over-reliance on a single region
          otherHalting = list(timeLimit = 86400, minUtility = 1e-2), # Criteria to stop after 1 day (86400 s) or if the improvement is negligible 
          parallel = TRUE,    # Parallel execution
          plotProgress = TRUE
        )
      } else {
        
        # Create a empty variables 
        binary_combinations <- list()
  
        # If a interim results file from a previous calculation exists than we dont need to calculate those combinations again
        if (file.exists(paste0(resultsdir, "interim_results.csv"))) {
          interim_results <- read.csv(paste0(resultsdir, "interim_results.csv"))
          already_calculated <- unique(as.character(interim_results$comb))  # Extract existing combinations
        } else {
          interim_results <- NULL
          already_calculated <- character(0)  # Empty character vector if nothing calculated yet
        }
        
        # Generate all combinations and convert to binary vectors
        for (i in 8:8) {
          combs <- combn(norm_feature_list$Molecule.Name, i, simplify = FALSE)
          
          # Convert each combination to a binary vector
          for (comb in combs) {
            # Create the "name" string as stored in interim_results$comb
            comb_string <- toString(comb)  # comma-separated names
            
            # Only add if not already calculated
            if (!(comb_string %in% already_calculated)) {
              binary_vector <- as.integer(norm_feature_list$Molecule.Name %in% comb)
              binary_combinations <- append(binary_combinations, list(binary_vector))
            }
          }
        }
        
        # Run the combinations in parallel
        opt_res <- foreach(weights = binary_combinations, .combine = rbind, .packages = c("data.table")) %dopar% { # dopar declares to run it parallel
          result_list <- list()
          
          for (method_decoded in norm_methods) {
            for (pqn_vsn_none in pqn_vsn_methods) {
              # Try-catch block to catch errors inside the parallel loop
              res <- tryCatch({
                evaluate_combination(method_decoded, unlist(weights), pqn_vsn_none)
              }, error = function(e) {
                # Capture error details
                message("Error in evaluate_combination: ", conditionMessage(e))
                message("Problematic inputs: method_decoded=", method_decoded, 
                        ", pqn_vsn_none=", pqn_vsn_none, ", weights=", paste(unlist(weights), collapse = ", "))
                
                # Return a structured row with NA for scores but preserving method details
                return(as.data.table(list(Score = NA,
                                         pqn_vsn_used = toString(pqn_vsn_none),
                                         method_name = toString(method_decoded),
                                         comb = toString(norm_feature_list$Molecule.Name[weights > 0.5]),
                                         norm_curve_error = NA,
                                         raw_curve_error = NA,
                                         avg_silhouette_diff = NA,
                                         sd_silhouette_diff = NA,
                                         corr_factor_diff = NA,
                                         significant_features_diff = NA,
                                         significant_features_raw = NA,
                                         significant_features_norm = NA
                                       )))
              })
              
              # Append results
              result_list[[length(result_list) + 1]] <- as.data.frame(res)
            }
          }
          
          # Convert results to data.table
          result_dt <- rbindlist(result_list, fill = TRUE)
          
          # Check if it's time to save results
          if (as.numeric(difftime(Sys.time(), last_save_time, units = "secs")) > next_save_interval) {
            fwrite(result_dt, interim_results_file, append = TRUE)
            
            last_save_time <- Sys.time() # Update the last save time
            next_save_interval <- runif(1, min = 580, max = 620) # Set a new random save interval
            
            message("Saved interim results at ", format(last_save_time), 
                    " - Next save in ", round(next_save_interval, 2), " seconds")
          }
          
          return(result_dt)
        }
        
        #opt_res <- rbindlist(opt_res, fill = TRUE)
      }
      
      print(paste('Calculation took:', Sys.time() - start_time))
      
      # Stop the parallel cluster after optimization
      stopCluster(cl)
      
      registerDoSEQ() # Deregister the backend
      
      # if there are already calculate results than bind it to the new ones
      if (!is.null(interim_results)) {
        opt_res <- rbind(interim_results, opt_res)
      }
      
      # Save the results
      if(bayesoptim) {
        write.csv(opt_res$scoreSummary, paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable', run, '.csv'))
      } else {
        fwrite(opt_res, paste0(resultsdir, 'bayesOpt_normalising_scoreSummaryTable_all_combinations.csv'))
      }
      
      # Return the best parameters
      return(opt_res$scoreSummary)
    }
  
  
  
  # Function to plot all significantly changed molecules
  plot_significant_changed_features <- function(all_needed_features, foldername, resultsdir, info_file_dir) {
    # Define the folder path
    folder_path <- file.path(resultsdir, foldername)
    
    # Get all CSV files from the folder
    csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    
    # Initialize an empty vector to store significant feature names
    significant_features <- data.frame(Annotation = character(), adj.P.Val = numeric())
    
    # Loop through each file and extract significant features
    for (file in csv_files) {
      data <- fread(file)
      
      # Filter based on logFC and adjusted p-value
      sig_data <- data %>%
        filter((logFC < -1 | logFC > 1) & adj.P.Val < 0.05) %>%
        select(Annotation, adj.P.Val)  # Keep only relevant columns
      
      # Append to significant_features dataframe
      significant_features <- rbind(significant_features, sig_data)
    }
    # Sort by adj.P.Val (smallest first)
    significant_features <- significant_features %>% arrange(adj.P.Val)
    
    # Keep only the first occurrence of each unique feature (preserving smallest adj.P.Val)
    significant_features <- significant_features %>% distinct(Annotation, .keep_all = TRUE)
    
    write.csv(significant_features, paste0(resultsdir, "all-significantly_regulated_features.csv"))
    
    # Get unique feature names
    unique_names <- unique(significant_features$Annotation)
    
    # Open a PDF to save plots
    pdf(file = file.path(resultsdir, "significantly_regulated_feature_plots.pdf"))
    
    # Loop through each unique feature and generate plots
    for (feature in unique_names) {
      # Filter data for the specific feature
      significant_feature <- all_needed_features %>% filter(Annotation == feature)
      
      # Prepare data for plotting
      plot_data <- prepare_data_for_plot(significant_feature, info_file_dir)
      
      # Generate plot
      p <- ggplot(plot_data, aes(x = Timepoint, y = log2(Area), 
                                 color = Intervention, linetype = Sex, 
                                 group = interaction(Donor, Intervention))) +
        geom_line() +
        theme_minimal() +
        labs(title = paste(feature, "over Time by Donor and Intervention"),
             x = "Timepoint",
             y = "log2(Area)")
      
      # Print plot to PDF
      print(p)
    }
    
    # Close the PDF
    dev.off()
  }
  
  
  
##############################################################
## Functions to Check for possible Normalization Parameters ##
##############################################################
   
  
# Function: Calculate and evaluate CV to check for possible Housekeeping Molecules
# Input: experiment_data, stable_threshold, min_sample_count
# Output: Dataframe of consistently stable molecules
  evaluate_stability_of_molecule_paper_list <- function(experiment_data, stable_threshold = 0.25, min_sample_count = 10) {
    # Computing mean, standard deviation, and count of areas for each molecule and paper
    molecule_stats_by_paper <- experiment_data %>%
      group_by(Molecule.Name, paper) %>%
      summarise(mean_Area = mean(Area, na.rm = TRUE),
                sd_Area = sd(Area, na.rm = TRUE),
                count = n(),
                .groups = 'drop')  # Clean-up the grouping structure
    
    # Calculate the coefficient of variation (CV) for each molecule within each paper
    molecule_stats_by_paper$CV <- molecule_stats_by_paper$sd_Area / molecule_stats_by_paper$mean_Area
    
    # Define a stability threshold for the CV
    stable_molecules_by_paper <- molecule_stats_by_paper %>%
      filter(CV < stable_threshold, count > min_sample_count)  # Filter out stable molecules with enough data points
    
    # Results can be examined here
    print(stable_molecules_by_paper)
    
    # If you want to consider a molecule stable only if it is consistent across both papers, you might do:
    stable_molecules_across_papers <- stable_molecules_by_paper %>%
      group_by(Molecule.Name) %>%
      filter(n_distinct(paper) == 2) %>%
      summarise(avg_CV = mean(CV), .groups = 'drop')  # Compute average CV across papers
    
    # Filter molecules that are stable in both papers under the specified threshold
    consistently_stable_molecules <- stable_molecules_across_papers %>%
      filter(avg_CV < stable_threshold)
    
    return(consistently_stable_molecules)
  }
  
  
  

  
  # Function: Check for correlations in the untargeted dataset and save to CSV file
  # Input: all_data (whole untargeted dataset), df_ids, resultsdir, threshold, at_least_n_samples, p_value_threshold
  # Output: CSV file of correlating features
  identify_correlating_features <- function(all_data, df_ids, datadir, resultsdir, threshold = 0.90, at_least_n_samples = 100, 
                                            p_value_threshold = 0.05, filename = "correlating_features.csv") {
    
    # get features and merge with ID
    SIRIUS <- process_SIRIUS_table(datadir)
    df_features <- merge(df_ids, SIRIUS, by = "id", all.x = TRUE)
    #df_features <- df_features[!duplicated(df_features$id), ] # remove duplicate annotated features
    
    # Get whole untargeted dataset as dataframe
    experiment_data_untargeted <- cbind(df_ids[rownames(all_data), ], all_data)
    experiment_data_untargeted <- na.omit(experiment_data_untargeted) # Remove NA values
    transformed_data <- t(experiment_data_untargeted[, -c(1:4)]) # Exclude columns 'id', 'rt', 'mz', 'charge' and transform so features are in column
    corr_matrix <- rcorr(as.matrix(transformed_data), type = "spearman")
    diag(corr_matrix$P)[is.na(diag(corr_matrix$P))] <- 1e-10 # Replace NA values in the diagonal with 1
    
    # Function to identify features with correlations meeting criteria
    identify_correlated_features <- function(corr_matrix, p_values, threshold, at_least_n_samples, p_value_threshold) {
      correlated_features <- logical(nrow(corr_matrix))
      
      # Iterate over rows of the correlation matrix and check if conditions are met
      for (i in seq_len(nrow(corr_matrix))) {
        above_threshold <- abs(corr_matrix[i, ]) >= threshold
        below_p_value_threshold <- p_values[i, ] <= p_value_threshold
        correlated_features[i] <- sum(above_threshold & below_p_value_threshold) >= at_least_n_samples
      }
      
      # Extract names of correlated features using which()
      correlated_indices <- which(correlated_features)
      correlated_features <- rownames(corr_matrix)[correlated_indices]
      
      return(correlated_features)
    }
    
    # Identify features with correlations meeting criteria
    correlated_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), threshold, at_least_n_samples, p_value_threshold)
    correlated_feature_ids <- df_ids[df_ids$id %in% correlated_features, "id"]
    #correlated_feature_names <- df_features[df_features$id %in% correlated_feature_ids, "Annotation"]
    #print(correlated_feature_names)
    
    # Look at the correlating features in more detail
    data_length <- ncol(experiment_data_untargeted)
    correlating_data <- experiment_data_untargeted[experiment_data_untargeted$id %in% correlated_feature_ids, ]
    correlating_data$mean_area <- rowMeans(correlating_data[, 5:data_length])
    correlating_data <- correlating_data[, -c(5:data_length)]
    correlating_data <- correlating_data[order(-correlating_data$mean_area), ]
    correlating_data$Annotation <- df_features[df_features$id %in% correlated_feature_ids, "Annotation"]
    #print(correlating_data)
    
    # Write correlating data to CSV file
    write.csv(correlating_data, file = paste(resultsdir, filename), row.names = FALSE)
  }
  
  
  
  # Function: Check for correlations in the untargeted dataset and save to CSV file
  # Input: datadir and resultsdir, where to search for datafiles and where to save them
  # Output: CSV file of correlating features
  identify_correlating_features_fixed_thresholds <- function(datadir, resultsdir) {
    
    # get all significant feature data
    significant_features_object <- get_significant_abundant_features(datadir)
    all_significant_feature_data <- significant_features_object$all_significant_feature_data
    
    # need the charge (could implement in previous function, yet is used often, so just do it ones here separatly)
    exp_data <- process_merged_files(datadir)
    df_ids <- exp_data$df_ids
    significant_feature_data <- merge(all_significant_feature_data, df_ids[, c('id', 'charge')], by='id', all.x=TRUE)
    
    # filter features
      # filter out to small features (row median < 5e3)
      significant_feature_data <- significant_feature_data %>%
                                  filter(apply(.[, 10:(ncol(.)-1)], 1, median) > 5e3)      
      # filter to get features that are in every sample
      significant_feature_data <- significant_feature_data %>% filter(if_all(10:(ncol(significant_feature_data)-1), ~. > 1e3))      
      # filter out features in injection peak (rt < 0.6)
      significant_feature_data <- significant_feature_data[significant_feature_data$rt >= 0.6, ]
      # filter duplicates, use the one with highest area
      significant_feature_data <- significant_feature_data %>%
                                  mutate(sum_column = rowSums(.[-(1:9)])) %>%  # Create sum column from column 10 to ncol - 1
                                  group_by(Annotation) %>%
                                  filter(sum_column == max(sum_column)) %>%
                                  ungroup() %>%
                                  select(-sum_column)
      # filter feature groups that occur several times (probably same feature with several fragmentations)
      # filter +1, -1 charge seperatly
        # Group by feature_group  
        grouped_df <- significant_feature_data %>%
                      group_by(feature_group)
        
        # Filter out rows occurring only once, keep only duplicate ones
        filtered_df <- grouped_df %>%
                       filter(n() > 1) %>%
                       ungroup()
        
        # Filter out rows with charge equal to 1 and -1 separately
        charge_1_df <- filtered_df %>%
                       filter(charge == 1) %>%
                       mutate(row_mean = rowMeans(select(., 10:(ncol(.) - 1)), na.rm = TRUE)) %>%
                       group_by(feature_group) %>%
                       arrange(feature_group, desc(row_mean)) %>%
                       filter(row_number() == 1) %>%
                       ungroup() %>%
                       select(-row_mean)
        
        charge_minus_1_df <- filtered_df %>%
                             filter(charge == -1) %>%
                             mutate(row_mean = rowMeans(select(., 10:(ncol(.) - 1)), na.rm = TRUE)) %>%
                             group_by(feature_group) %>%
                             arrange(feature_group, desc(row_mean)) %>%
                             filter(row_number() == 1) %>%
                             ungroup() %>%
                             select(-row_mean)
        
        # Combine the results
        final_df <- bind_rows(grouped_df %>% filter(n() == 1),  # Rows occurring only once
                              charge_1_df,  # Rows with charge 1
                              charge_minus_1_df  # Rows with charge -1
                              )
      
      
      #unique_feature_groups <- unique(final_df$feature_group)
    
    # transform and prepare data for correlation matrix
      final_df <- as.data.frame(final_df)
      rownames(final_df) <- final_df$id
      final_df <- final_df[, -c(1:9)]
      final_df <- final_df[, -ncol(final_df)]
      transformed_data <- t(final_df) # Exclude all non-area columns and transform so features are in column
    
    
    corr_matrix <- rcorr(as.matrix(transformed_data), type = "spearman")
    diag(corr_matrix$P)[is.na(diag(corr_matrix$P))] <- 1e-10 # Replace NA values in the diagonal with 1
    
    
    # Function to identify features with correlations meeting criteria
      identify_correlated_features <- function(corr_value_matrix, p_values, threshold, at_least_n_samples, p_value_threshold) {
        correlated_features <- logical(nrow(corr_value_matrix))
        
        # Iterate over rows of the correlation matrix and check if conditions are met
        for (i in seq_len(nrow(corr_value_matrix))) {
          above_threshold <- abs(corr_value_matrix[i, ]) >= threshold
          below_p_value_threshold <- p_values[i, ] <= p_value_threshold
          correlated_features[i] <- sum(above_threshold & below_p_value_threshold) >= at_least_n_samples
        }
        
        # Extract names of correlated features using which()
        correlated_indices <- which(correlated_features)
        correlated_features <- rownames(corr_value_matrix)[correlated_indices]
        
        return(correlated_features)
      }
    
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                        threshold = 0.8, at_least_n_samples = 5, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.8_p0.05_n5.csv"), row.names = FALSE)
    
    
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                         threshold = 0.9, at_least_n_samples = 5, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.9_p0.05_n5.csv"), row.names = FALSE)
    
    
    # do the same with other parameters
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                         threshold = 0.9, at_least_n_samples = 10, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.9_p0.05_n10.csv"), row.names = FALSE)
    
    
    # do the same with other parameters
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                        threshold = 0.95, at_least_n_samples = 10, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.95_p0.05_n10.csv"), row.names = FALSE)
    
    
    # do the same with other parameters
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                        threshold = 0.98, at_least_n_samples = 10, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.98_p0.05_n10.csv"), row.names = FALSE)
    
    
    # do the same with other parameters
    # Identify features with correlations meeting criteria
    correlating_features <- identify_correlated_features(round(corr_matrix$r, 4), round(corr_matrix$P, 4), 
                                                         threshold = 0.5, at_least_n_samples = 100, p_value_threshold = 0.05)
    correlating_feature_details <- all_significant_feature_data[all_significant_feature_data$id %in% correlating_features, ]
    
    # Write correlating data to CSV file
    write.csv(correlating_feature_details, file = paste(resultsdir, "correlating_features_t0.5_p0.05_n100.csv"), row.names = FALSE)
    
  } 
  
  
  
  
  # Function to show the distribution of normalization features as weighted network
  print_norm_molecule_network_to_pdf <- function(results, cutoff = NA, remove_duplicates = FALSE, abbrevate = FALSE, abbrevate_len = 10, run = 1, molecule_list = NULL) {
    print(paste('Showing distribution of cutoff:', cutoff, 'and removing duplicates:', remove_duplicates))
    
    # Select specific columns
    results_filtered <- results[, .(Score, pqn_vsn_used, method_name, comb)]
    
    # Filter rows based on Score condition
    min_score <- min(abs(results_filtered$Score))
    if (!is.na(cutoff) ) {
      results_filtered <- results_filtered[abs(results_filtered$Score) < min_score * cutoff]
    }
    
    # Remove duplicte points 
    if (remove_duplicates) {
          results_filtered <- results_filtered[!duplicated(results_filtered, by = c("comb", "method_name", "pqn_vsn_used")), ]
    }
    
    # Count occurrences in the 'method_name' column
    method_counts <- results_filtered[, .N, by = method_name]
    print(method_counts)
    
    # Count occurrences in the 'pqn_vsn' column
    pqn_vsn_counts <- results_filtered[, .N, by = pqn_vsn_used]
    print(pqn_vsn_counts)
    
    results_filtered$comb <- iconv(results_filtered$comb, from = "latin1", to = "UTF-8", sub = "") # to remove some imcompatible signs
    
    # Split 'comb' into individual molecules
    results_filtered[, molecule_list := strsplit(comb, ", ")]
    # Get a unique list of all molecules
    all_molecules <- unique(unlist(results_filtered$molecule_list))
    # Initialize a co-occurrence matrix with 0s
    co_occurrence_matrix <- matrix(0, nrow = length(all_molecules), ncol = length(all_molecules))
    rownames(co_occurrence_matrix) <- all_molecules
    colnames(co_occurrence_matrix) <- all_molecules
    # Initialize the molecule count vector with zeros
    molecule_count <- rep(0, length(all_molecules))
    names(molecule_count) <- all_molecules
    # Fill the matrix with co-occurrence counts
    for (mol_list in results_filtered$molecule_list) {
      # Count how often each molecule appears in this row
      for (mol in mol_list) {
        molecule_count[mol] <- molecule_count[mol] + 1
      }
      
      # Create pairs and update the co-occurrence matrix
      if (length(mol_list) >= 2) {
        pairs <- combn(mol_list, 2, simplify = FALSE)
        for (pair in pairs) {
          co_occurrence_matrix[pair[1], pair[2]] <- co_occurrence_matrix[pair[1], pair[2]] + 1
          co_occurrence_matrix[pair[2], pair[1]] <- co_occurrence_matrix[pair[2], pair[1]] + 1
        }
      }
    }
    # Convert to a percentage of the maximum possible occurrences
    total_percentage <- (molecule_count / nrow(results_filtered)) * 100
    # Add molecule information to a data.frame for the graph
    molecule_info <- data.frame(
      molecule = all_molecules,
      percentage = total_percentage
    )
    # Adjust node sizes based on percentage (e.g., scale up for visualization)
    molecule_info$node_size <- molecule_info$percentage * 0.5
    
    if(!is.null(molecule_list)) {
      molecule_info <- merge(molecule_info, molecule_list, 
                             by.x = "molecule", by.y = "Molecule.Name", 
                             all.x = TRUE)
      
      fixed_colors <- c(
        "selected feature" = "skyblue2",
        "linear sweat feature" = "green4"
      )
      
      # Identify any other origins
      other_origins <- setdiff(unique(molecule_info$origin), names(fixed_colors))
      
      # Create a palette for the other origins
      if(length(other_origins) > 0) {
        palette_colors <- brewer.pal(max(3, length(other_origins)), "Set2")[1:length(other_origins)]
        names(palette_colors) <- other_origins
        # Combine fixed and palette-assigned colors
        origin_colors <- c(fixed_colors, palette_colors)
      
      } else {
        origin_colors <- fixed_colors
      }
      
      
    }
    
    
    # Create a network plot
    # Convert the co-occurrence matrix to an edge list
    edge_list <- which(co_occurrence_matrix > 0, arr.ind = TRUE)
    pair_counts <- data.frame(
      molecule1 = rownames(co_occurrence_matrix)[edge_list[, 1]],
      molecule2 = colnames(co_occurrence_matrix)[edge_list[, 2]],
      count = co_occurrence_matrix[edge_list]
    )
    
    # Remove duplicate edges (as the matrix is symmetric)
    pair_counts <- pair_counts[pair_counts$molecule1 < pair_counts$molecule2, ]
    
    # Create the network graph
    g <- graph_from_data_frame(d = pair_counts, directed = FALSE, vertices = molecule_info)
    
    # Set node sizes and edge widths
    V(g)$size <- molecule_info$node_size
    if (abbrevate == TRUE) {
      abbreviated_names <- abbreviate(V(g)$name, minlength = abbrevate_len)
      V(g)$label <- paste(abbreviated_names, "\n",
                          round(molecule_info$percentage[match(V(g)$name, molecule_info$molecule)], 1), "%")
    } else {
      V(g)$label <- paste(V(g)$name, "\n",
                          round(molecule_info$percentage[match(V(g)$name, molecule_info$molecule)], 1), "%")
    }
    E(g)$width <- pair_counts$count / max(pair_counts$count) * 5  # Scale edge width
    E(g)$label <- pair_counts$count                              # Add edge labels
    
    # Set the output file
    if (is.na(cutoff) ) {
      if(remove_duplicates) {
        pdf(paste0(resultsdir,"molecule-combinations_network_plot_all_without-duplicates_run", run,".pdf"))
      } else {
        pdf(paste0(resultsdir,"molecule-combinations_network_plot_all_with-duplicates_run", run,".pdf"))
      }
    } else {
      if(remove_duplicates) {
        pdf(paste0(resultsdir,"molecule-combinations_network_plot_min-", cutoff, "_without-duplicates_run", run,".pdf"))
      } else {
        pdf(paste0(resultsdir,"molecule-combinations_network_plot_min-", cutoff, "_with-duplicates_run", run,".pdf"))
      }
    }
    
    
    if(is.null(molecule_list)) {
      # Plot the graph
      plot(
        g,
        vertex.label = V(g)$label,          # Use the new labels with percentages
        vertex.label.cex = 0.75,
        edge.label = E(g)$label,           # Co-occurrence counts as edge labels
        edge.label.cex = 0.7,
        edge.color = "grey",
        vertex.color = "skyblue",
        main = "Molecule Co-Occurrence Network"
      )
    } else {
      # Assign node colors using merged `origin` data
      V(g)$color <- origin_colors[molecule_info$origin[match(V(g)$name, molecule_info$molecule)]]
      
      
      # Before plotting the network, set up the plotting area to add a legend
      layout(matrix(1:2, ncol = 2), widths = c(4, 2))  # Make space for legend
      
      # Plot the graph
      plot(
        g,
        vertex.label = V(g)$label,
        vertex.label.cex = 0.75,
        edge.label = E(g)$label,
        edge.label.cex = 0.7,
        edge.color = "grey",
        vertex.color = V(g)$color,
        main = "Molecule Co-Occurrence Network"
      )
      
      # Add legend
      par(mar = c(0, 0, 0, 0))
      plot.new()
      legend("center", legend = names(origin_colors), fill = origin_colors, cex = 0.8, bty = "n", title = "Origin")
    }
    
    
    
    # Save the plot and close the device
    dev.off()
  }
  
  
  