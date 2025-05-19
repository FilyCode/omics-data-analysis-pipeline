#################################################################
# THIS IS A FILE WITH VARIOUS FUNCTIONS USED IN A METABOLOMIC FEATURE ANALYSIS 
# THIS FILE SHOULD CONTAIN ANYTHING NEEDED TO MAKE A FULL DATA ANALYSIS
# FUNCTIONS RANGE FROM DATA PREPARATION (READ IN, WORK WITH IT, SAVE)
# REMOVING NOISE, VARIOUS NORMALIZATIONS, OUTLIER DETECTION, STATISTICS
# DIFFERENT PLOTTING FUNCTIONS TO VISUALIZE THE DATA
#
# File currently needs to be used together with 'MS-data-analysis_functions.R'
#################################################################


##############
# Section to import all needed libraries
##############

  # Data Manipulation and Wrangling
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(data.table)
  library(readr)
  
  # Statistical Analysis and Modeling
  library(limma)
  library(sva)
  library(qmtools)
  library(pls)
  library(caret)
  library(matrixStats)
  library(cluster)
  library(imputeLCMD)
  
  # Visualization
  library(ggplot2)
  library(ggrepel)
  library(ggfortify)
  library(ggtext)
  library(factoextra)
  library(EnhancedVolcano)
  library(gridExtra)
  library(viridis)
  library(sjPlot)
  library(igraph)
  
  # Bioinformatics
  library(Biobase)
  library(SummarizedExperiment)
  library(MSnbase)
  library(MetaboCoreUtils)
  library(MsCoreUtils)
  library(cliqueMS)
  library(webchem)
  library(vsn)
  
  # (Nonlinear) Regression Analysis
  library(minpack.lm)
  library(nlstools)
  library(nlsMicrobio)
  library(glmnet)
  library(ParBayesianOptimization)

  # Spreadsheet Management
  library(openxlsx)
  
  # Miscellaneous
  library(igraph)
  library(Hmisc)
  library(pcaMethods)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(respR)
  library(httr)
  library(jsonlite)
  library(RColorBrewer)


##############
# Section to read in and prepare data, filter, normalize features, check for outliers
##############

# Function: Process SIRIUS annotations table to extract important information and prepare for downstream analysis
# Input: datadir (directory containing SIRIUS compound_identifications.tsv file)
# Output: Data frame containing processed SIRIUS annotations, with columns representing essential information such as feature ID, annotation, retention time, ion mass, molecular formula, confidence score, SMILES, and PubChem IDs.
process_SIRIUS_table <- function(datadir) {
  # Read in SIRIUS table
  SIRIUS_raw <- read.delim(paste0(datadir, "SIRIUS/compound_identifications.tsv"), header = TRUE)
  SIRIUS_raw <- SIRIUS_raw %>%
    mutate(pubchemids = str_split(pubchemids, ";") %>% sapply(`[`, 1))
  
  # Reduce df to important columns
  KeepColumns <- c("id", "featureId", "molecularFormula", "name", "InChIkey2D", "ConfidenceScore", "ionMass", "retentionTimeInSeconds", "smiles", "pubchemids")
  SIRIUS <- SIRIUS_raw[, KeepColumns]
  
  # Parse "name" column and give empty column a name
  for(m in 1:dim(SIRIUS)[1]){
    if(SIRIUS$name[m] == "null"){
      SIRIUS$name[m] <- paste0("InChIKey: ", SIRIUS$InChIkey2D[m])
    }
  }
  
  # Reorder df for export
  neworder <- c("featureId", "name", "retentionTimeInSeconds", "ionMass", "molecularFormula", "ConfidenceScore", "smiles", "pubchemids")
  newnames <- c("id", "Annotation", "RT", "MZ", "Formula", "Confidence", "SMILES", "CID")
  SIRIUS <- SIRIUS[, neworder]
  colnames(SIRIUS) <- newnames
  SIRIUS$RT <- round(SIRIUS$RT / 60, 3)
  SIRIUS$Confidence <- as.numeric(SIRIUS$Confidence)
  
  # Use webchem package for better annotation
  pc.properties <- c("Title", "XLogP", "IUPACName")
  pc.info <- tryCatch({
    pc_prop(SIRIUS$CID, properties = pc.properties)
  }, error = function(e) {
    warning("Failed to retrieve data: ", e$message)
    return(data.frame(CID = SIRIUS$CID, Title = NA, IUPACName = NA))
  }) %>%
    mutate(Title = ifelse(is.na(Title), IUPACName, Title)) %>%
    select(-IUPACName) %>%
    dplyr::rename(Name = Title)
  
  pc.info$CID <- as.double(pc.info$CID)
  
  SIRIUS <- SIRIUS %>%
    mutate(CID = as.numeric(CID)) %>%
    left_join(pc.info, by = "CID", multiple = "first") %>%
    mutate(Annotation = Name) %>%
    select(-Name)
  
  return(SIRIUS)
}



# Function: Extract targeted features from raw data
# Input: data (raw data), df_ids (dataframe with IDs), inputFile (CSV file with features)
# Output: Dataframe with extracted features
extract_feature_list <- function(data, df_ids, inputFile, ppm_range = 3, rt_range = 0.07) {
  
  # Read the input CSV file
  TL <- read.csv(inputFile, header = TRUE)
  if (ncol(TL) < 2) {
    TL <- read.csv2(inputFile, header = TRUE)
  }
  
  # check if mz is already in the list otherwise calculate mz from Molecular Formula
  if (any(tolower(colnames(TL)) %in% c("m/z", "mz", "m.z"))) {
    colnames(TL)[which(tolower(colnames(TL)) %in% c("m/z", "mz", "m.z"))] <- "MZ"
    
  } else {
    # Convert Formula to m/z depending on specified adduct
    electron_mass <- 0.000548579909
    adduct_types <- unique(TL$Precursor.Adduct)
    for(add in adduct_types){
      if(add == "[M]+"){
        TL[,add] <- as.numeric(formula2mz(TL$Molecular.Formula, adduct = add)) - electron_mass
      }else{
        TL[,add] <- as.numeric(formula2mz(TL$Molecular.Formula, adduct = add))
      }
    }
    TL <- TL %>%
      pivot_longer(cols = all_of(adduct_types), names_to = "AddType", values_to = "tempMZ") %>%
      mutate(MZ = ifelse(Precursor.Adduct == AddType, tempMZ, NA)) %>%
      filter(!is.na(MZ)) %>%
      select(-tempMZ,-AddType)
  }
  
  output_df <- data.frame(matrix(nrow = 0, ncol = length(df_ids)+1))
  
  for(x in 1:dim(TL)[1]){
    # Define values of the current molecule in loop
    mol <- as.character(TL[x, "Molecule.Name"])
    RealMZ <- as.numeric(TL[x, "MZ"])
    RealRT <- as.numeric(TL[x, "Explicit.Retention.Time"])
    
    # 1. find all mzs within 5 ppm of standard mz
    featuredata.ppm <- df_ids %>%
      mutate(ppm = ((RealMZ - mz) / RealMZ) * 10**6,
             Molecule.Name = mol) %>%
      filter(ppm > -ppm_range, ppm < ppm_range) %>%
      distinct(id, .keep_all = TRUE) %>%
      select(-ppm)
    
    # 2. compare RTs and choose the one feature with least deviation to true RT, threshold of 0.1 max. deviation
    featuredata.rt <- featuredata.ppm %>%
      mutate(rt.abs.dev = abs(rt - RealRT)) %>%
      filter(rt.abs.dev == min(rt.abs.dev, na.rm = TRUE)) %>%
      mutate(id = ifelse(rt.abs.dev > rt_range, NA, id)) %>%
      select(-rt.abs.dev)
    
    # 3. add to output data frame
    output_df <- rbind(output_df, featuredata.rt)
  }
  
  output_df <- output_df %>%
    filter(!is.na(id)) %>%
    mutate(id = as.character(id))
  
  data_df <- data[output_df$id, ]
  data_df <- cbind(output_df[, 1:4], data_df)
  data_df$Molecule.Name <- output_df$Molecule.Name
  data_df <- data_df[!str_detect(rownames(data_df), "NA"),]
  
  return(data_df)
}


extract_feature_list_advanced <- function(data, df_ids, inputFile, ppm_range = 3, rt_range = 0.07) {
  
  # Read the input CSV file
  TL <- read.csv(inputFile, header = TRUE)
  if (ncol(TL) < 2) {
    TL <- read.csv2(inputFile, header = TRUE)
  }
  
  # check if mz is already in the list otherwise calculate mz from Molecular Formula
  if (any(tolower(colnames(TL)) %in% c("m/z", "mz", "m.z"))) {
    colnames(TL)[which(tolower(colnames(TL)) %in% c("m/z", "mz", "m.z"))] <- "MZ"
    
  } else {
    # Convert Formula to m/z depending on specified adduct
    electron_mass <- 0.000548579909
    adduct_types <- unique(TL$Precursor.Adduct)
    for(add in adduct_types){
      if(add == "[M]+"){
        TL[,add] <- as.numeric(formula2mz(TL$Molecular.Formula, adduct = add)) - electron_mass
      }else{
        TL[,add] <- as.numeric(formula2mz(TL$Molecular.Formula, adduct = add))
      }
    }
    TL <- TL %>%
      pivot_longer(cols = all_of(adduct_types), names_to = "AddType", values_to = "tempMZ") %>%
      mutate(MZ = ifelse(Precursor.Adduct == AddType, tempMZ, NA)) %>%
      filter(!is.na(MZ)) %>%
      select(-tempMZ,-AddType)
  }
  
  output_df <- data.frame(matrix(nrow = 0, ncol = length(df_ids)+1))
  
  for(x in 1:dim(TL)[1]){
    # Define values of the current molecule in loop
    mol <- as.character(TL[x, "Molecule.Name"])
    MZ_x <- as.numeric(TL[x, "MZ"])
    RT <- as.numeric(TL[x, "Explicit.Retention.Time"])
    
    # 1. find all mzs within 5 ppm of standard mz
    featuredata.ppm <- df_ids %>%
      mutate(ppm = ((MZ_x - mz) / MZ_x) * 10**6,
             Molecule.Name = mol) %>%
      filter(ppm > -ppm_range, ppm < ppm_range) %>%
      select(-ppm)
    
    # 2. compare RTs and give back all features below the RT threshold
    featuredata.rt <- featuredata.ppm %>%
      distinct(id, .keep_all = TRUE) %>%
      mutate(rt.abs.dev = abs(rt - RT)) %>%
      filter(rt.abs.dev < rt_range) %>%
      select(-rt.abs.dev) %>%
      filter(!is.na(id)) %>%
      mutate(id = as.character(id))
    
    # 3. Extract the relevant data
    data_df <- data[featuredata.rt$id, ]
    data_df <- cbind(featuredata.rt[, 1:4], data_df)  # Add metadata
    data_df$Molecule.Name <- featuredata.rt$Molecule.Name
    
    # 4. Handle cases with multiple features having the same MZ (±ppm_range) in TF
    similar_mz_features <- TL %>%
      filter((abs(MZ_x - TL$MZ) / MZ_x * 10^6) < ppm_range) %>%
      filter(abs(Explicit.Retention.Time - RT) < rt_range * 1.2) %>%
      arrange(Explicit.Retention.Time)  # Order by RT
    
    if (nrow(similar_mz_features) > 1) {
      # 5. Sort features by RT
      data_df <- data_df %>% arrange(rt)
      
      nr_in_seq <- which(similar_mz_features$Molecule.Name == mol)
      nr_similar_features <- nrow(similar_mz_features)
      
      # 6. Assign features correctly in data_df
      # Keep only the correct feature for each sample
      for (col in colnames(data_df)[5:(ncol(data_df) - 1)]) {
        
        non_na_indices <- which(!is.na(data_df[[col]]))  # Find non-NA values
        
        if (length(non_na_indices) == 0) next  # Skip if there are no values in this column
        
        if (length(non_na_indices) >= nr_similar_features) {
          # Keep only the highest values by RT order
          sorted_indices <- non_na_indices[order(data_df[[col]][non_na_indices], decreasing = TRUE)]
          keep_indices <- order(sorted_indices[1:nr_similar_features], decreasing = TRUE)  # Keep the top values but sort in the original order
          data_df[[col]][setdiff(non_na_indices, keep_indices)] <- NA  # Set the rest to NA
          
          # Now select the value corresponding to nr_in_seq
          data_df[[col]][setdiff(keep_indices, keep_indices[nr_in_seq])] <- NA  # Keep only the selected row and set the others to NA
        } else {
          # If there are fewer non-NA values than needed, choose the one closest to RT
          closest_idx <- which.min(abs(data_df$rt[non_na_indices] - RT))
          data_df[[col]][setdiff(non_na_indices, non_na_indices[closest_idx])] <- NA # Keep only the selected row and set the others to NA
        }
      }
      
    } else {
      # No RT conflicts, just pick the highest intensity feature per sample
      for (col in colnames(data_df)[5:(ncol(data_df)-1)]) {
        max_idx <- which.max(data_df[[col]])
        data_df[[col]][-max_idx] <- NA
      }
    }
    
    # 7. Sum up and remove duplicate rows and keep Molecule.Name column as last column
    data_df <- data_df %>%
      group_by(Molecule.Name) %>%
      summarise(across(everything(), ~ ifelse(all(is.na(.)), NA, max(., na.rm = TRUE))), .groups = "drop") %>%
      relocate(Molecule.Name, .after = last_col())  # Move Molecule.Name to last column
    
    # 8. Append to output
    output_df <- rbind(output_df, data_df)
  }
  
  rownames(output_df) <-  make.unique(as.character(output_df$id))
  
  return(output_df)
}


# Function: Process raw file and return data frames for further steps
# Input: datadir (directory containing raw data, already prepared csv file from mzMine, includes hpos, npos and neg)
# Output: List containing 'data' (processed data) and 'df_ids' (IDs dataframe)
process_merged_files <- function(datadir) {
  file_merged <- paste0(datadir, "areas.csv")
  
  raw <- if (file.info(file_merged)$size < 2^31 - 1) {
    # Normal fread for files < 2 GB
    data.table::fread(file_merged, data.table = FALSE)
  } else {
    print("Area file is too big, need a different reading approach, this could take some time.")
    
    # Chunked reading for large files, with temporary files
    (function(file, chunk_lines = 100000) {
      chunks <- list()
      con <- file(file, open = "r")
      on.exit(close(con))
      
      header <- readLines(con, n = 1)
      col_names <- strsplit(header, ",")[[1]]
      
      repeat {
        lines <- readLines(con, n = chunk_lines)
        if (length(lines) == 0) break
        
        tmp <- tempfile()
        writeLines(c(header, lines), tmp)
        
        chunk <- read_csv(tmp, col_names = col_names, show_col_types = FALSE)
        chunks[[length(chunks) + 1]] <- chunk
        
        unlink(tmp)
      }
      
      bind_rows(chunks)
    })(file_merged)
  }

  area_columns <- sort(grep("^datafile.*area$", names(raw), value = TRUE))
  df_ids <- raw[,c("id", "rt", "mz")]
  df_areas <- raw[,c("id", area_columns)] %>% column_to_rownames("id")
  newcolnames <- sub("^datafile.(.*).*.mzML.area$", "\\1", names(df_areas))
  newcolnames <- gsub("_1$", "", newcolnames) # Remove the "_1" suffix (used in some old file names, mostly not relevant)
  newcolnames <- gsub("^[0-9]{6}_|FiS[0-9]+_", "", newcolnames) # Remove the date and "FiS4/5/6/7" prefix (used in some new file names so the hpos and npos from the same sample get put together)
  colnames(df_areas) <- newcolnames
  
  # This collapses the columns belonging to the same sample with different measurement conditions into a single column to reduce dimensionality
  # Requires there to be only ONE valid value across all measurement conditions per row and sample!!
  measurement_conditions <- unique(sub(".*_", "_", newcolnames))   # Extract the part of the name after the last underscore to get the measurement conditions
  gsub.pattern <- paste0("(", paste(measurement_conditions, collapse = "|"), ")")

  data <- as.data.table(df_areas, keep.rownames = "ID") # Convert to data.table
  data <- melt(data, id.vars = "ID", variable.name = "Sample", value.name = "Area") # Reshape from wide to long format
  data[, Sample := gsub(gsub.pattern, "", Sample)] # Apply gsub replacement in one step
  data[, Area := as.numeric(Area)]
  data <- data[, .(Area = sum(Area, na.rm = TRUE)), by = .(ID, Sample)] # Aggregate to sum areas per ID and Sample
  data[, Area := fifelse(Area == 0, NA_real_, Area)] # Replace 0 with NA
  data <- dcast(data, ID ~ Sample, value.var = "Area") # Reshape from long to wide format
  
  # Convert back to rownames
  setDF(data) # Convert back to data.frame if needed
  rownames(data) <- data$ID
  data$ID <- NULL
  
  # Check for ProcBlank columns, meaning it was not filterd in mzMine, put all procblanks in one col., filter the features of the samples
  if (any(grepl("ProcBlank", colnames(data)))) {
    # Use this approach with sweep instead of rowwise of dplyr as this is way faster
    procblank_cols <- grep("ProcBlank", colnames(data), value = TRUE)
    
    # Calculate average ProcBlank per row
    procblank_avg <- rowMeans(select(data, all_of(procblank_cols)), na.rm = TRUE)
    
    # Identify sample columns to apply filtering
    sample_cols <- colnames(data)[
      sapply(data, is.numeric) &
        !grepl("Plasma", colnames(data)) &
        !colnames(data) %in% procblank_cols
    ]
    
    # Create a matrix of the data to filter
    sample_matrix <- as.matrix(data[, sample_cols])
    
    # Apply vectorized filtering
    filter_mask <- sweep(sample_matrix, 1, 3 * procblank_avg, `<`)
    sample_matrix[filter_mask & !is.na(sample_matrix)] <- NA
    
    # Replace filtered columns in original data
    data[, sample_cols] <- sample_matrix
    
    # Remove ProcBlank columns
    data <- data %>% select(-all_of(procblank_cols))
  }
  
  
  Polarity <- as.data.table(df_areas, keep.rownames = "id") # Convert to data.table
  Polarity <- melt(Polarity, id.vars = "id", variable.name = "Sample", value.name = "Area") # Reshape to long format
  Polarity <- Polarity[!is.na(Area)] # Remove NA values
  Polarity[, charge := gsub(".*.([a-z]{3})$", "\\1", Sample)] # Extract charge using regex
  Polarity <- unique(Polarity[, .(id, charge)]) # Keep distinct id-charge pairs
  Polarity[, charge := fifelse(charge == "pos", 1, -1)] # Convert charge to numeric (1 for pos, -1 for others)
  df_ids <- merge(df_ids, Polarity, by = "id", all.x = TRUE) # Merge with df_ids, keeping all original rows
  df_ids$rt <- as.numeric(df_ids$rt)
  df_ids$mz <- as.numeric(df_ids$mz)

  return(list(data = data, df_ids = df_ids))
}


# Use PubChem and ChatGPT to assist in assining common names instead of complicated chemical UPAC names
update_annotations_ai_assisted <- function(df) {
  openai_api_key <- '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'
  
  # Function to query PubChem and get common names and synonymes
  get_pubchem_name <- function(name, smiles) {
    # Store potential identifiers
    identifiers <- list(name = name, smiles = smiles)
    
    # Remove NULL or empty identifiers
    identifiers <- identifiers[!is.na(identifiers) & nchar(identifiers) > 0]
    
    if (length(identifiers) == 0) {
      return(NA)  # No valid identifiers
    }
    
    # Query PubChem for all identifiers
    cid_list <- lapply(names(identifiers), function(source) {
      get_cid(identifiers[[source]], from = source)
    })
    
    # Flatten and remove NAs
    cid_list <- unlist(cid_list)
    cid_list <- cid_list[!is.na(cid_list)]
    
    if (length(cid_list) == 0) {
      return(NA)  # No results found
    }
    
    # Get properties from PubChem
    pubchem_data <- pc_prop(unique(cid_list), properties = c("IUPACName", "Title", "Synonyms"))
    
    if (is.na(pubchem_data)) {
      return(NA)
    }
    
    # Collect all names
    all_names <- c(
      if (!is.null(pubchem_data$Title[[1]])) pubchem_data$Title[[1]] else NULL,
      if (!is.null(pubchem_data$Synonyms[[1]])) pubchem_data$Synonyms[[1]] else NULL,
      if (!is.null(pubchem_data$IUPACName[[1]])) pubchem_data$IUPACName[[1]] else NULL
    )
    
    if (length(all_names) == 0) {
      return(NA)  # No valid names found
    }
    
    return(all_names)  # Return all names for OpenAI processing
  }
  
  # Function to query ChatGPT and return the best Annotation for the feature
  get_chatgpt_name <- function(original_name, pubchem_names, openai_api_key) {
    if (is.na(pubchem_names) || length(pubchem_names) == 0) {
      return(original_name)  # Keep original if no PubChem results
    }
    
    # Format PubChem names properly
    pubchem_list <- paste(pubchem_names, collapse = ", ")
    
    # Structured prompt
    prompt <- paste(
      "You are an expert in metabolomics and chemical nomenclature. Your task is to determine the most commonly used name for a metabolite in biological and chemical studies.",
      "\n\n### Provided Information:",
      "\n- **Original Name:**", original_name,
      "\n- **PubChem Suggested Names:**", pubchem_list,
      "\n\n### Instructions:",
      "1. If one of the PubChem names is already the best and most commonly used in metabolomics, return it as is.",
      "2. If another name is more commonly used, return that name instead.",
      "3. If unsure or no better name exists, return the original name.",
      "4. For dipeptides and tripeptides, always return the name in the format 'Amino1-Amino2' (e.g., Tyr-Gly, Gly-Ala-Ser), unless another format is significantly more common.",
      "5. Avoid long IUPAC names unless they are the standard name used in metabolomics. For example, return 'citric acid' instead of '2-hydroxy-1,2,3-propane-tricarboxylic acid'.",
      "6. Prefer commonly used names over overly detailed systematic names. For example, use 'choline' instead of '2-hydroxy-N,N,N-trimethylethanaminium'.",
      "7. If a compound is a well-known drug or metabolite, do not change its widely accepted name. For example, keep 'Acetaminophen' instead of 'N-(4-hydroxyphenyl)acetamide'.",
      "8. If multiple names are equally valid, choose the one most frequently used in metabolomics literature."
    )
    
    # API request
    response <- POST(
      url = "https://api.openai.com/v1/chat/completions",
      add_headers(Authorization = paste("Bearer", openai_api_key)),
      content_type_json(),
      body = toJSON(list(
        model = "gpt-4-turbo",
        messages = list(
          list(role = "system", content = "You are an expert in chemical nomenclature, particularly in metabolomics."),
          list(role = "user", content = prompt)
        ),
        temperature = 0.3
      ), auto_unbox = TRUE)
    )
    
    # Process API response
    chatgpt_response <- content(response, as = "parsed", simplifyVector = TRUE)
    if (!is.null(chatgpt_response$choices) && length(chatgpt_response$choices) > 0) {
      return(chatgpt_response$choices[[1]]$message$content)
    } else {
      return(original_name)  # Keep original if ChatGPT fails
    }
  }
  
  # Apply functions to the dataframe
  df <- df %>%
    rowwise() %>%
    mutate(
      Refined_Annotation = get_chatgpt_name(Annotation, get_pubchem_name(Annotation, Formula, SIMILES, mz), Formula, SMILES, mz, openai_api_key)
    ) %>%
    ungroup()
  
  # Replace old Annotation column
  df <- df %>%
    select(-Annotation) %>%
    rename(Annotation = Refined_Annotation)
  
  return(df)
}




# Function: Process data with SIRIUS annotations to obtain significant and abundant features
# Input: datadir (directory containing raw data and SIRIUS results)
# Output: List containing 'se' (SummarizedExperiment object), 'significant_features' (data frame of significant features),
#         'significant_abundant_features' (data frame of significant and abundant features), 'mean_abundance_each_col' (mean abundance of each column)get_significant_abundant_features <- function(datadir) {
get_significant_abundant_features <- function(datadir, exp_data, info_file_dir, resultsdir, figuredir, exp, tlists, normalization_list = NULL, remove_outliers = TRUE) {
  if (!dir.exists(resultsdir)) {
    dir.create(resultsdir, recursive = TRUE)
  }
  
  # Process SIRIUS table
  SIRIUS <- process_SIRIUS_table(datadir)
  
  # Get the whole processed data set
  all_data <- exp_data$data
  df_ids <- exp_data$df_ids
  df_data <- all_data # could be used to already select here subdataframe
  

  # Merge feature data from Sirius with data ids
  df_features <- merge(df_ids, SIRIUS, by = "id", all.x = TRUE)
  df_features <- df_features[!duplicated(df_features$id), ] # remove duplicate annotated features
  
  # Read in the metadata/phenodata from file 
  add_info <- read.csv(info_file_dir)
  if (length(add_info) < 2) {
    add_info <- read.csv2(info_file_dir)
  }
    
  # Create phenodata dataframe from imported metadata, check if the 'Intervention' column exists
  if ("Sample" %in% names(add_info)) {
    if ("Intervention" %in% names(add_info)) {
      df_phenodata <- add_info %>%
        mutate(Intervention = factor(Intervention, levels = unique(Intervention)),
               Donor = factor(Donor, levels = unique(Donor))) %>%
        column_to_rownames("Sample")
    } else {
      df_phenodata <- add_info %>%
        mutate(Intervention = case_when(
          Timepoint == 't0' ~ "pre_Intervention",
          Timepoint == 'T0' ~ "pre_Intervention",
          Timepoint == '0h' ~ "pre_Intervention",
          Timepoint == '0' ~ "pre_Intervention",
          Timepoint != 't0' & Timepoint != 'T0' & Timepoint != '0h' & Timepoint != '0' ~ "post_Intervention",
          is.na(Timepoint) ~ "post_Intervention"
        )) %>%
        mutate(Intervention = factor(Intervention, levels = c("pre_Intervention", "post_Intervention")),
               Donor = factor(Donor, levels = unique(Donor))) %>%
        column_to_rownames("Sample")
    }
  } else {
    # if the metadata file doesnt contain sample names, we need to create them ourself, we need them to create the SummarizedExperiment
    df_phenodata <- data.frame(Sample = colnames(df_data), Intervention = rep("post_Intervention", length(colnames(df_data)))) %>%
      mutate(Intervention = factor(Intervention, levels = c("post_Intervention"))) %>%
      column_to_rownames("Sample") %>%
      rename("Intervention")
  }
  
  # Standardize the rownames and colnames to ensure compatibility
  normalize_names <- function(names) {
    # Remove 6-digit numbers followed by an underscore (if present)
    gsub("^[0-9]{6}_", "", names)
  }
  
  # Apply normalization to rownames of df_phenodata and colnames of df_data
  rownames(df_phenodata) <- normalize_names(rownames(df_phenodata))
  colnames(df_data) <- normalize_names(colnames(df_data))
  
  
  # Ensure rownames of df_phenodata are identical to colnames of df_data by reordering df_phenodata
  df_data <- df_data[, colnames(df_data) %in% rownames(df_phenodata)]
  df_phenodata <- df_phenodata[match(colnames(df_data), rownames(df_phenodata)), ]
  

  # SummarizedExperiment: Assembly
  se <- SummarizedExperiment(assays = list("raw" = df_data), 
                             rowData = df_features,
                             colData = df_phenodata)
  
  # Feature Reduction to reduce noise
  se <- removeFeatures_advanced(se, i = "raw", group = se$Intervention, cut = 0.7, min_groups = 1, method = "missing") # remove features that have less than 70% valid values in at least on group/condition
  
  # find outlier samples and remove them before the normalization steps (currently only show which is an outlier and dont remove)
  raw_assay <- as.data.frame(assay(se, "raw"))
  modified_assay <- check_samples_for_outliers_with_powerfunction(raw_assay, resultsdir, figuredir, exp)
  trunc_df <- modified_assay$truncated_df
  assay(se, i = "cut_of_features") <- trunc_df
  
  #removed_outlier_samples <- modified_assay$outlier_sample_names
  # Remove outlier samples by checking QC molecules
  removed_outlier_samples <- find_outliers_with_QC_molecules(df_data, df_ids, df_features, df_phenodata, tlists, resultsdir)
  if (remove_outliers) {
    se <- se[, !(colnames(se) %in% removed_outlier_samples)] 
  }
  se <- removeFeatures_advanced(se, i = "cut_of_features", group = se$Intervention, cut = 0.7, min_groups = 1, method = "missing") # remove features that have less than 70% valid values in at least on group/condition
  
  # Delete duplicate annotations by identifying most data dense and intense features
  clean <- extract.non.duplicates.vector(se, 'cut_of_features')
  se <- se[rownames(se) %in% as.character(clean),] 
  
  # Normalization (Sample size error correction)
  if (!is.null(normalization_list)) {
    assay(se, i = "cut_of_features") <- normalize_data(assay(se, i = "cut_of_features"), normalization_list, df_f)[,-c(1:9)]
  }
    
  # Imputation, Normalization and Clustering
  assay(se, i = "log2") <- assay(se, i = "cut_of_features") %>% mutate_all(~log2(.) + 20)
  
  # vsn with this, currently not used
  se <- normalizeIntensity(se, i = "cut_of_features", name = "vsn", method = "vsn") # should do the same as normalizePQN
  # pqn with this, currently not used
  se <- normalizeIntensity(se, i = "log2", name = "log2_pqn", method = "pqn") # should do the same as normalizePQN
  # introduce final sample size normalization HERE
  se <- imputeIntensity(se, i = "log2_pqn", name = "log2_pqn_MinProb", method = "MinProb")
  se <- clusterFeatures(se, i = "log2_pqn_MinProb", rtime_var = "rt", rt_cut = 0.003)
  # Imputing without pqn or vsn, currently used
  se <- imputeIntensity(se, i = "log2", name = "log2_MinProb", method = "MinProb")
  se <- clusterFeatures(se, i = "log2_MinProb", rtime_var = "rt", rt_cut = 0.003)
  assay(se, i = "MinProb") <- 2 ** (assay(se, i = "log2_MinProb") - 20)
  
  # Get significant features
  feature_groups <- as.data.frame(rowData(se)[, c("id", "rt", "mz", "rtime_group", "feature_group", "Annotation", "Formula", "Confidence", "SMILES")])
  #significant_features <- feature_groups[!is.na(feature_groups$Confidence) & feature_groups$Confidence > 0.25, ]
  significant_features <- feature_groups

  
  
  # check which rows and columns to append, only use the rows of significant features and no outliers
  #columns_to_append <- df_data[significant_features$id, !(colnames(df_data) %in% removed_outlier_samples)] 
  #columns_to_append <- assay(se, i = "log2_pqn_MinProb")[significant_features$id, !(colnames(df_data) %in% removed_outlier_samples)] 
  columns_to_append <- assay(se, i = "MinProb")
  all_significant_feature_data <- cbind(significant_features, columns_to_append)

  
  
  # additional normalization to adjust for size effects (can already normalize here and dont need to do it in the figures)
    #all_significant_feature_data <- normalize_data(all_significant_feature_data, normalization_list)
  
  
  # try to see if we can find samples that have a strong PEG contamination peak (big one at 4.5 to 5.1 min)
  # filter_for_ethylenglycol <- all_significant_feature_data[grepl("ethylene glycol", all_significant_feature_data$Annotation, ignore.case = TRUE), ]
  # sumup_ethylenglycol <- colSums(filter_for_ethylenglycol[, 10:ncol(filter_for_ethylenglycol)], na.rm = TRUE)
  # filter_futher_for_ethylenglycol <- filter_for_ethylenglycol %>% filter(rt > 4.6) %>% filter(rt < 5.2)
  # filter_featrues_for_ethylenglycol <- all_significant_feature_data %>% filter(rt > 4.7) %>% filter(rt < 5.1)
  # sumup_ethylenglycol <- colSums(filter_featrues_for_ethylenglycol[, 10:ncol(filter_for_ethylenglycol)], na.rm = TRUE)
  # 
  # all_features_with_charge <- merge(all_significant_feature_data, raw[, 1:6], by = "id", all.x = TRUE)
  # all_features_with_charge$charge[is.na(all_features_with_charge$charge)] <- 2
  # filter_featrues_for_ethylenglycol_by_charge <- all_features_with_charge %>% 
  #                                                 filter(rt.x > 4.9) %>% filter(rt.x < 5.1) %>% 
  #                                                 filter(charge == 2) %>% filter(is.na(Annotation)) %>%
  #                                                 filter(mz > 400) #%>% filter(mz < 800)
  # sumup_ethylenglycol <- colSums(filter_featrues_for_ethylenglycol_by_charge[, 10:ncol(filter_featrues_for_ethylenglycol_by_charge)], na.rm = TRUE)
  
  
  data_length <- ncol(all_significant_feature_data)
  #all_significant_feature_data[, 10:data_length][is.na(all_significant_feature_data[, 10:data_length])] <- 1 #make Na 1
  #num_cols_abundant <- rowSums(all_significant_feature_data[, 10:data_length] > 1e5, na.rm = TRUE)
  #significant_abundant_features <- all_significant_feature_data[num_cols_abundant >= (data_length * 0.1), ] # features must be abundant in at least 10% of samples to still get features of control samples
  significant_abundant_features <- all_significant_feature_data
  
  mean_abundance_each_col <- colMeans(significant_abundant_features[, 10:data_length], na.rm = TRUE)
  
  return(list(se = se, significant_features = significant_features, significant_abundant_features = significant_abundant_features, 
              mean_abundance_each_col = mean_abundance_each_col, all_significant_feature_data = all_significant_feature_data))
}




########
# section with various small functions used in other functions
# for example: normalize data, calculate stats (mean, std dev), check for outliers
#              help with figure plotting (size adjustments) and figure saving
#              prepare data for plots



# Function to normalize the areas of each column with the sum of all areas of the according column
normalize_data <- function(df, normalization_list = NULL, df_f = NULL) {
  if (!'Annotation' %in% colnames(df)) {
    if (!is.null(df_f)) {
      df <- merge(df_f[,1:9], df, by=0) %>% mutate(Row.names = as.numeric(Row.names)) %>% arrange(Row.names) %>% column_to_rownames("Row.names")
    } else {
      print('No Annotation column found and no Feature dataframe given, please provide necessary data.')
      stop()
    }
  } 
  
  df_normalized <- df
  sample_names <- colnames(df)[10:ncol(df)]
  norm_factors <- numeric(length(sample_names))
  
  # if no molecule list for normalization is provided use all features for normalization
  if (is.null(normalization_list)) {
    for (col_index in 10:ncol(df)) {
      col_sum <- sum(df[, col_index], na.rm = TRUE)
      df_normalized[, col_index] <- (df[, col_index] / col_sum)
      colnames(df_normalized)[col_index] <- paste0(colnames(df)[col_index], "_norm")
      norm_factors[col_index - 9] <- col_sum
    }
  } else {  # use the molecules from the list for normalization of each sample
    normalization_molecules_all_data <- read.csv(normalization_list)
    if (! "Molecule.Name" %in% colnames(normalization_molecules_all_data)) { #check for Molecule.Name as otherwise it is the wrong delimiter need different read in
      normalization_molecules_all_data <- read.csv2(normalization_list)
    }
    normalization_molecules <- as.vector(normalization_molecules_all_data$Molecule.Name)
    
    # Loop over the column indices from 10 to the last column
    for (col_index in 10:ncol(df)) {
      # Initialize col_sum to zero
      col_sum <- 0
      
      # Loop over the names in the normalization_list
      for (i in seq_along(normalization_molecules)) {
        # Sum the values in column with index col_index for matching rows
        col_sum <- col_sum + sum(df[df$Annotation == normalization_molecules[i], col_index], na.rm = TRUE)
      }
      
      # Normalize the values in the column with index col_index
      df_normalized[, col_index] <- (df[, col_index] / col_sum)
      
      # Update column names of df_normalized
      #colnames(df_normalized)[col_index] <- paste0(colnames(df)[col_index], "_norm")
      
      # Store the col_sum
      norm_factors[col_index - 9] <- col_sum
    }
  }
  
  # Create df_norm_sums as a dataframe
  df_norm_sums <- data.frame(ID = sample_names, NormFactor = norm_factors)

  write.csv(df_norm_sums, file= paste0(resultsdir, "sample_normalization_sums.csv"), row.names = FALSE )
  
  return(df_normalized)
}



# Calculate mean, std dist, and rank
calculate_stats <- function(df) {
  df$mean <- rowMeans(df[, 10:ncol(df)], na.rm = TRUE)
  df$std_dev <- apply(df[, 10:(ncol(df)-1)], 1, sd, na.rm = TRUE)
  df$rel_std_dev <- (df$std_dev / df$mean) * 100
  df$mean_rank <- rank(-df$mean)
  return(df)
}



# Function to find duplicates and give back a list of the rows we want to keep
extract.non.duplicates.vector <- function(sumex, assay_name){
  # Extract only necessary columns from features and combine with areas
  features <- as.data.frame(rowData(sumex)) %>% select(Annotation, id)
  areas <- as.data.frame(assay(sumex, i = assay_name))
  combined <- cbind(features, areas)
  
  # Identify non-unique Annotations
  non_unique <- unique(na.omit(combined$Annotation[duplicated(combined$Annotation) | duplicated(combined$Annotation, fromLast = TRUE)]))
  
  # Separate clean and non-unique annotations
  clean_features <- combined %>%
    filter(!Annotation %in% non_unique) %>%
    pull(id)
  
  # Process non-unique annotations
  stats <- combined %>%
    filter(Annotation %in% non_unique) %>%
    pivot_longer(cols = -c(Annotation, id), names_to = "variable", values_to = "value") %>%
    group_by(Annotation, id) %>%
    dplyr::summarize(
      data.density = round(sum(!is.na(value)) / n(), 4) * 100,
      avg = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(Annotation) %>%
    mutate(max_density = max(data.density)) %>%
    filter(
      data.density > max_density - 15,
      avg == max(avg)
    ) %>%
    ungroup()
  
  # Combine clean and filtered IDs
  rows_to_keep <- sort(c(clean_features, stats$id))
  
  return(rows_to_keep)
}



# Function to find potential outlier with QC molecules by comparing them within the data set
find_outliers_with_QC_molecules <- function(df_a, df_ids, df_f, df_p, tlists, resultsdir) {
  # use the usual screening transition list to find QC molecules and compare them within the data set
  inputFile <- paste0(tlists,"Screening_TransitionList.csv")
  Screen <- extract_feature_list(df_a, df_ids, inputFile)
  Screen <- Screen %>%
    distinct(id, .keep_all = TRUE)
  
  # define QC molecules
  AAs <- c("Arginine", "Asparagin", "Aspartic Acid", "Glutamic Acid", "Glutamine", "Proline", "Phenylalanine",
           "Pyroglutamic acid", "Serine", "Threonine", "Tyrosine", "Valine", "Tryptophan", "Lysine", "Methionine")
  
  # pivots the df to turn into plotable format
  AA_plots <- Screen[Screen$Molecule.Name %in% AAs,] %>%
    select(-rt,-mz,-charge) %>%
    left_join(df_f %>% mutate(id = as.character(id)), by = "id") %>%
    pivot_longer(cols = -c(colnames(df_f), Molecule.Name), names_to = "SampleID", values_to = "Area") %>%
    mutate(SampleID = factor(SampleID, levels = unique(SampleID)),
           Area = ifelse(is.na(Area), 1, Area),
           log2 = log2(Area)) %>%
    left_join(df_p %>% rownames_to_column("SampleID"), join_by("SampleID"))
  
  # uses quantile calculations to identify samples which are outside the norm
  # assumes that the majority of the data is of decent quality!!
  AA_stats <- AA_plots %>%
    group_by(Molecule.Name) %>%
    mutate(Q1 = quantile(log2, 0.25, na.rm = TRUE),
           Q3 = quantile(log2, 0.75, na.rm = TRUE),
           IQR = Q3 - Q1,
           lb_IQR = Q1 - 1.5 * IQR,
           percentile_cutoff = quantile(log2, 0.05, na.rm = TRUE),
           final_lower_bound = pmax(lb_IQR, percentile_cutoff),
           outlier = log2 < final_lower_bound)
  
  # calculates how often a sample is identified as an outlier across all tested molecules
  # samples above 70% outlier rate are flagged
  outlier_summary <- AA_stats %>%
    group_by(SampleID) %>%
    summarise(total_molecules = n(),
              outlier_count = sum(outlier),
              outlier_percentage = outlier_count / total_molecules) %>%
    filter(outlier_percentage > 0.7)
  outliers <- outlier_summary$SampleID
  
  print("Samples identified and exluded as low quality:\n")
  print(outliers)
  
  df_a <- df_a[,!colnames(df_a) %in% outliers]
  df_p <- df_p[colnames(df_a),]
  
  # Adjust label size dynamically
  sample_count <- length(unique(AA_plots$SampleID))
  label_size <- ifelse(sample_count < 50, 12, ifelse(sample_count < 100, 8, 5))
  
  
  # plots and saves amino acid area value profiles across all samples
  ggplot(AA_plots, aes(x = SampleID, y = Area, color = Molecule.Name)) +
    geom_line(aes(group = Molecule.Name)) +
    geom_vline(data = outlier_summary, aes(xintercept = SampleID), alpha = 0.5, color = "red") +
    theme_minimal() +
    theme(
      axis.text.x = element_markdown(angle = 90, vjust = 0.5, size = label_size)  # Dynamic font size, Enable colored text
    ) +
    scale_x_discrete(labels = function(x) {
      ifelse(
        x %in% outliers,
        sprintf("<span style='color:red;'>%s</span>", x),  # Red for outliers
        sprintf("<span style='color:black;'>%s</span>", x) # Black for others
      )
    }) +
    labs(title = "QC Molecule Overview for Outlier Detection")
  
  ggsave(paste0(resultsdir, "AminoAcids.png"), width = 30, height = 10)
  
  return(outliers)
}



# Function to fit power function, cut off features, and identify outlier samples
check_samples_for_outliers_with_powerfunction <- function(df, resultsdir, figuredir, exp) {
  # Helper function to fit a power function and cut off features
  fit_and_cut <- function(x) {
    # Apply log2 transformation
    log_x <- log2(x)
    
    # Sort the log-transformed values and store the sorting indices
    sorted_indices <- order(-log_x, na.last = TRUE)
    sorted_log_x <- log_x[sorted_indices]
    
    # Replace NA values with 1 so cutoff works
    sorted_log_x_NAto1 <- sorted_log_x
    sorted_log_x_NAto1[is.na(sorted_log_x_NAto1)] <- 1
    
    # Fit power function y = a * x^b using non-linear least squares
    # nlsLM ignores NAs 
    #power_fit <- nlsLM(sorted_log_x ~ a * (1:length(sorted_log_x))^b, start = list(a = 1, b = 1)) # first simple formula, tested more complex one, adjusted some other following lines as well
    power_fit <- nlsLM(sorted_log_x ~ a * (1 + b * exp(-k * (1:length(sorted_log_x)))) * (1:length(sorted_log_x))**n, 
                       start = list(a = 2, b = 10, k = 4e-5, n = -5e-2), # choose the most often occuring values as start values
                       upper = c(a = 10, b = 100, k = 1e-2, n = -1e-5)) # n must be smaller than 0 otherwise the curve at the beginning won't fit, then have to give all parameters a limit
    #power_fit <- nlsLM(sorted_log_x ~ (a * (1 + b * exp(-k * (1:length(sorted_log_x)))) * (1:length(sorted_log_x))**n) + c, 
    #                   start = list(a = 1, b = 1, k = 0.1, n = 1, c = 0))
    #power_fit <- nlsLM(sorted_log_x ~ -a * log(b * (1:length(sorted_log_x)) + c), 
    #                   start = list(a = 1, b = 1, c = 0))
    
    # Get coefficients
    coefs <- coef(power_fit)
    a <- coefs['a']
    b <- coefs['b']
    k <- coefs['k']
    n <- coefs['n']
    #c <- coefs['c']
    
    # Calculate R-squared
    res_std_err <- summary(power_fit)$sigma
    
    # Calculate the gradient of the power function
    #power_gradient <- a * b * (1:length(sorted_log_x))^(b - 1)
    power_gradient <- a * (1:length(sorted_log_x))**(n-1) * (-b*k*(1:length(sorted_log_x)) + b*n + n*exp(k*(1:length(sorted_log_x)))) * exp(-k*(1:length(sorted_log_x)))

    # Calculate function values
    # predict() would we easier but it only predicted till the length of sorted_log_x without NAs
    #function_values <- a * (1:length(sorted_log_x))^b
    function_values <- (a * (1 + b * exp(-k * (1:length(sorted_log_x)))) * (1:length(sorted_log_x))**n)
    #function_values <- (a * (1 + b * exp(-k * (1:length(sorted_log_x)))) * (1:length(sorted_log_x))**n) + c
    
    suppressWarnings( # the lengths dont all fit so warnings are possible, suppress them 
    # Find cut-off index where gradient condition is met
    cutoff_index <- which.max((lead(sorted_log_x_NAto1) == 1) | # check if next value would be 1 (=NA), can prevent problems with other conditions
                                ((abs(diff(sorted_log_x_NAto1)) > 4 * power_gradient) &
                                (seq_along(sorted_log_x) > 0.5 * length(sorted_log_x)) & # for simple formula 0.1, for complex one need 0.5, otherwise cut off is to early
                                (sorted_log_x < function_values)))
    )
    
    # Mark cut-off features as NA if there are any to cut off
    if (cutoff_index < length(sorted_log_x)) {
      sorted_log_x[(cutoff_index):length(sorted_log_x)] <- NA
    }
    
    # Reorder the truncated feature vector back to the original order
    truncated_log_x <- rep(NA, length(x))
    truncated_log_x[sorted_indices] <- sorted_log_x
    
    # Convert log2-transformed values back to original scale
    original_ordered_truncated_x <- 2^truncated_log_x
    
    # Return the truncated feature vector with NA for cut-off values in original order
    return(list(original_ordered_truncated_x = as.data.frame(original_ordered_truncated_x), 
                cutoff_index = as.list(cutoff_index),
                res_std_err = as.list(res_std_err),
                function_values = function_values,
                coefs =  coefs))
  }

  
  # for bigger datasets parallelize the process to speed it up, create the cluster accordingly, 
  # too big clusters have more overhead so create dynamically
  if (ncol(df) > 20) {
    # if there are lots of samples than parallelize, check how many cores are available and how many are needed
    maxCores <- min(c(40, detectCores()-1))
    coresNeeded <-  ncol(df) %/% 5
    coresUsed <- min(c(coresNeeded, maxCores))
    
    cl <- makeCluster(coresUsed)
    
    clusterEvalQ(cl, {
      library(dplyr)
      library(minpack.lm)
      library(nlstools)
      library(tidyr)
      library(ggplot2)
    })
    
    clusterExport(cl, varlist = c('fit_and_cut'), envir = environment())
    
    # Apply the fit_and_cut function to each column
    results  <- parLapply(cl, df, fit_and_cut)
    stopCluster(cl)
    
  } else {
    # Apply the fit_and_cut function to each column
    results  <- lapply(df, fit_and_cut)
  }
  
  # Extract truncated data and cutoff indices
  truncated_df <- as.data.frame(lapply(results, function(res) res$original_ordered_truncated_x))
  cutoff_indices <- unlist(lapply(results, function(res) res$cutoff_index))
  res_std_err_values <- unlist(lapply(results, function(res) res$res_std_err))
  function_values_list <- as.data.frame(lapply(results, function(res) res$function_values))
  coefs <- unlist(lapply(results, function(res) res$coefs))
  
  colnames(truncated_df) <- colnames(df)
  rownames(truncated_df) <- rownames(df)
  
  
  # Calculate the number of remaining features after cutoff for each sample
  remaining_features <- colSums(!is.na(truncated_df))
  
  # Calculate the median of the remaining features
  median_features <- median(remaining_features)
  # Calculate the Median Absolute Deviation (MAD)
  mad <- median(abs(remaining_features - median_features))
  
  # Identify outlier samples using the modified z-score
  modified_z_scores <- 0.6745 * abs((remaining_features - median_features) / mad) # use modified z-score (uses median and more robust than normal z-score)
  outlier_samples <- which(modified_z_scores > 3.5)  # Threshold from literature
  outlier_sample_names <- names(df)[outlier_samples]
  
  # Generate plots
  plot_log2_with_power_overlay(df, outlier_samples, cutoff_indices, function_values_list, res_std_err_values, resultsdir, figuredir, exp)
  
  # Save remaining features and residual standard errors to csv
  df_features <- data.frame(Sample_name = names(remaining_features), remaining_features = remaining_features)
  write.csv(df_features, file = paste0(resultsdir,"remaining-features-for-pqn_each-sample.csv"), row.names = FALSE)
  df_res_std_err <- data.frame(Sample_name = names(res_std_err_values), res_std_err = res_std_err_values)
  write.csv(df_res_std_err, file = paste0(resultsdir,"res-std-err_each-sample.csv"), row.names = FALSE)
  write.csv(coefs, file = paste0(resultsdir,"coeffients-of-powerfunction-per-sample.csv"))
  write.csv(outlier_sample_names, file = paste0(resultsdir,"outlier-samples.csv"), row.names = FALSE)
  
  # Return the original dataframe, the truncated dataframe, and the list of outlier samples
  list(original_df = df, truncated_df = truncated_df, remaining_features_list = remaining_features, 
       outlier_samples = outlier_samples, outlier_sample_names = outlier_sample_names, 
       res_std_err_values = res_std_err_values)
}



# Function to get an overview over the dataset by calculating CV and plotting it
get_CV_overview_of_dataset <- function(data, group = 'Donor', df_p = NULL, by_sample = FALSE) {
  # Helper function to calculate median, standard deviation, and CV
  calculate_stats <- function(df) {
    apply(df, 1, function(row) {
      median_val <- median(row, na.rm = TRUE)
      sd_val <- sd(row, na.rm = TRUE)
      cv_val <- ifelse(median_val != 0, sd_val / median_val, NA)
      c(median = median_val, sd = sd_val, cv = cv_val)
    })
  }
  
  # Transform based on the boolean flag
  if (by_sample) {
    if (is.null(df_p)) {
      stop("Metadata (df_p) must be provided when calculating stats for samples.")
    }
    
    # Calculate statistics for samples (transpose data)
    stats <- calculate_stats(t(data))
    
    # Create a summary dataframe
    stats_df <- data.frame(t(stats))
    stats_df$Sample <- rownames(stats_df)
    
    # Merge with metadata
    final_df <- merge(stats_df, df_p, by.x = "Sample", by.y = "row.names")
    
    # Plot the CV distribution with Donor (or other metadata column) on the x-axis
    ggplot(final_df, aes(x = factor(!!sym(group)), y = cv, fill = factor(!!sym(group)))) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # Suppress outliers in boxplot for clarity
      geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +  # Add jittered points
      labs(title = paste("CV Distribution by", group), x = group, y = "Coefficient of Variation (CV)", fill = group) +
      theme_minimal()
    ggsave(paste0(resultsdir, 'CV-overview-over-samples-by-', group,'.png'))
  } else {
    # Calculate statistics for features
    stats <- calculate_stats(data)
    
    # Create a summary dataframe
    stats_df <- data.frame(t(stats))
    stats_df$Feature <- rownames(stats_df)
    
    # Plot the CV distribution for features
    ggplot(stats_df, aes(x = cv)) +
      geom_histogram(binwidth = 0.01, fill = "blue", alpha = 0.7) +
      labs(title = "CV Distribution of Features", x = "Coefficient of Variation (CV)", y = "Frequency", fill = group) +
      theme_minimal()
    ggsave(paste0(resultsdir, 'CV-overview-over-samples-by-feature.png'))
  }
  
  return(stats_df)
}



# Function to improve the removeFeatures function from the qmtools package
# adapted to change min_group number for kept features
  
  # Function: Remove features based on different filtering criteria
  # Input: SummarizedExperiment or matrix-like x, method (missing, blankratio, rsd, icc), additional parameters
  # Output: Indices of features to retain
  removeFeatures_advanced <- function(x, i = NULL, method = c("missing", "blankratio", "rsd", "icc"), 
                              group, levels = NULL, cut = 0.7, min_groups = 1, 
                              blank_samples = NULL, qc_samples = NULL, bio_samples = NULL, blank_min_n = NULL, type = c("median", "mean")) {
    #method <- match.arg(method)
    
    if (!is.null(i) && inherits(x, "SummarizedExperiment")) {
      m <- assay(x, i)
      m <- as.matrix(m)
    } else if (!is.matrix(x)) {
      m <- as.matrix(x)
    }
    
    # Function: Remove features based on missing values
    # Input: matrix-like x, group vector, levels to consider, cutoff, min_groups
    # Output: Indices of features to retain
    .removeMiss <- function(m, group, levels = NULL, cut = 0.7, min_groups = 1) {
      levels <- levels %||% unique(group)
      stopifnot(all(levels %in% group))
      keep_count <- integer(nrow(m))  # Initialize a count vector for all features

      for (level in levels) {
        m_ <- m[, group == level, drop = FALSE]
        non_missing_frac <- rowSums(!is.na(m_)) / ncol(m_)
        keep_idx <- which(non_missing_frac >= cut)
        keep_count[keep_idx] <- keep_count[keep_idx] + 1  # Increment the count for features meeting the criterion
      }

      # Retain features that meet the cutoff in at least `min_groups` groups
      final_idx <- which(keep_count >= min_groups)
      sort(final_idx)
    }
    
    # Function: Remove features based on RSD with min_groups consideration
    # Input: matrix-like x, QC samples, RSD cutoff, min_groups
    # Output: Indices of features to retain
    .removeRSD <- function(m, qc_samples, cut = 0.3, min_groups = 1) {
      stopifnot(!missing(qc_samples))
      rsd <- apply(m[, qc_samples, drop = FALSE], 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
      keep_idx <- rsd <= cut
      which(rowSums(matrix(keep_idx, nrow = nrow(m))) >= min_groups)
    }
    
    # Function: Remove features based on ICC with min_groups consideration
    # Input: matrix-like x, QC samples, bio samples, ICC cutoff, min_groups
    # Output: Indices of features to retain
    .removeICC <- function(m, qc_samples, bio_samples, cut = 0.4, min_groups = 1) {
      .verify_package("nlme")
      stopifnot(!missing(qc_samples), !missing(bio_samples))
      m_ <- m[, c(bio_samples, qc_samples)]
      f <- factor(c(seq_len(length(bio_samples)), rep("Q", length(qc_samples))))
      icc_vec <- vapply(seq_len(nrow(m)), function(i) {
        d <- data.frame(y = m[i, c(bio_samples, qc_samples)], f = f)
        fit <- tryCatch(nlme::lme(y ~ 1, random = ~ 1 | f, data = d, na.action = na.omit), error = function(e) NA_real_)
        vv <- tryCatch(as.numeric(nlme::VarCorr(fit)[, "Variance"]), error = function(e) NA_real_)
        vv[1] / sum(vv)
      }, numeric(1))
      keep_idx <- icc_vec >= cut
      which(rowSums(matrix(keep_idx, nrow = nrow(m))) >= min_groups)
    }
    
    # Function: Remove features based on QC/blank ratio with min_groups consideration
    # Input: matrix-like x, blank samples, QC samples, cutoff, type, min blanks, min_groups
    # Output: Indices of features to retain
    .removeBlankRatio <- function(m, blank_samples, qc_samples, cut = 2, type = c("median", "mean"), blank_min_n = 3, min_groups = 1) {
      type <- match.arg(type)
      qc_summary <- apply(m[, qc_samples, drop = FALSE], 1, match.fun(type), na.rm = TRUE)
      blank_summary <- apply(m[, blank_samples, drop = FALSE], 1, match.fun(type), na.rm = TRUE)
      ratio <- qc_summary / blank_summary
      keep_idx <- ratio >= cut
      if (blank_min_n > 1) {
        blank_non_missing <- apply(m[, blank_samples, drop = FALSE], 1, function(x) sum(!is.na(x)))
        keep_idx <- keep_idx & (blank_non_missing >= blank_min_n)
      }
      which(rowSums(matrix(keep_idx, nrow = nrow(x))) >= min_groups)
    }
    
    idx_to_keep <- switch(method,
                          "missing" = .removeMiss(m, group, levels, cut, min_groups),
                          "blankratio" = .removeBlankRatio(m, qc_samples, cut, min_groups),
                          "rsd" = .removeRSD(m, qc_samples, bio_samples, cut, min_groups),
                          "icc" = .removeICC(m, blank_samples, qc_samples, cut, type, blank_min_n, min_groups)
                          )
    
    x[idx_to_keep, ]
  }



# Function to make calculations for plotting (normalize, log2, mean, std_dev, rank) 
# to prepare data for the plot and give them back in a list
prepare_data <- function(raw_sign_abund_features, distinguish_sample_names, normalization_list = NULL) {
  
  # multiply Chrom columns by 2, as they have a lower concentration than Kim (chrom 120uL extraction, Kim 60uL extraction)
  #chrom_columns <- grep("Chrom", names(raw_sign_abund_features))
  #raw_sign_abund_features[, chrom_columns] <- raw_sign_abund_features[, chrom_columns] * 2
  
  # Convert log2 of raw data and get stats for raw and raw log2 data
  raw_sign_abund_features_log <- raw_sign_abund_features
  raw_sign_abund_features <- calculate_stats(raw_sign_abund_features) # get stats only after saving to _log df, to dont copy columns beforehand, otherwise more columns in calculate_stats, which would then also include the mean for the log2 stats
  data_length <- ncol(raw_sign_abund_features_log)
  raw_sign_abund_features_log[, 10:data_length] <- apply(raw_sign_abund_features_log[, 10:data_length], 2, log2)
  raw_sign_abund_features_log <- calculate_stats(raw_sign_abund_features_log)
  
  # Convert log2 of normalized data and get stats for normalized and normalized log2 data
  norm_sign_abund_features <- normalize_data(raw_sign_abund_features, normalization_list)
  norm_sign_abund_features_log <- norm_sign_abund_features
  norm_sign_abund_features <- calculate_stats(norm_sign_abund_features) # get stats only after saving to _log df, to dont copy columns beforehand, otherwise more columns in calculate_stats, which would then also include the mean for the log2 stats
  data_length <- ncol(norm_sign_abund_features_log)
  norm_sign_abund_features_log[, 10:data_length] <- apply(norm_sign_abund_features_log[, 10:data_length], 2, log2)
  norm_sign_abund_features_log <- calculate_stats(norm_sign_abund_features_log)
  
  # Separate data for each sample, make log2, get stats
  raw_sample_data <- list()
  raw_sample_data_log <- list()
  for (sample_name in distinguish_sample_names) {
    info_columns <- raw_sign_abund_features[, 1:9]
    sub_data <- raw_sign_abund_features[, grep(sample_name, colnames(raw_sign_abund_features))]
    raw_sample_data[[sample_name]] <- cbind(info_columns, sub_data)
    raw_sample_data[[sample_name]] <- calculate_stats(raw_sample_data[[sample_name]])
    raw_sample_data_log[[sample_name]] <- cbind(info_columns, log2(sub_data))
    raw_sample_data_log[[sample_name]] <- calculate_stats(raw_sample_data_log[[sample_name]])
  }
  
  # Normalize data, and make log2 of it
  norm_sample_data <- list()
  norm_sample_data_log <- list()
  for (sample_name in distinguish_sample_names) {
    info_columns <- norm_sign_abund_features[, 1:9]
    sub_data <- norm_sign_abund_features[, grep(sample_name, colnames(norm_sign_abund_features))]
    norm_sample_data[[sample_name]] <- cbind(info_columns, sub_data)
    norm_sample_data[[sample_name]] <- calculate_stats(norm_sample_data[[sample_name]])
    norm_sample_data_log[[sample_name]] <- cbind(info_columns, log2(sub_data))
    norm_sample_data_log[[sample_name]] <- calculate_stats(norm_sample_data_log[[sample_name]])
  }
  
  return(list(raw_sample_data = raw_sample_data,
              raw_sample_data_log = raw_sample_data_log,
              norm_sample_data = norm_sample_data,
              norm_sample_data_log = norm_sample_data_log,
              raw_sign_abund_features = raw_sign_abund_features,
              raw_sign_abund_features_log = raw_sign_abund_features_log,
              norm_sign_abund_features = norm_sign_abund_features,
              norm_sign_abund_features_log = norm_sign_abund_features_log))
}



# a function to get lists of dataframes and combine and merge them for plotting
combine_and_merge <- function(raw_list, norm_list, sample_names) {
  combined_list <- list()
  for (i in seq_along(sample_names)) {
    # Combine columns for raw data
    raw_combined <- raw_list[[i]][c("id", "mean", "std_dev", "mean_rank")]
    # Rename columns for raw data
    raw_suffix <- paste0("_", sample_names[i])
    names(raw_combined)[2:4] <- paste0(names(raw_combined)[2:4], raw_suffix)
    
    # Combine columns for normalized data
    norm_combined <- norm_list[[i]][c("id", "mean", "std_dev", "mean_rank")]
    # Rename columns for normalized data
    norm_suffix <- paste0("_norm_", sample_names[i])
    names(norm_combined)[2:4] <- paste0(names(norm_combined)[2:4], norm_suffix)
    
    # Merge based on ID
    if (i == 1) {
      raw_combined <- cbind(raw_combined, Annotation = raw_list[[i]]$Annotation)
      combined_list <- merge(raw_combined, norm_combined, by = "id", all = TRUE)
    } else {
      combined_list <- merge(combined_list, raw_combined, by = "id", all = TRUE)
      combined_list <- merge(combined_list, norm_combined, by = "id", all = TRUE)
    }
  }
  return(combined_list)
}




# Function to adjust size of plot
fig <- function(width, heigth){
  options(repr.plot.width = width, repr.plot.height = heigth)
}



# Function to save plot as SVG
save_plot_as_svg <- function(plot, figuredir, experiment_name, plot_name) {
  # Create the directory if it doesn't exist
  if (!dir.exists(figuredir)) {
    dir.create(figuredir, recursive = TRUE)
  }
  
  # Define the file path
  filename <- file.path(figuredir, paste0(experiment_name, "_", plot_name, ".svg"))
  
  # Save the plot as SVG
  #ggsave(filename, plot = plot) # save with ggplot2, svglite needed, systemfonts macht probleme
  save_plot(filename, fig = plot, width = 30, height = 20) # save with sjPlot
}



# Function to save a plot with lots of figures in it (grid.arrange) as SVG
save_grid_plot_as_svg <- function(list_of_plots, figuredir, experiment_name, plot_name) {
  # Create the directory if it doesn't exist
  if (!dir.exists(figuredir)) {
    dir.create(figuredir, recursive = TRUE)
  }
  
  # Define the file path
  filename <- file.path(figuredir, paste0(experiment_name, "_", plot_name, ".svg"))
  
  # Save the plot as SVG
  svg(file = filename, width = 150, height = 80)
  do.call(grid.arrange, list_of_plots) 
  dev.off()
}


# Function to save a scatter plot as SVG
save_scatter_plot_as_svg <- function(plot, figuredir, experiment_name, plot_name) {
  # Create the directory if it doesn't exist
  if (!dir.exists(figuredir)) {
    dir.create(figuredir, recursive = TRUE)
  }
  
  # Define the file path
  filename <- file.path(figuredir, paste0(experiment_name, "_", plot_name, ".svg"))
  
  # Save the plot as SVG
  #ggsave(filename, plot = plot) # save with ggplot2, svglite needed, systemfonts macht probleme
  save_plot(filename, fig = plot, width = 20, height = 20) # save with sjPlot
}





#######
# all the big functions for plotting various data
#######


# Function to plot the ranked log2 of each sample and overlay the power function, marking outliers
plot_log2_with_power_overlay <- function(original_df, outlier_samples, cutoff_indices, function_values_list, 
                                         res_std_err_values, resultsdir, figuredir, exp) {
  # Function to generate the plot for each sample
  generate_sample_plot <- function(original_data, sample_name, cutoff_indices, function_values, res_std_err, outlier) {
    # Replace NA values with 1
    original_data[is.na(original_data)] <- 1
    
    # Log2 of data
    original_log2_data <- log2(original_data)
    
    # Create data frame for plotting
    df_plot <- data.frame(rank = rank(-original_log2_data, na.last = NA,  ties.method = 'first'), 
                          log2_value = original_log2_data)
    
    # Create a data frame for the function values
    function_df <- data.frame(rank = 1:length(function_values), function_values = function_values)
    
    # Plot
    p <- ggplot(df_plot, aes(x = rank, y = log2_value)) +
      geom_point(aes(color = "Experimental Data")) +  # Add color mapping for points
      geom_vline(aes(xintercept = cutoff_indices, color = "Feature Cutoff"), linetype = "dashed") +  # Add linetype mapping for cutoff lines
      geom_line(data = function_df, aes(x = rank, y = function_values, color = "Fitted Power Function", linetype = "Fitted Power Function"), linetype = "dotted") +  # Add color and linetype mapping for fitted line
      annotate("text", x = Inf, y = Inf, label = paste("Res. Std. Err.:", round(res_std_err, 3)), 
               hjust = 1.1, vjust = 1.1, size = 5) +
      labs(x = "Rank", y = "log2(AUC)", 
           title = paste("Power Curve fitting to detected features of ", sample_name, " (log2 Abundance)")) +
      theme(legend.position = "top",
            axis.text = element_text(size = 12),  # Adjust size of axis labels
            axis.title = element_text(size = 14),  # Adjust size of axis titles
            plot.title = element_text(size = 16),  # Adjust size of plot title
            legend.title = element_text(size = 12),  # Change legend title font size
            legend.text = element_text(size = 11)) +  # Change legend text font size
      ylim(0, 1.2 * max(df_plot$log2_value, na.rm = TRUE)) +  # Adjust lower limit to 0, upper limit to the calculated value
      scale_color_manual(values = c("Experimental Data" = "black", "Fitted Power Function" = "blue", "Feature Cutoff" = "red"), name = "Legend")  # Define colors and legend title

    # Add "IS OUTLIER" text if outlier is TRUE
    if (outlier) {
      p <- p + annotate("text", x = Inf, y = Inf, label = "IS OUTLIER", color = "red", size = 16, hjust = 1.1, vjust = 1.5) +
               theme(panel.background = element_rect(fill = rgb(1, 0.75, 0.8, alpha = 0.7), colour = NA))
    }
    
    return(p)
  }
  
  # Generate plots for each sample
  plots <- lapply(seq_along(original_df), function(i) {
    generate_sample_plot(original_df[[i]], colnames(original_df[i]), cutoff_indices[i], function_values_list[[i]],
                         res_std_err_values[i], i %in% outlier_samples)
  })
  
  # Save the plots to a PDF file
  pdf_file <- file.path(resultsdir, "check_quality_of_each_sample_with_powerfunction.pdf")
  pdf(file = pdf_file, width = 150, height = 80)
  grid.arrange(grobs = plots)
  dev.off()
  
  plot_name <- c("all_samples_individual_raw_check_cutoffs_and_outliers")
  
  # Save plot as SVG
  #save_grid_plot_as_svg(plots, figuredir, exp, plot_name) #just takes really long, only use when needed
}




# Function to generate mean plots to show the features over all different samples
get_all_mean_plots <- function(raw_sign_abund_features, distinguish_sample_names, resultsdir, figuredir, normalization_list = NULL) {
  
  # First have all functions to use later
  
  # Function to plot mean and standard deviation for one sample
  plot_mean_std <- function(df, title, names_to_label = NULL) {
    top_20 <- head(df[order(-df$mean), ], 20)
    num_features <- nrow(df)
    sum_mean_top20 <- sum(top_20$mean)
    sum_mean_rest <- sum(df$mean) - sum_mean_top20
    
    p <- ggplot(df, aes(x = mean_rank, y = mean), fig(10,4)) +
            geom_ribbon(aes(ymin = mean - std_dev, ymax = mean + std_dev, fill = "Standard Deviation"), alpha = 0.4) +
            geom_point(aes(color = "Mean"), size = 1) +
            geom_text_repel(data = top_20, aes(label = Annotation), color = "black", max.overlaps = Inf, 
                            min.segment.length = 0.01, box.padding = 0.5, nudge_x = 4) +
            geom_text_repel(data = df[df$Annotation %in% names_to_label, ], aes(label = Annotation), color = "darkgrey", 
                            max.overlaps = Inf, vjust = 5, min.segment.length = 0.01)
            if (grepl("norm",title)) {
              if (grepl("log",title)) {
                p <- p + labs(x = "Mean Rank", y = "log2(nAUC)", color = NULL, fill = NULL, title = title)
              } else {
                p <- p + labs(x = "Mean Rank", y = "nAUC", color = NULL, fill = NULL, title = title)
              }
            } else {
              if (grepl("log",title)) {
                p <- p + labs(x = "Mean Rank", y = "log2(AUC)", color = NULL, fill = NULL, title = title)
              } else {
                p <- p + labs(x = "Mean Rank", y = "AUC", color = NULL, fill = NULL, title = title)
              }
            }
    p <- p + scale_color_manual(values = c("Mean" = "darkblue"), guide = guide_legend(title = "Mean")) +
             scale_fill_manual(values = "skyblue", guide = guide_legend(title = "Standard Deviation")) +
             theme(legend.position = "top",
                   axis.text = element_text(size = 12),  # Adjust size of axis labels
                   axis.title = element_text(size = 14),  # Adjust size of axis titles
                   plot.title = element_text(size = 16),  # Adjust size of plot title
                   legend.title = element_text(size=12), #change legend title font size
                   legend.text = element_text(size=11)) +  #change legend text font size
          if (grepl("log",title)) {
             annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste("Number of Features:", num_features))  # Add text with number of features
          } else {
             annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste("Number of Features:", num_features, "\n",
                                                                                   "Sum of mean_norm (Top 20):", round(sum_mean_top20,2), "\n",
                                                                                   "Sum of mean_norm (Rest):", round(sum_mean_rest,2)))  # Add text with number of features and sum of mean_norm to top right corner
          }
    
    return(p)
  }
  
  
  # Function to merge data with sample_name as factor for plotting
  merge_dataframes <- function(df_list, distinguish_sample_names) {
    # Create an empty list to store modified data frames
    merged_df_list <- list()
    
    # Iterate through each sample name and corresponding data frame
    for (sample_name in distinguish_sample_names) {
      # Retrieve the data frame for the current sample name
      df <- df_list[[sample_name]]
      
      # Add only the needed columns (mean_rank, mean, std_dev)
      df <- df[, c("id", "rt", "mz", "rtime_group", "feature_group", "Annotation", "Confidence", "SMILES", "mean", "std_dev", "rel_std_dev", "mean_rank")]
      
      # Add a group column with the sample name
      df$group <- as.factor(sample_name)
      
      # Append the modified data frame to the list
      merged_df_list[[sample_name]] <- df
    }
    
    # Combine all data frames into one
    merged_df <- do.call(rbind, merged_df_list)
    
    return(merged_df)
  }
  
  
  
  
  
  # Function to plot mean and standard deviation from several samples in one figure
  plot_mean_std_combine <- function(df_list, distinguish_sample_names, title, names_to_label = NULL) {
    num_features <- nrow(df_list[[1]])  # Assuming all dataframes in df_list have the same number of features
    
    # Merge data frames
    merged_df <- merge_dataframes(df_list, distinguish_sample_names)
    
    p <- ggplot(merged_df, aes(x = mean_rank, y = mean), fig(10,4)) +
      geom_ribbon(aes(x = mean_rank, ymin = mean - std_dev, ymax = mean + std_dev, fill = group), alpha = 0.3) +  
      geom_point(aes(x = mean_rank, y = mean, color = group), size = 1)
      
    if (grepl("norm",title)) {
      if (grepl("log",title)) {
        p <- p + labs(x = "Mean Rank", y = "log2(nAUC)", color = NULL, fill = NULL, title = title)
      } else {
        p <- p + labs(x = "Mean Rank", y = "nAUC", color = NULL, fill = NULL, title = title)
      }
    } else {
      if (grepl("log",title)) {
        p <- p + labs(x = "Mean Rank", y = "log2(AUC)", color = NULL, fill = NULL, title = title)
      } else {
        p <- p + labs(x = "Mean Rank", y = "AUC", color = NULL, fill = NULL, title = title)
      }
    }
    
    p <- p + theme(legend.position = "top") +
      annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste("Number of Features:", num_features)) +
      scale_color_manual(name = "Mean", values = turbo(length(distinguish_sample_names)), labels = distinguish_sample_names) +  # Set manual color scale
      scale_fill_manual(name = "Standard Deviation", values = turbo(length(distinguish_sample_names)), labels = distinguish_sample_names) +  # Set manual fill scale
      theme(legend.position = "top",
            axis.text = element_text(size = 12),  # Adjust size of axis labels
            axis.title = element_text(size = 14),  # Adjust size of axis titles
            plot.title = element_text(size = 16),  # Adjust size of plot title
            legend.title = element_text(size=12), #change legend title font size
            legend.text = element_text(size=11))  #change legend text font size
    
    return(p)
  }
  
  # start of data preparation and then plotting
  
  # Prepare data
  prepared_data <- prepare_data(raw_sign_abund_features, distinguish_sample_names, normalization_list)
  
  raw_sample_data <- prepared_data$raw_sample_data
  raw_sample_data_log <- prepared_data$raw_sample_data_log
  norm_sample_data <- prepared_data$norm_sample_data
  norm_sample_data_log <- prepared_data$norm_sample_data_log
  raw_sign_abund_features <- prepared_data$raw_sign_abund_features
  raw_sign_abund_features_log <- prepared_data$raw_sign_abund_features_log
  norm_sign_abund_features <- prepared_data$norm_sign_abund_features
  norm_sign_abund_features_log <- prepared_data$norm_sign_abund_features_log
  
  
  # Generate plots
  pdf(file = file.path(resultsdir, "feature_ranking_sample-mean.pdf"), width = 15, height = 15)
  
  # Plot absolute values and log2 values for all raw data
  all_data_plot <- plot_mean_std(raw_sign_abund_features, "Feature Ranking of raw Abundance (All Data)")
  all_data_log_plot <- plot_mean_std(raw_sign_abund_features_log, "Feature Ranking of raw log2 Abundance (All Data)")
  
  plot_list_raw <- list(all_data_plot, all_data_log_plot)
  sample_plot_names <- c("all_data_mean-rank_raw_plot", "all_data_mean-rank_raw_log_plot")
  
  # Plot absolute raw values and log2 values for each sample
  for (sample_name in distinguish_sample_names) {
    sample_data_plot <- plot_mean_std(raw_sample_data[[sample_name]], paste0("Feature Ranking of raw Abundance (", sample_name, ")"))
    sample_data_log_plot <- plot_mean_std(raw_sample_data_log[[sample_name]], paste0("Feature Ranking of raw log2 Abundance (", sample_name, ")"))
    
    # Assign names to plots
    sample_plot_names <- c(sample_plot_names, paste0(sample_name, "mean-rank_raw_plot"),
                           paste0(sample_name, "mean-rank_raw_log_plot"))
    
    plot_list_raw <- c(plot_list_raw, list(sample_data_plot, sample_data_log_plot))
  }

    # Arrange and display plots
  do.call(grid.arrange, plot_list_raw)
  
  
  # Plot absolute values and log2 values for all normalized data
  all_data_plot <- plot_mean_std(norm_sign_abund_features, "Feature Ranking of normalized Abundance (All Data)")
  all_data_log_plot <- plot_mean_std(norm_sign_abund_features_log, "Feature Ranking of normalized log2 Abundance (All Data)")
  
  plot_list_norm <- list(all_data_plot, all_data_log_plot)
  sample_plot_names <- c(sample_plot_names, "all_data_mean-rank_norm_plot", "all_data_mean-rank_norm_log_plot")
  
  # Plot absolute normalized values and log2 values for each sample
  for (sample_name in distinguish_sample_names) {
    sample_data_plot <- plot_mean_std(norm_sample_data[[sample_name]], paste("Feature Ranking of normalized Abundance (", sample_name, ")", sep = ""))
    sample_data_log_plot <- plot_mean_std(norm_sample_data_log[[sample_name]], paste("Feature Ranking of normalized log2 Abundance (", sample_name, ")", sep = ""))
    
    # Assign names to plots
    sample_plot_names <- c(sample_plot_names, paste0(sample_name, "mean-rank_norm_plot"),
                           paste0(sample_name, "mean-rank_norm_log_plot"))
    
    plot_list_norm <- c(plot_list_norm, list(sample_data_plot, sample_data_log_plot))
  }
  
  # Arrange and display plots
  do.call(grid.arrange, plot_list_norm)
  

  # plot all raw sample types in one plot
  all_data_plot <- plot_mean_std_combine(raw_sample_data, distinguish_sample_names, "Compare Feature Ranking of raw Abundance between various Sample Types")
  all_data_log_plot <- plot_mean_std_combine(raw_sample_data_log, distinguish_sample_names, "Compare Feature Ranking of raw log2 Abundance between various Sample Types")
  
  plot_list_combine <- list(all_data_plot, all_data_log_plot)
  sample_plot_names <- c(sample_plot_names, "samples-combined_mean-rank_raw_plot", "samples-combined_mean-rank_raw_log_plot")
  
  # Arrange and display plots
  do.call(grid.arrange, plot_list_combine)
  
   
  # plot all normalized sample types in one plot
  all_data_plot <- plot_mean_std_combine(norm_sample_data, distinguish_sample_names, "Compare Feature Ranking of normalized Abundance between various Sample Types")
  all_data_log_plot <- plot_mean_std_combine(norm_sample_data_log, distinguish_sample_names, "Compare Feature Ranking of normalized log2 Abundance between various Sample Types")
  
  plot_list_combine_norm <- list(all_data_plot, all_data_log_plot)
  sample_plot_names <- c(sample_plot_names, "samples-combined_mean-rank_norm_plot", "samples-combined_mean-rank_norm_log_plot")
  
  # Arrange and display plots
  do.call(grid.arrange, plot_list_combine_norm)
  
  
  dev.off() # Close PDF device
  
  
  # save plots as svg in directory
    # List of plots
    plot_list <- c(
      plot_list_raw,
      plot_list_norm,
      plot_list_combine,
      plot_list_combine_norm
    )
    
    # Save each plot as SVG
    for (i in seq_along(plot_list)) {
      save_plot_as_svg(plot_list[[i]], figuredir, exp, sample_plot_names[i])
    }
}






# Function: Plot the ranking of features by abundance for each sample and save to a PDF file
# Input: 
#   - significant_abundant_features: Data frame containing significant and abundant features with samples as columns and features as rows.
#   - resultsdir: Directory where the PDF file will be saved.
# Output:
#   - PDF file containing plots of feature rankings for each sample.
plot_feature_ranking_each_sample_to_pdf <- function(significant_abundant_features, resultsdir, figuredir, normalization_list = NULL) {
  # Open a PDF device for plotting
  pdf_file <- file.path(resultsdir, "feature_ranking_each-sample.pdf")
  pdf(file = pdf_file, width = 150, height = 80)
  
  # plot raw abundance value
  plot_list_raw <- list()
  data_length <- ncol(significant_abundant_features)
  
  # Loop over columns 9 to data_length (all samples from study)
  for (col_index in seq(10, data_length)) {
    df_rank <- data.frame(name = significant_abundant_features$Annotation, sample_area = significant_abundant_features[, col_index])
    df_rank <- na.omit(df_rank)
    
    # Rank the rows based on the values in the current column
    df_rank$rank <- rank(-df_rank$sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
    
    # Identify the top 20 rows
    top20 <- head(df_rank[order(-df_rank$sample_area), ], 20)
    
    # Create the plot for the current column
    p <- ggplot(df_rank, aes(x = rank, y = sample_area)) +
      geom_point() +
      geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
      labs(x = "Rank", y = "AUC",
           title = paste("Feature Ranking of", colnames(significant_abundant_features)[col_index], "(Raw Abundance)"))
    
    # Add the plot to the list
    plot_list_raw[[col_index - 9]] <- p
  }
  
  # Arrange the raw abundance plots in a grid
  do.call(grid.arrange, plot_list_raw)
  
  # plot log2 abundance value
  plot_list_raw_log <- list()
  
  # make raw abundance to log abundance
  significant_abundant_features_log <- significant_abundant_features[, 1:data_length]
  significant_abundant_features_log[, 10:data_length] <- apply(significant_abundant_features_log[, 10:data_length], 2, log2)
  
  # Loop over columns 9 to data_length (all samples from study)
  for (col_index in seq(10, data_length)) {
    df_rank <- data.frame(name = significant_abundant_features_log$Annotation, sample_area = significant_abundant_features_log[, col_index])
    df_rank <- na.omit(df_rank)
    
    # Rank the rows based on the values in the current column
    df_rank$rank <- rank(-df_rank$sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
    
    # Identify the top 20 rows
    top20 <- head(df_rank[order(-df_rank$sample_area), ], 20)
    
    # Create the plot for the current column
    p <- ggplot(df_rank, aes(x = rank, y = sample_area)) +
      geom_point() +
      geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
      labs(x = "Rank", y = "log(AUC)",
           title = paste("Feature Ranking of", colnames(significant_abundant_features)[col_index], "(log2 Abundance)"))
    
    # Add the plot to the list
    plot_list_raw_log[[col_index - 9]] <- p
  }
  
  # Arrange the log2 abundance plots in a grid
  do.call(grid.arrange, plot_list_raw_log)
  
  
  # Plot normalized abundance value
  plot_list_norm <- list()
  
  # Normalize data
  significant_abundant_features_norm <- normalize_data(significant_abundant_features, normalization_list)
  
  # Loop over columns 9 to data_length (all samples from study)
  for (col_index in seq(10, data_length)) {
    df_rank <- data.frame(name = significant_abundant_features_norm$Annotation, sample_area = significant_abundant_features_norm[, col_index])
    df_rank <- na.omit(df_rank)
    
    # Rank the rows based on the values in the current column
    df_rank$rank <- rank(-df_rank$sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
    
    # Identify the top 20 rows
    top20 <- head(df_rank[order(-df_rank$sample_area), ], 20)
    
    # Create the plot for the current column
    p <- ggplot(df_rank, aes(x = rank, y = sample_area)) +
      geom_point() +
      geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
      labs(x = "Rank", y = "nAUC",
           title = paste("Feature Ranking of", colnames(significant_abundant_features)[col_index], "(normalized Abundance)"))
    
    # Add the plot to the list
    plot_list_norm[[col_index - 9]] <- p
  }
  
  # Arrange the normalized abundance plots in a grid
  do.call(grid.arrange, plot_list_norm)
  
  # Plot normalized log2 abundance value
  plot_list_norm_log <- list()
  
  # Loop over columns 9 to data_length (all samples from study)
  for (col_index in seq(10, data_length)) {
    df_rank <- data.frame(name = significant_abundant_features_norm$Annotation, sample_area = significant_abundant_features_norm[, col_index])
    df_rank <- na.omit(df_rank)
    
    # Convert normalized abundance to log2 abundance
    df_rank$log2_sample_area <- log2(df_rank$sample_area)
    
    # Rank the rows based on the values in the current column
    df_rank$rank <- rank(-df_rank$log2_sample_area, na.last = NA,  ties.method = 'first')  # Use -df_rank[, col_index] for descending order
    
    # Identify the top 20 rows
    top20 <- head(df_rank[order(-df_rank$log2_sample_area), ], 20)
    
    # Create the plot for the current column
    p <- ggplot(df_rank, aes(x = rank, y = log2_sample_area)) +
      geom_point() +
      geom_text_repel(data = top20, aes(label = name, color = "Top20"), max.overlaps = Inf, min.segment.length = 0.01) + 
      labs(x = "Rank", y = "log2(nAUC)",
           title = paste("Feature Ranking of", colnames(significant_abundant_features)[col_index], "(normalized log2 Abundance)"))
    
    # Add the plot to the list
    plot_list_norm_log[[col_index - 9]] <- p
  }
  
  # Arrange the normalized log2 abundance plots in a grid
  do.call(grid.arrange, plot_list_norm_log)
  
  
  # Close the PDF device
  dev.off()
  
  
  # save plots as svg in directory
  # List of plots
  plot_list <- list(
    plot_list_raw,
    plot_list_raw_log,
    plot_list_norm,
    plot_list_norm_log
  )
  
  # List of plot names
  plot_names <- c("all_samples_individual_raw", "all_samples_individual_raw_log", "all_samples_individual_norm", "all_samples_individual_norm_log")
  
  # Save each plot as SVG
  for (i in seq_along(plot_list)) {
    save_grid_plot_as_svg(plot_list[[i]], figuredir, exp, plot_names[i])
  }
  
}



# Function: Plot the ranking of features by abundance for each sample and save to a PDF file
# Input: 
#   - all_needed_features: Data frame containing significant and abundant features with samples as columns and features as rows.
#   - resultsdir: Directory where the PDF file will be saved.
# Output:
#   - PDF file containing plots of feature rankings for each sample.
plot_feature_ranking_overview_to_pdf <- function(all_needed_features, resultsdir) {
  # Open a PDF device for plotting
  pdf_file <- file.path(resultsdir, "feature_ranking_data-overview.pdf")
  pdf(file = pdf_file, width = 15, height = 8)
  
  # plot mean area to rank
  data_length <- ncol(all_needed_features)
  
  # Compute rank for each sample
  rank_long <- all_needed_features %>%
    pivot_longer(cols = all_of(10:data_length), names_to = "Sample", values_to = "Sample_Area") %>%
    group_by(Sample) %>%
    mutate(Rank = rank(-Sample_Area, ties.method = "first")) %>%
    ungroup()
  
  # Compute mean and standard deviation of ranks
  rank_summary <- rank_long %>%
    group_by(id) %>%
    summarise(Mean_Rank = mean(Rank, na.rm = TRUE),
              SD_Rank = sd(Rank, na.rm = TRUE))
  
  # Compute mean and standard deviation of areas
  area_summary <- all_needed_features %>%
    rowwise() %>%
    mutate(Mean_Area = mean(c_across(all_of(10:data_length)), na.rm = TRUE),
           SD_Area = sd(c_across(all_of(10:data_length)), na.rm = TRUE), 
           log2_Mean_Area = mean(log2(c_across(all_of(10:data_length))), na.rm = TRUE),
           log2_SD_Area = sd(log2(c_across(all_of(10:data_length))), na.rm = TRUE)) %>%
    ungroup()
  
  # Merge summaries
  summary_df <- left_join(rank_summary, area_summary, by = "id")
  
  # Plot Mean Area vs. Mean Rank with SD ribbons
  p <- ggplot(summary_df, aes(x = Mean_Area, y = Mean_Rank)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Rank - SD_Rank, ymax = Mean_Rank + SD_Rank), alpha = 0.2) +
    labs(x = "Mean AUC", y = "Mean Rank",
         title = "Feature Ranking Based on Mean Abundance") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Plot Mean Area vs. log2(Mean Rank) with SD ribbons
  p <- ggplot(summary_df, aes(x = log2(Mean_Area), y = Mean_Rank)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Rank - SD_Rank, ymax = Mean_Rank + SD_Rank), alpha = 0.2) +
    labs(x = "log2(Mean AUC)", y = "Mean Rank",
         title = "Feature Ranking Based on Mean Abundance (Log Scale)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Rank features based on mean area
  area_summary <- area_summary %>%
    mutate(Mean_Rank = rank(-Mean_Area, ties.method = "first"))
  
  # Plot Mean Rank vs. Mean Area with SD ribbons
  p <- ggplot(area_summary, aes(x = Mean_Rank, y = Mean_Area)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Area - SD_Area, ymax = Mean_Area + SD_Area), alpha = 0.2) +
    labs(x = "Mean Rank", y = "Mean AUC",
         title = "Feature Ranking Based on Mean Abundance") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Plot Mean Rank vs. log2(Mean Area) with SD ribbons
  p <- ggplot(area_summary, aes(x = Mean_Rank, y = log2_Mean_Area)) +
    geom_point() +
    geom_ribbon(aes(ymin = log2_Mean_Area - log2_SD_Area, ymax = log2_Mean_Area + log2_SD_Area), alpha = 0.2) +
    labs(x = "Mean Rank", y = "log2(Mean AUC)",
         title = "Feature Ranking Based on Mean Abundance (Log Scale)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Close the PDF device
  dev.off()
}


plot_feature_ranking_overview_sample_filtering_to_pdf <- function(all_needed_features, resultsdir, add_info, Group, select) {
  # Open a PDF device for plotting
  pdf_file <- paste0(resultsdir, "feature_ranking_data-overview", Group, "-", select ,".pdf")
  pdf(file = pdf_file, width = 15, height = 8)
  
  # plot mean area to rank
  data_length <- ncol(all_needed_features)
  
  # Filter metadata based on Group and select value
  selected_samples <- add_info %>%
    filter(.data[[Group]] == select) %>%  # Dynamically select column
    pull(Sample)  # Extract relevant Sample names
  
  # Ensure selected_samples exist in the dataset
  selected_samples <- intersect(selected_samples, colnames(all_needed_features)[10:data_length])
  
  # Compute rank for each sample in the filtered set
  rank_long <- all_needed_features %>%
    select(id, all_of(selected_samples)) %>%  # Keep only selected samples
    pivot_longer(cols = all_of(selected_samples), names_to = "Sample", values_to = "Sample_Area") %>%
    group_by(Sample) %>%
    mutate(Rank = rank(-Sample_Area, ties.method = "first")) %>%
    ungroup()
  
  # Compute mean and standard deviation of ranks
  rank_summary <- rank_long %>%
    group_by(id) %>%
    summarise(Mean_Rank = mean(Rank, na.rm = TRUE),
              SD_Rank = sd(Rank, na.rm = TRUE))
  
  # Compute mean and standard deviation of areas
  area_summary <- all_needed_features %>%
    rowwise() %>%
    mutate(Mean_Area = mean(c_across(all_of(selected_samples)), na.rm = TRUE),
           SD_Area = sd(c_across(all_of(selected_samples)), na.rm = TRUE), 
           log2_Mean_Area = mean(log2(c_across(all_of(selected_samples))), na.rm = TRUE),
           log2_SD_Area = sd(log2(c_across(all_of(selected_samples))), na.rm = TRUE)) %>%
    ungroup()
  
  # Merge summaries
  summary_df <- left_join(rank_summary, area_summary, by = "id")
  
  # Plot Mean Area vs. Mean Rank with SD ribbons
  p <- ggplot(summary_df, aes(x = Mean_Area, y = Mean_Rank)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Rank - SD_Rank, ymax = Mean_Rank + SD_Rank), alpha = 0.2) +
    labs(x = "Mean AUC", y = "Mean Rank",
         title = "Feature Ranking Based on Mean Abundance") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Plot Mean Area vs. log2(Mean Rank) with SD ribbons
  p <- ggplot(summary_df, aes(x = log2(Mean_Area), y = Mean_Rank)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Rank - SD_Rank, ymax = Mean_Rank + SD_Rank), alpha = 0.2) +
    labs(x = "log2(Mean AUC)", y = "Mean Rank",
         title = "Feature Ranking Based on Mean Abundance (Log Scale)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Rank features based on mean area
  area_summary <- area_summary %>%
    mutate(Mean_Rank = rank(-Mean_Area, ties.method = "first"))
  
  # Plot Mean Rank vs. Mean Area with SD ribbons
  p <- ggplot(area_summary, aes(x = Mean_Rank, y = Mean_Area)) +
    geom_point() +
    geom_ribbon(aes(ymin = Mean_Area - SD_Area, ymax = Mean_Area + SD_Area), alpha = 0.2) +
    labs(x = "Mean Rank", y = "Mean AUC",
         title = "Feature Ranking Based on Mean Abundance") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Plot Mean Rank vs. log2(Mean Area) with SD ribbons
  p <- ggplot(area_summary, aes(x = Mean_Rank, y = log2_Mean_Area)) +
    geom_point() +
    geom_ribbon(aes(ymin = log2_Mean_Area - log2_SD_Area, ymax = log2_Mean_Area + log2_SD_Area), alpha = 0.2) +
    labs(x = "Mean Rank", y = "log2(Mean AUC)",
         title = "Feature Ranking Based on Mean Abundance (Log Scale)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14), 
      plot.title = element_text(size = 16, face = "bold") 
    )
  
  print(p)
  
  # Close the PDF device
  dev.off()
}



# Function to generate scatter plots and filter features based on conditions for any two sample names
generate_scatter_plots_and_get_filterd_featrues <- function(raw_sign_abund_features, sample_names, resultsdir, add_further_labels = FALSE,
                                   names_for_mean_comparison = NULL, names_for_std_dev_comparison = NULL, names_for_rank_comparison = NULL,
                                   normalization_list = NULL) {
  
  # Check if the sample_names list contains exactly two names
  if (length(sample_names) != 2) {
    print("Sample list needs 2 names. Please choose only 2 for comparison.")
    return()
  }
  
  # Prepare data 
  prepared_data <- prepare_data(raw_sign_abund_features, sample_names, normalization_list)
  
  #raw_sample_data <- prepared_data$raw_sample_data
  raw_sample_data_log <- prepared_data$raw_sample_data_log
  #norm_sample_data <- prepared_data$norm_sample_data
  norm_sample_data_log <- prepared_data$norm_sample_data_log
  
  
  # merge dataframes to for comparison plots
  combined_data <- combine_and_merge(raw_sample_data_log, norm_sample_data_log, sample_names)
  
  
  # Define conditions for selecting features that are different in the comparison
  condition_mean <- with(combined_data, (get(paste0("mean_", sample_names[[1]])) > 1.5 * get(paste0("mean_", sample_names[[2]])) |
                           get(paste0("mean_", sample_names[[1]])) * 1.5 < get(paste0("mean_", sample_names[[2]])) ) )#&
                           #get(paste0("mean_", sample_names[[2]])) > 13)
  
  condition_std_dev <- with(combined_data, 
                            (get(paste0("std_dev_", sample_names[[1]])) > 2 * get(paste0("std_dev_", sample_names[[2]])) |
                             get(paste0("std_dev_", sample_names[[1]])) * 2 < get(paste0("std_dev_", sample_names[[2]]))) &
                            (get(paste0("std_dev_", sample_names[[1]])) > 2.5 | 
                             get(paste0("std_dev_", sample_names[[2]])) > 2.5))
  
  condition_rank <- with(combined_data, (get(paste0("mean_rank_", sample_names[[1]])) > 1.8 * get(paste0("mean_rank_", sample_names[[2]])) |
                           get(paste0("mean_rank_", sample_names[[1]])) * 1.8 < get(paste0("mean_rank_", sample_names[[2]])) ) )#&
                           #get(paste0("mean_rank_", sample_names[[1]])) > 190 & 
                           #get(paste0("mean_rank_", sample_names[[2]])) < 350)
  
  # Filter annotations based on conditions (if labels in the plot are wanted)
  filtered_annot_mean <- combined_data$Annotation[condition_mean & !is.na(combined_data$Annotation) & combined_data$Annotation != ""]
  filtered_annot_std_dev <- combined_data$Annotation[condition_std_dev & !is.na(combined_data$Annotation) & combined_data$Annotation != ""]
  filtered_annot_rank <- combined_data$Annotation[condition_rank & !is.na(combined_data$Annotation) & combined_data$Annotation != ""]
  
  # Update names for comparison with filtered annotations
  if (add_further_labels == TRUE) {
    names_for_mean_comparison <- c(names_for_mean_comparison, filtered_annot_mean)
    names_for_std_dev_comparison <- c(names_for_std_dev_comparison, filtered_annot_std_dev)
    names_for_rank_comparison <- c(names_for_rank_comparison, filtered_annot_rank)
  }
    
  # Generate scatter plots
  pdf(file = file.path(resultsdir, "scatter_plots_comparison.pdf"), width = 18, height = 15)
  

  # Scatter plot for raw mean comparison
  plot_mean_comparison <- ggplot(combined_data, aes(x = get(paste0("mean_", sample_names[[1]])), 
                                                    y = get(paste0("mean_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_mean, ], aes(x = get(paste0("mean_", sample_names[[1]])), y = get(paste0("mean_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("log2(AUC)", sample_names[[1]]), y = paste("log2(AUC)", sample_names[[2]]), title = paste0("Comparison of Mean of raw mean Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("mean_", sample_names[[1]] )]], combined_data[[paste0("mean_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_", sample_names[[1]] )]], combined_data[[paste0("mean_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("mean_", sample_names[[1]] )]], combined_data[[paste0("mean_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_", sample_names[[1]] )]], combined_data[[paste0("mean_", sample_names[[2]] )]]))) +
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Scatter plot for norm mean comparison
  plot_mean_comparison_norm <- ggplot(combined_data, aes(x = get(paste0("mean_norm_", sample_names[[1]])), 
                                                         y = get(paste0("mean_norm_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_mean, ], aes(x = get(paste0("mean_norm_", sample_names[[1]])), y = get(paste0("mean_norm_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("log2(nAUC)", sample_names[[1]]), y = paste("log2(nAUC)", sample_names[[2]]), title = paste0("Comparison of Mean of normalized Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("mean_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_norm_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("mean_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_norm_", sample_names[[2]] )]]))) +
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Scatter plot for raw std dev comparison
  plot_std_dev_comparison <- ggplot(combined_data, aes(x = get(paste0("std_dev_", sample_names[[1]])), 
                                                       y = get(paste0("std_dev_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_std_dev, ], aes(x = get(paste0("std_dev_", sample_names[[1]])), y = get(paste0("std_dev_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("log2(AUC)", sample_names[[1]]), y = paste("log2(AUC)", sample_names[[2]]), title = paste0("Comparison of Std. dev. of raw Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("std_dev_", sample_names[[1]] )]], combined_data[[paste0("std_dev_", sample_names[[2]] )]]),
           max(combined_data[[paste0("std_dev_", sample_names[[1]] )]], combined_data[[paste0("std_dev_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("std_dev_", sample_names[[1]] )]], combined_data[[paste0("std_dev_", sample_names[[2]] )]]),
           max(combined_data[[paste0("std_dev_", sample_names[[1]] )]], combined_data[[paste0("std_dev_", sample_names[[2]] )]]))) +
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Scatter plot for norm std dev comparison
  plot_std_dev_comparison_norm <- ggplot(combined_data, aes(x = get(paste0("std_dev_norm_", sample_names[[1]])), 
                                                            y = get(paste0("std_dev_norm_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_std_dev, ], aes(x = get(paste0("std_dev_norm_", sample_names[[1]])), y = get(paste0("std_dev_norm_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("log2(nAUC)", sample_names[[1]]), y = paste("log2(nAUC)", sample_names[[2]]), title = paste0("Comparison of Std. dev. of normalized Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("std_dev_norm_", sample_names[[1]] )]], combined_data[[paste0("std_dev_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("std_dev_norm_", sample_names[[1]] )]], combined_data[[paste0("std_dev_norm_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("std_dev_norm_", sample_names[[1]] )]], combined_data[[paste0("std_dev_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("std_dev_norm_", sample_names[[1]] )]], combined_data[[paste0("std_dev_norm_", sample_names[[2]] )]]))) +
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Scatter plot for raw rank comparison
  plot_rank_comparison <- ggplot(combined_data, aes(x = get(paste0("mean_rank_", sample_names[[1]])), 
                                                    y = get(paste0("mean_rank_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_rank, ], aes(x = get(paste0("mean_rank_", sample_names[[1]])), y = get(paste0("mean_rank_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("rank of log2(mean AUC)", sample_names[[1]]), y = paste("rank of log2(mean AUC)", sample_names[[2]]), title = paste0("Comparison of Rank of raw Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("mean_rank_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_rank_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("mean_rank_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_rank_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_", sample_names[[2]] )]]))) +
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Scatter plot for norm rank comparison
  plot_rank_comparison_norm <- ggplot(combined_data, aes(x = get(paste0("mean_rank_norm_", sample_names[[1]])), 
                                                         y = get(paste0("mean_rank_norm_", sample_names[[2]]))), fig(10,10)) +
    geom_abline(intercept = 0, slope = 1, color = "gray80", linewidth = 0.4) +  # 45? line
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue", linetype = "solid", linewidth = 0.8) +  # Regression line
    geom_point(data = combined_data[condition_rank, ], aes(x = get(paste0("mean_rank_norm_", sample_names[[1]])), y = get(paste0("mean_rank_norm_", sample_names[[2]]))), color = "red") +  # Red points for filtered data
    labs(x = paste("rank of log2(mean nAUC)", sample_names[[1]]), y = paste("rank of log2(mean nAUC)", sample_names[[2]]), title = paste0("Comparison of Rank of normalized Feature Abundance (", sample_names[[1]], " vs ", sample_names[[2]], ")")) +
    coord_fixed() +
    xlim(c(min(combined_data[[paste0("mean_rank_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_rank_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_norm_", sample_names[[2]] )]]))) +
    ylim(c(min(combined_data[[paste0("mean_rank_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_norm_", sample_names[[2]] )]]),
           max(combined_data[[paste0("mean_rank_norm_", sample_names[[1]] )]], combined_data[[paste0("mean_rank_norm_", sample_names[[2]] )]]))) + 
    theme(axis.text = element_text(size = 12),  # Adjust size of axis labels
          axis.title = element_text(size = 14),  # Adjust size of axis titles
          plot.title = element_text(size = 16))  # Adjust size of plot title
  
  # Add all plots to the plot_list
  plot_list <- list(plot_mean_comparison, plot_mean_comparison_norm, plot_std_dev_comparison, plot_std_dev_comparison_norm, plot_rank_comparison, plot_rank_comparison_norm)
  
  # Arrange and display plots
  do.call(grid.arrange, plot_list)  
  
  dev.off() # Close the PDF device
  
  
  # save plots as svg in directory
  # List of plots
  plot_list <- list(
    plot_mean_comparison,
    plot_mean_comparison_norm,
    plot_std_dev_comparison,
    plot_std_dev_comparison_norm,
    plot_rank_comparison,
    plot_rank_comparison_norm
  )
  
  # List of plot names
  plot_names <- c("scatter_mean_comparison", "scatter_mean_comparison_norm", 
                  "scatter_std_dev_comparison", "scatter_std_dev_comparison_norm", 
                  "scatter_rank_comparison", "scatter_rank_comparison_norm")
  
  # Save each plot as SVG
  for (i in seq_along(plot_list)) {
    save_scatter_plot_as_svg(plot_list[[i]], figuredir, exp, plot_names[i])
  }
  
  
  
  
  # Filter features based on conditions and return
  filtered_features <- list(
    #mean_comparison = raw_sign_abund_features[condition_mean, ],
    mean_comparison = combined_data[condition_mean, ],
    std_dev_comparison = combined_data[condition_std_dev, ],
    rank_comparison = combined_data[condition_rank, ]
  )
  
  return(filtered_features)
}



# Function to categorize values into three groups based on median and standard deviation
categorize_groups <- function(column) {
  median_val <- median(column, na.rm = TRUE)
  std_dev <- sd(column, na.rm = TRUE)
  
  if (median_val > 5) {
    lower_bound <- ceiling(median_val - std_dev)
    upper_bound <- floor(median_val + std_dev)
  } else {
    lower_bound <- round(median_val - std_dev, 2)
    upper_bound <- round(median_val + std_dev, 2)
  }
  
  factor(ifelse(column < lower_bound, paste0("<", lower_bound),
                ifelse(column <= upper_bound, paste0(lower_bound, " - ", upper_bound),
                       paste0(">", upper_bound))))
}



# Function: Plot data to get an overview of the data
# Input: resultsdir (directory for result files), all_needed_features, add_info, list of Groups
# Output: PDF with plots
plot_to_pdf_objective_data_overview <- function(resultsdir, all_needed_features, add_info, Group_list, Groups_to_group, color_donor_only = FALSE) {
  
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
  
  reshaped_data <- reshape2::dcast(reshaped_data, variable ~ id, value.var = "value")
  colnames(reshaped_data)[1] <- "Sample"
  rownames(reshaped_data) <- reshaped_data[,1]
  pure_data <- reshaped_data[,-1]
  
  # Merge reshaped data with add_info
  merged_data <- merge(reshaped_data, add_info, by.x = "Sample", by.y = "Sample")
  
  merged_data[, 2:(length(merged_data) - (ncol(add_info) -1))] <- as.data.frame(apply(merged_data[, 2:(length(merged_data) - (ncol(add_info) -1))], 2, replace_na_and_Inf))
  pure_data <- as.data.frame(apply(pure_data, 2, replace_na_and_Inf))
  rownames(pure_data) <- reshaped_data[,1]
  
  # Group the values in the column into 3 'normal distributed' groups
  for (col in Groups_to_group) {
    if (col %in% colnames(merged_data)) {
      merged_data[[col]] <- categorize_groups(merged_data[[col]])
    }
  }
  
  merged_data_all <- merged_data
  pure_data_all <- pure_data
  
  # Set up PDF device
  if (color_donor_only) {
    pdf(file = paste0(resultsdir, "data_overview_for_different_Groups_Donor-coloring.pdf"), width = 15, height = 5 * (length(Group_list) %/% 2 + 1)) 
    
  } else {
    pdf(file = paste0(resultsdir, "data_overview_for_different_Groups.pdf"), width = 15, height = 5 * (length(Group_list) %/% 2 + 1)) 
  }
  
  # List of comparison methods
  comparison_methods <- c("PCA", "PLSDA", "RLA")
  
  # Loop over each comparison method
  for (method in comparison_methods) {
    plot_list <- list()
    
    for(Group in Group_list) {
      merged_data <- merged_data_all[!is.na(merged_data_all[[Group]]), ] # remove samples that have NA for this group
      pure_data <- pure_data_all[rownames(pure_data_all) %in% merged_data$Sample,]
      merged_data[[Group]] <- as.factor(merged_data[[Group]])
      merged_data$Timepoint <- as.factor(merged_data$Timepoint)
      merged_data$Donor <- as.factor(merged_data$Donor)
      
      # Calculate sample counts
      sample_counts <- table(merged_data[[Group]])
      sample_labels <- paste0(names(sample_counts), "\n(n = ", sample_counts, ")")
      
      # Plot raw data for the current comparison method
      if (method == "PCA") {
        pca_raw <- prcomp(pure_data, center = TRUE, scale. = TRUE)
        
        if (color_donor_only) {
          plot <- autoplot(pca_raw, data = merged_data, colour = 'Donor', frame = TRUE, frame.type = 'norm') +
            labs(title = paste0("PCA (Raw Data), ", Group)) +
            theme_minimal()
        
        } else if (Group %in% c("Timepoint", "Donor")) {
          plot <- autoplot(pca_raw, data = merged_data, colour = Group, frame = TRUE, frame.type = 'norm') +
            labs(title = paste0("PCA (Raw Data), ", Group)) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")
            )
        } else {
          plot <- autoplot(pca_raw, data = merged_data, colour = Group, shape = 'Timepoint', frame = TRUE, frame.type = 'norm') +
            labs(title = paste0("PCA (Raw Data), ", Group, " (Shaped by Timepoint)")) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")
            )
        }
        plot_list[[length(plot_list) + 1]] <- plot
        
      } else if (method == "PLSDA") {
        plsda_raw <- caret::plsda(pure_data, as.factor(merged_data[[Group]]), ncomp = 2)
        plsda_scores_raw <- data.frame(Comp.1 = numeric(nrow(merged_data)),
                                       Comp.2 = numeric(nrow(merged_data)),
                                       Group = factor(nrow(merged_data)),
                                       stringsAsFactors = FALSE)
        plsda_scores_raw$Comp.1 <- plsda_raw$scores[,1]
        plsda_scores_raw$Comp.2 <- plsda_raw$scores[,2]
        plsda_scores_raw$Group <- merged_data[[Group]]
        plsda_scores_raw$Timepoint <- as.factor(merged_data$Timepoint)
        plsda_scores_raw$Donor <- as.factor(merged_data$Donor)
        
        # Calculate the comp%
        # Get scores matrix (n x 2)
        scores <- plsda_raw$scores
        
        # Total variance in X
        total_var <- sum(apply(pure_data, 2, var))
        
        # Variance of each PLS component (across all samples)
        comp_var <- apply(scores, 2, var)
        
        # Percent variance explained by each component
        expl_var <- 100 * comp_var / total_var
        
        xlab <- sprintf("Comp. 1 (%.1f%%)", expl_var[1])
        ylab <- sprintf("Comp. 2 (%.1f%%)", expl_var[2])
        
        if (color_donor_only) {
          plot <- ggplot(plsda_scores_raw, aes(x = Comp.1, y = Comp.2, color = !!sym("Donor"))) +
            geom_point(size = 3) +
            labs(title = paste0("PLSDA (Raw Data), ", Group),
                 x = xlab,
                 y = ylab) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")
            )
          
        } else if (Group %in% c("Timepoint", "Donor")) {
          plot <- ggplot(plsda_scores_raw, aes(x = Comp.1, y = Comp.2, color = !!sym("Group"))) +
            geom_point(size = 3) +
            labs(title = paste0("PLSDA (Raw Data), ", Group),
                 x = xlab,
                 y = ylab) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")
            )
        } else {
          plot <- ggplot(plsda_scores_raw, aes(x = Comp.1, y = Comp.2, color = !!sym("Group"), shape = !!sym("Timepoint"))) +
            geom_point(size = 3) +
            labs(title = paste0("PLSDA (Raw Data), ", Group, " (Shaped by Timepoint)"),
                 x = xlab,
                 y = ylab) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),
              axis.text = element_text(size = 12),
              plot.title = element_text(size = 16, face = "bold")
            )
          
        }
        plot_list[[length(plot_list) + 1]] <- plot
        
      } else if (method == "RLA") {
        rla_data <- t(pure_data) # need to transform to get features to rows
        rla_raw <- rowRla(rla_data, f = merged_data[[Group]]) #makes rla row wise
        rla_long <- as.data.frame(rla_raw) %>%
          rownames_to_column(var = "Feature") %>%
          pivot_longer(cols = -Feature, names_to = "Sample", values_to = "Value")

        plot <- ggplot(rla_long, aes(x = Sample, y = Value)) +
          geom_boxplot(outlier.size = 0.1) +
          labs(title = paste0("RLA (Raw Data), ", Group), x = "Samples", y = "Relative Log Abundance") +
          theme_minimal() +
          theme(axis.text.x = element_blank(),
                axis.title = element_text(size = 14),
                axis.text = element_text(size = 12),
                plot.title = element_text(size = 16, face = "bold"))
        plot_list[[length(plot_list) + 1]] <- plot

      }
    }
    
    # Combine all plots for the current comparison method into one plot
    do.call(grid.arrange, c(plot_list, ncol = 2))
  }
  
  # Close PDF device
  dev.off()
}





stat_group_comparison <- function(data, resultsdir, group_var, variable_list, get_distinct_data = TRUE, name = 'statistical_tests_results') {
  results <- data.frame(
    Variable = character(),
    `Variance Test (p)` = character(),
    `Test Type` = character(),
    `Mean (M) \\ Mean (F)` = character(),
    `t (df)` = character(),
    `p-value` = character(),
    `Normality (p)` = character(),
    Interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  if (get_distinct_data) {
    distinct_data <- data %>% distinct(Donor, .keep_all = TRUE)
  } else {
    distinct_data <- data
  }
  
  
  for (var in variable_list) {
    df <- distinct_data[, c(var, group_var)]
    df <- df[complete.cases(df), ]
    group_levels <- unique(df[[group_var]])
    
    if (length(group_levels) != 2) {
      warning(paste("Skipping", var, "- needs exactly two groups"))
      next
    }
    
    g1 <- df[df[[group_var]] == group_levels[1], var]
    g2 <- df[df[[group_var]] == group_levels[2], var]
    
    
    # CHECK IF BOTH g1 and g2 ARE NUMERIC
    is_numeric <- is.numeric(g1) && is.numeric(g2)
    
    if (is_numeric) {
      # For numerical groups t-test
      var_p <- tryCatch(var.test(g1, g2)$p.value, error = function(e) NA)
      var_p_formatted <- ifelse(var_p < 0.001, format(var_p, scientific = TRUE, digits = 2), round(var_p, 3))
      test_type <- ifelse(!is.na(var_p) && var_p > 0.05, "t-test", "Welch")
      
      ttest <- t.test(g1, g2, var.equal = (test_type == "t-test"))
      t_val <- round(ttest$statistic, 3)
      df_val <- round(ttest$parameter, 2)
      t_df <- paste0(t_val, " (", df_val, ")")
      p_val <- ifelse(ttest$p.value < 0.001, format(ttest$p.value, scientific = TRUE, digits = 2),
                      format(round(ttest$p.value, 4), nsmall = 4))
      
      m1 <- round(mean(g1), 2)
      m2 <- round(mean(g2), 2)
      mean_string <- paste0("\\makecell{", m1, " \\\\ ", m2, "}")
      
      shap1 <- if (length(g1) >= 3) signif(shapiro.test(g1)$p.value, 3) else NA
      shap2 <- if (length(g2) >= 3) signif(shapiro.test(g2)$p.value, 3) else NA
      normality <- paste0("\\makecell{", 
                          ifelse(!is.na(shap1), paste0(shap1, " (", group_levels[1], ")"), "-"), 
                          " \\\\ ", 
                          ifelse(!is.na(shap2), paste0(shap2, " (", group_levels[2], ")"), "-"), 
                          "}")
      if (ttest$p.value < 0.05) {
        direction <- ifelse(m1 > m2, paste("Higher in", group_levels[1]), paste("Higher in", group_levels[2]))
      } else {
        direction <- "No significant difference"
      }
      normality_note <- if (all(c(shap1, shap2) > 0.05, na.rm = TRUE)) "Follows normal distribution" else "Non-normal"
      interpretation <- paste0("\\makecell{", direction, " \\\\ ", normality_note, "}")
      
      
    } else {
      # For categorical variable, Fisher's Exact Test (2x2) 
      test_type <- "Fisher"
      # Build the 2x2 table
      tbl <- table(df[[var]], df[[group_var]])
      
      if (nrow(tbl) != 2 || ncol(tbl) != 2) {
        warning(paste("Skipping", var, "- Fisher's test needs exactly two levels in both variable and group"))
        next
      }
      
      fisher <- fisher.test(tbl)
      var_p_formatted <- "-"
      t_df <- "-"
      m1 <- paste0(sum(df[[var]] == rownames(tbl)[1] & df[[group_var]] == group_levels[1]), 
                   " / ", sum(df[[group_var]] == group_levels[1]))
      m2 <- paste0(sum(df[[var]] == rownames(tbl)[1] & df[[group_var]] == group_levels[2]), 
                   " / ", sum(df[[group_var]] == group_levels[2]))
      mean_string <- paste0("\\makecell{", m1, " \\\\ ", m2, "}")
      
      p_val <- ifelse(fisher$p.value < 0.001, format(fisher$p.value, scientific=TRUE, digits=2), 
                      format(round(fisher$p.value, 4), nsmall=4))
      normality <- "-"
      # Interpret "higher" as higher proportion
      prop1 <- sum(df[[var]] == rownames(tbl)[1] & df[[group_var]] == group_levels[1]) / sum(df[[group_var]] == group_levels[1])
      prop2 <- sum(df[[var]] == rownames(tbl)[1] & df[[group_var]] == group_levels[2]) / sum(df[[group_var]] == group_levels[2])
      if (fisher$p.value < 0.05) {
        direction <- ifelse(prop1 > prop2, paste("More", rownames(tbl)[1], "in", group_levels[1]), 
                            paste("More", rownames(tbl)[1], "in", group_levels[2]))
      } else {
        direction <- "No significant difference"
      }
      interpretation <- paste0("\\makecell{", direction, " \\\\ ", "Not tested" , "}")
    }
    
    # Add row
    results <- rbind(results, data.frame(
      Variable = var,
      `Variance Test (p)` = var_p_formatted,
      `Test Type` = test_type,
      `Mean (M) \\ Mean (F)` = mean_string,
      `t (df)` = t_df,
      `p-value` = p_val,
      `Normality (p)` = normality,
      Interpretation = interpretation,
      stringsAsFactors = FALSE
    ))
  }
    
    # Write CSV
    write.csv(results, paste0(resultsdir, name, '.csv'), row.names = FALSE)
    return(results)
  }




# Function: Plot Tyrosine to Tryptophane over time
# Input: resultsdir (directory for result files), targeted_experiment_data, add_info, list of Groups
# Output: PDF with plots
plot_to_pdf_Tyr_to_Trp <- function(resultsdir, targeted_experiment_data, add_info, Group_list) {
  
  Donors <- unique(add_info$Donor)
  
  needed_data <- targeted_experiment_data[grepl("Tyrosine|Tryptophan", targeted_experiment_data$Molecule.Name, ignore.case = TRUE), ]
  needed_data$Molecule.Name <- sub("^DL-", "", needed_data$Molecule.Name) # Remove "DL-" if present
  rownames(needed_data) <- needed_data$Molecule.Name
  needed_data <- t(needed_data) 
  needed_data <- as.data.frame(needed_data)
  needed_data <- needed_data[-c(1:4, nrow(needed_data)),]
  needed_data$`Tyrosine`<- as.numeric(needed_data$`Tyrosine`)
  needed_data$`Tryptophan`<- as.numeric(needed_data$`Tryptophan`)
  needed_data$Tyr.to.Trp <- needed_data$`Tyrosine` / needed_data$`Tryptophan`
  needed_data$`Tyrosine.log2`<- log2(needed_data$`Tyrosine`)
  needed_data$`Tryptophan.log2`<- log2(needed_data$`Tryptophan`)
  needed_data$Sample <- rownames(needed_data)
  
  merged_data <- merge(needed_data, add_info, by.x = "Sample", by.y = "Sample")
  merged_data <- merged_data %>%   # Calculate logFC between consecutive time points for each donor
    group_by(Donor, Intervention) %>%
    arrange(Timepoint) %>%
    mutate(
      logFC_Tyr = log2(`Tyrosine` / lag(`Tyrosine`)),
      logFC_Trp = log2(`Tryptophan` / lag(`Tryptophan`))
    ) %>%
    ungroup()
  
  
  # Set up PDF device
  pdf(file = paste0(resultsdir, "Tyrosine_to_Tryptophane_comparison_plots.pdf"), width = 20, height = 6 * (length(Donors) %/% 4 + 1)) 
  
  
  # Plot the Tyrosine to Tryptophan ratio for each donor
  plot_list <- list()

  for (donor in Donors) {
    filtered_data <- merged_data %>% filter(Donor == donor)
    
    plot <- ggplot(filtered_data, aes(x = Timepoint, y = Tyr.to.Trp, color = Intervention, group = Intervention)) +
            geom_line() +               
            geom_point() +              
            labs(x = "Time Point", y = "Tyrosine/Tryptophan Ratio", title = paste0("Tyr/Trp Ratio over Time by Intervention  \nDonor ", donor)) +
            theme_minimal() +           
            scale_color_manual(values = c("high c/p" = "blue", "low c/p" = "green")) +
      theme(
        axis.title = element_text(size = 14),   # axis labels
        axis.text = element_text(size = 12)    # tick labels
      )
    
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  # call grid arrange later to plot to page after the overview plots
  # yet need it to be able to see the scaling for the overview plots
  
  
  
  # Plot data per Interventions and stratify it for some plots
  plot_list2 <- list()
  
  
  # Plot the mean Tyrosine to Tryptophan ratio for the Interventions
    for(Group in Group_list) {
      merged_data[[Group]] <- as.factor(merged_data[[Group]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group)) %>%
        summarise(
          Median = median(Tyr.to.Trp, na.rm = TRUE),
          SD = sd(Tyr.to.Trp, na.rm = TRUE)
        )
      
      plot <- ggplot(summary_data, aes(x = Timepoint, y = Median, color = !!sym(Group))) +
        geom_point(size = 4) +  # Median points
        geom_errorbar(aes(ymin = Median - SD, ymax = Median + SD), width = 0.2) +  # Error bars
        geom_line() + 
        labs(
          x = "Time Point",
          y = "Tyr/Trp Ratio",
          title = paste("Tyr/Trp Ratio over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list2[[length(plot_list2) + 1]] <- plot
    }
    
    
    if (length(Group_list) > 1) {
      merged_data[[Group_list[[1]]]] <- as.factor(merged_data[[Group_list[[1]]]])
      merged_data[[Group_list[[2]]]] <- as.factor(merged_data[[Group_list[[2]]]])
      
      # Summarize the data with two grouping variables
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%
        summarise(
          Median = median(Tyr.to.Trp, na.rm = TRUE),
          SD = sd(Tyr.to.Trp, na.rm = TRUE)
        )
      
      # Create the plot with two grouping variables
      plot <- ggplot(summary_data, aes(x = Timepoint, y = Median, color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +  # Median points
        geom_errorbar(aes(ymin = Median - SD, ymax = Median + SD), width = 0.2) +  # Error bars
        geom_line() + 
        labs(
          x = "Time Point",
          y = "Tyr/Trp Ratio",
          title = paste("Tyr/Trp Ratio over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list2[[length(plot_list2) + 1]] <- plot
      
      
      plot_ratio <- list()
      
      for(Group in unique(merged_data[[Group_list[[1]]]])) {
        filtered_data <- merged_data %>%
          filter(!!sym(Group_list[[1]]) == Group) %>%
          filter(!!sym(Group_list[[2]]) == unique(merged_data[[Group_list[[2]]]])[1])
        
        summary_data <- filtered_data %>%
          group_by(Timepoint) %>%
          summarise(
            Median = median(Tyr.to.Trp, na.rm = TRUE),
            SD = sd(Tyr.to.Trp, na.rm = TRUE)
          )
        
        # Create the plot
        plot <- ggplot() +
          # Plot individual donor data as light gray lines
          geom_line(data = filtered_data, aes(x = Timepoint , y = Tyr.to.Trp, group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          # Plot summary data (median and error bars) in blue
          geom_line(data = summary_data, aes(x = Timepoint , y = Median), color = "blue", linewidth = 1) +
          #geom_errorbar(data = summary_data, aes(x = Timepoint , ymin = Median - SD, ymax = Median + SD), width = 0.2, color = "blue") +  # Error bars
          # Labels and theme
          labs(
            x = "Time Point",
            y = "Tyr/Trp Ratio",
            title = paste("Tyr/Trp Ratio over Time for", Group, unique(merged_data[[Group_list[[2]]]])[1])
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_ratio[[length(plot_ratio) + 1]] <- plot
        
        
        filtered_data <- merged_data %>%
          filter(!!sym(Group_list[[1]]) == Group) %>%
          filter(!!sym(Group_list[[2]]) == unique(merged_data[[Group_list[[2]]]])[2])
        
        summary_data <- filtered_data %>%
          group_by(Timepoint) %>%
          summarise(
            Median = median(Tyr.to.Trp, na.rm = TRUE),
            SD = sd(Tyr.to.Trp, na.rm = TRUE)
          )
        
        # Create the plot
        plot <- ggplot() +
          # Plot individual donor data as light gray lines
          geom_line(data = filtered_data, aes(x = Timepoint , y = Tyr.to.Trp, group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          # Plot summary data (median and error bars) in blue
          geom_line(data = summary_data, aes(x = Timepoint , y = Median), color = "blue", linewidth = 1) +
          #geom_errorbar(data = summary_data, aes(x = Timepoint , ymin = Median - SD, ymax = Median + SD), width = 0.2, color = "blue") +  # Error bars
          # Labels and theme
          labs(
            x = "Time Point",
            y = "Tyr/Trp Ratio",
            title = paste("Tyr/Trp Ratio over Time for", Group, unique(merged_data[[Group_list[[2]]]])[2])
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_ratio[[length(plot_ratio) + 1]] <- plot
      }
    }
  
  
  # Plot the mean Tyrosine and Tryptophan values for the Interventions
    plot_list_Tyr_error_bars <- list()
    plot_list_Trp_error_bars <- list()
    plot_list_Tyr <- list()
    plot_list_Trp <- list()
    
    for(Group in Group_list) {
      merged_data[[Group]] <- as.factor(merged_data[[Group]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group)) %>%
        summarise(
          Median_Tyr = median(!!sym('Tyrosine.log2'), na.rm = TRUE),
          SD_Tyr = sd(!!sym('Tyrosine.log2'), na.rm = TRUE),
          Median_Trp = median(!!sym('Tryptophan.log2'), na.rm = TRUE),
          SD_Trp = sd(!!sym('Tryptophan.log2'), na.rm = TRUE)
        )
      
      # Plot Tyrosine
      plot_Tyr <- ggplot(summary_data, aes(x = Timepoint , y = Median_Tyr, color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_Tyr - SD_Tyr, ymax = Median_Tyr + SD_Tyr), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2(Area)",
          title = paste("Tyrosine over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Tyr_error_bars[[length(plot_list_Tyr_error_bars) + 1]] <- plot_Tyr
      
      
      # Plot Tryptophan
      plot_Trp <- ggplot(summary_data, aes(x = Timepoint , y = Median_Trp, color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_Trp - SD_Trp, ymax = Median_Trp + SD_Trp), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2(Area)",
          title = paste("Tryptophan over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Trp_error_bars[[length(plot_list_Trp_error_bars) + 1]] <- plot_Trp
    }
    
    
    if (length(Group_list) > 1) {
      merged_data[[Group_list[[1]]]] <- as.factor(merged_data[[Group_list[[1]]]])
      merged_data[[Group_list[[2]]]] <- as.factor(merged_data[[Group_list[[2]]]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%
        summarise(
          Median_Tyr = median(!!sym('Tyrosine.log2'), na.rm = TRUE),
          SD_Tyr = sd(!!sym('Tyrosine.log2'), na.rm = TRUE),
          Median_Trp = median(!!sym('Tryptophan.log2'), na.rm = TRUE),
          SD_Trp = sd(!!sym('Tryptophan.log2'), na.rm = TRUE)
        )
      
      # Plot Tyrosine
      plot_Tyr <- ggplot(summary_data, aes(x = Timepoint , y = Median_Tyr, color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_Tyr - SD_Tyr, ymax = Median_Tyr + SD_Tyr), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2(Area)",
          title = paste("Tyrosine over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Tyr_error_bars[[length(plot_list_Tyr_error_bars) + 1]] <- plot_Tyr
      
      
      # Plot Tryptophan
      plot_Trp <- ggplot(summary_data, aes(x = Timepoint , y = Median_Trp, color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_Trp - SD_Trp, ymax = Median_Trp + SD_Trp), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2(Area)",
          title = paste("Tryptophan over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Trp_error_bars[[length(plot_list_Trp_error_bars) + 1]] <- plot_Trp
    }
    
    
    for(Group in unique(merged_data[[Group_list[[1]]]])) {
      for (Subgroup in unique(merged_data[[Group_list[[2]]]])) {
        filtered_data <- merged_data %>%
          filter(!!sym(Group_list[[1]]) == Group) %>%
          filter(!!sym(Group_list[[2]]) == Subgroup)
        
        summary_data <- filtered_data %>%
          group_by(Timepoint) %>%
          summarise(
            Median_Tyr = median(!!sym('Tyrosine.log2'), na.rm = TRUE),
            SD_Tyr = sd(!!sym('Tyrosine.log2'), na.rm = TRUE),
            Median_Trp = median(!!sym('Tryptophan.log2'), na.rm = TRUE),
            SD_Trp = sd(!!sym('Tryptophan.log2'), na.rm = TRUE)
          )
        
        # Plot Tyrosine
        plot_Tyr <- ggplot() +
          geom_line(data = filtered_data, aes(x = Timepoint , y = !!sym('Tyrosine.log2'), group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          geom_line(data = summary_data, aes(x = Timepoint , y = Median_Tyr), color = "blue", linewidth = 1) +
          labs(
            x = "Time Point",
            y = "log2(Area)",
            title = paste("Tyrosine over Time for", Group, Subgroup)
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_list_Tyr[[length(plot_list_Tyr) + 1]] <- plot_Tyr
        
        # Plot Tryptophan
        plot_Trp <- ggplot() +
          geom_line(data = filtered_data, aes(x = Timepoint , y = !!sym('Tryptophan.log2'), group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          geom_line(data = summary_data, aes(x = Timepoint , y = Median_Trp), color = "blue", linewidth = 1) +
          labs(
            x = "Time Point",
            y = "log2(Area)",
            title = paste("Tryptophan over Time for", Group, Subgroup)
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_list_Trp[[length(plot_list_Trp) + 1]] <- plot_Trp
      }
    }
    
    
    plot_list_Tyr_error_bars_logFC <- list()
    plot_list_Trp_error_bars_logFC <- list()
    plot_list_Tyr_logFC <- list()
    plot_list_Trp_logFC <- list()
    
    # Plot the logFC Tyrosine and Tryptophan values for the Interventions
    for(Group in Group_list) {
      merged_data[[Group]] <- as.factor(merged_data[[Group]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group)) %>%
        summarise(
          Median_logFC_Tyr = median(logFC_Tyr, na.rm = TRUE),
          SD_logFC_Tyr = sd(logFC_Tyr, na.rm = TRUE),
          Median_logFC_Trp = median(logFC_Trp, na.rm = TRUE),
          SD_logFC_Trp = sd(logFC_Trp, na.rm = TRUE)
        )
      
      # Plot Tyrosine logFC
      plot_logFC_Tyr <- ggplot(summary_data, aes(x = Timepoint , y = Median_logFC_Tyr, color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_logFC_Tyr - SD_logFC_Tyr, ymax = Median_logFC_Tyr + SD_logFC_Tyr), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2 Fold Change",
          title = paste("Tyrosine log2 Fold Change over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Tyr_error_bars_logFC[[length(plot_list_Tyr_error_bars_logFC) + 1]] <- plot_logFC_Tyr
      
      # Plot Tryptophan logFC
      plot_logFC_Trp <- ggplot(summary_data, aes(x = Timepoint , y = Median_logFC_Trp, color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_logFC_Trp - SD_logFC_Trp, ymax = Median_logFC_Trp + SD_logFC_Trp), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2 Fold Change",
          title = paste("Tryptophan log2 Fold Change over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Trp_error_bars_logFC[[length(plot_list_Trp_error_bars_logFC) + 1]] <- plot_logFC_Trp
    }
    
    
    if (length(Group_list) > 1) {
      merged_data[[Group_list[[1]]]] <- as.factor(merged_data[[Group_list[[1]]]])
      merged_data[[Group_list[[2]]]] <- as.factor(merged_data[[Group_list[[2]]]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%
        summarise(
          Median_logFC_Tyr = median(logFC_Tyr, na.rm = TRUE),
          SD_logFC_Tyr = sd(logFC_Tyr, na.rm = TRUE),
          Median_logFC_Trp = median(logFC_Trp, na.rm = TRUE),
          SD_logFC_Trp = sd(logFC_Trp, na.rm = TRUE)
        )
      
      # Plot Tyrosine logFC
      plot_logFC_Tyr <- ggplot(summary_data, aes(x = Timepoint , y = Median_logFC_Tyr, color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_logFC_Tyr - SD_logFC_Tyr, ymax = Median_logFC_Tyr + SD_logFC_Tyr), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2 Fold Change",
          title = paste("Tyrosine log2 Fold Change over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Tyr_error_bars_logFC[[length(plot_list_Tyr_error_bars_logFC) + 1]] <- plot_logFC_Tyr
      
      # Plot Tryptophan logFC
      plot_logFC_Trp <- ggplot(summary_data, aes(x = Timepoint , y = Median_logFC_Trp, color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_logFC_Trp - SD_logFC_Trp, ymax = Median_logFC_Trp + SD_logFC_Trp), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2 Fold Change",
          title = paste("Tryptophan log2 Fold Change over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_Trp_error_bars_logFC[[length(plot_list_Trp_error_bars_logFC) + 1]] <- plot_logFC_Trp
    }
    
    
    
    for(Group in unique(merged_data[[Group_list[[1]]]])) {
      for (Subgroup in unique(merged_data[[Group_list[[2]]]])) {
        filtered_data <- merged_data %>%
          filter(!!sym(Group_list[[1]]) == Group) %>%
          filter(!!sym(Group_list[[2]]) == Subgroup)
        
        summary_data <- filtered_data %>%
          group_by(Timepoint) %>%
          summarise(
            Median_logFC_Tyr = median(logFC_Tyr, na.rm = TRUE),
            SD_logFC_Tyr = sd(logFC_Tyr, na.rm = TRUE),
            Median_logFC_Trp = median(logFC_Trp, na.rm = TRUE),
            SD_logFC_Trp = sd(logFC_Trp, na.rm = TRUE)
          )
        
        # Plot Tyrosine logFC
        plot_logFC_Tyr <- ggplot() +
          geom_line(data = filtered_data, aes(x = Timepoint , y = logFC_Tyr, group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          geom_line(data = summary_data, aes(x = Timepoint , y = Median_logFC_Tyr), color = "blue", linewidth = 1) +
          labs(
            x = "Time Point",
            y = "log2 Fold Change",
            title = paste("Tyrosine log2 Fold Change over Time for", Group, Subgroup)
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_list_Tyr_logFC[[length(plot_list_Tyr_logFC) + 1]] <- plot_logFC_Tyr
        
        # Plot Tryptophan logFC
        plot_logFC_Trp <- ggplot() +
          geom_line(data = filtered_data, aes(x = Timepoint , y = logFC_Trp, group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
          geom_line(data = summary_data, aes(x = Timepoint , y = Median_logFC_Trp), color = "blue", linewidth = 1) +
          labs(
            x = "Time Point",
            y = "log2 Fold Change",
            title = paste("Tryptophan log2 Fold Change over Time for", Group, Subgroup)
          ) +
          theme_minimal() +
          theme(
            axis.title = element_text(size = 14),   # axis labels
            axis.text = element_text(size = 12)    # tick labels
          )
        
        plot_list_Trp_logFC[[length(plot_list_Trp_logFC) + 1]] <- plot_logFC_Trp
      }
    }
    
  
  # Add plot list together and arrange accordingly to put on one page
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
    plot_list2 <- append(plot_list2, plot_list_Tyr_error_bars)
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
    plot_list2 <- append(plot_list2, plot_list_Trp_error_bars)
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
    plot_list2 <- append(plot_list2, plot_list_Tyr_error_bars_logFC)
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
    plot_list2 <- append(plot_list2, plot_list_Trp_error_bars_logFC)
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
    plot_list2 <- append(plot_list2, plot_ratio)
    plot_list2 <- append(plot_list2, plot_list_Tyr)
    plot_list2 <- append(plot_list2, plot_list_Trp)
    plot_list2 <- append(plot_list2, plot_list_Trp_logFC)
    plot_list2 <- append(plot_list2, plot_list_Trp_logFC)
    
    
  if (length(plot_list) > length(plot_list2)) {
    for (i in 1:(length(plot_list) - length(plot_list2)) ){ #fill in empty plots so the plots dont get stretched
      plot_list2[[length(plot_list2) + 1]] <- ggplot()
    }
  }
  do.call(grid.arrange, c(plot_list2, ncol = 4))
  
  
  do.call(grid.arrange, c(plot_list, ncol = 4))
  
  
  # Plot Tyrosine and Tryptophan as separate lines for each donor
  plot_list <- list()
  
  for (donor in Donors) {
    filtered_data <- merged_data %>% filter(Donor == donor)
    
    plot <- ggplot(filtered_data, aes(x = Timepoint , color = Intervention)) +
      geom_line(aes(y = !!sym('Tyrosine.log2'), linetype = "Tyrosine")) +  # Plot Tyrosine line with a specific linetype
      geom_line(aes(y = !!sym('Tryptophan.log2'), linetype = "Tryptophan")) +  # Plot Tryptophan line with a specific linetype
      geom_point(aes(y = !!sym('Tyrosine.log2'))) +  # Add points for Tyrosine
      geom_point(aes(y = !!sym('Tryptophan.log2'))) +  # Add points for Tryptophan
      labs(x = "Time Point", y = "log2(Area)", title = paste0("Tyrosine and Tryptophan over Time by Intervention  \nDonor ", donor)) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),   # axis labels
        axis.text = element_text(size = 12)    # tick labels
      ) +           
      scale_color_manual(values = c("high c/p" = "blue", "low c/p" = "green")) +  # Use blue and green for interventions
      scale_linetype_manual(values = c("Tyrosine" = "solid", "Tryptophan" = "dashed"))  # Use solid for Tyrosine and dashed for Tryptophan
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  do.call(grid.arrange, c(plot_list, ncol = 4))
  
  
  # Plot Tyrosine and Tryptophan log fold change as separate lines for each donor
  plot_list <- list()
  
  for (donor in Donors) {
    filtered_data <- merged_data %>% filter(Donor == donor)
    
    plot <- ggplot(filtered_data, aes(x = Timepoint , color = Intervention)) +
      geom_line(aes(y = logFC_Tyr, linetype = "Tyrosine")) +  # Plot Tyrosine logFC with a specific linetype
      geom_line(aes(y = logFC_Trp, linetype = "Tryptophan")) +  # Plot Tryptophan logFC with a specific linetype
      geom_point(aes(y = logFC_Tyr)) +  # Add points for Tyrosine logFC
      geom_point(aes(y = logFC_Trp)) +  # Add points for Tryptophan logFC
      labs(x = "Time Point", y = "log2 Fold Change", title = paste0("Log2 Fold Change of Tyr and Trp over Time by Intervention  \nDonor ", donor)) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),   # axis labels
        axis.text = element_text(size = 12)    # tick labels
      ) +           
      scale_color_manual(values = c("high c/p" = "blue", "low c/p" = "green")) +  # Use blue and green for interventions
      scale_linetype_manual(values = c("Tyrosine" = "solid", "Tryptophan" = "dashed"))  # Use solid for Tyrosine and dashed for Tryptophan
    
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  do.call(grid.arrange, c(plot_list, ncol = 4))
  
  
  # Close PDF device
  dev.off()
}



# Function: Plot amino acid ratios and individual values over time, stratified by groups
plot_to_pdf_Tyr_Trp_to_LNAA <- function(resultsdir, targeted_experiment_data, add_info, Group_list) {
  
  Donors <- unique(add_info$Donor)
  amino_acids <- c("Phenylalanine", "Tyrosine", "Tryptophan", "Leucine", "Isoleucine", "Valine", "Methionine", "Histidine", "Threonine")
  ratios <- c("Tyr_to_LNAA_excl_Tyr", "Trp_to_LNAA_excl_Trp", "Tyr_Trp_to_LNAA_excl_both")
  
  # Filter for required amino acids
  needed_data <- targeted_experiment_data[grepl(paste0("(", paste(amino_acids, collapse = "|"), ")$"), targeted_experiment_data$Molecule.Name, ignore.case = TRUE), ]
  needed_data$Molecule.Name <- sub("^DL-", "", needed_data$Molecule.Name) # Remove "DL-" if present
  rownames(needed_data) <- needed_data$Molecule.Name
  needed_data <- t(needed_data)
  needed_data <- as.data.frame(needed_data)
  needed_data <- needed_data[-c(1:4, nrow(needed_data)), ]
  
  # Convert amino acid values to numeric
  for (aa in amino_acids) {
    needed_data[[aa]] <- as.numeric(needed_data[[aa]])
  }
  
  # Calculate ratios
  needed_data$LNAA_sum_excl_Tyr <- rowSums(needed_data[, amino_acids[!amino_acids %in% "Tyrosine"]], na.rm = TRUE)
  needed_data$Tyr_to_LNAA_excl_Tyr <- needed_data$Tyrosine / needed_data$LNAA_sum_excl_Tyr
  needed_data$LNAA_sum_excl_Trp <- rowSums(needed_data[, amino_acids[!amino_acids %in% "Tryptophan"]], na.rm = TRUE)
  needed_data$Trp_to_LNAA_excl_Trp <- needed_data$Tryptophan / needed_data$LNAA_sum_excl_Trp
  needed_data$LNAA_sum_excl_Tyr_Trp <- rowSums(needed_data[, amino_acids[!amino_acids %in% c("Tyrosine", "Tryptophan")]], na.rm = TRUE)
  needed_data$Tyr_Trp_to_LNAA_excl_both <- (needed_data$Tyrosine + needed_data$Tryptophan) / needed_data$LNAA_sum_excl_Tyr_Trp
  
  # Log2 transformation
  for (aa in amino_acids) {
    needed_data[[paste0(aa, ".log2")]] <- log2(needed_data[[aa]])
  }
  needed_data$Sample <- rownames(needed_data)
  
  # Add/create additional information
  merged_data <- merge(needed_data, add_info, by.x = "Sample", by.y = "Sample")
  merged_data <- merged_data %>%   
    group_by(Donor, Intervention) %>%
    arrange(Timepoint) %>%
    mutate(across(all_of(amino_acids), list(logFC = ~log2(. / lag(.))), .names = "logFC_{col}")) %>%
    ungroup()
  
  write.csv(needed_data, paste0(resultsdir, "LNAA_data_and_ratios.csv"))
  write.csv(merged_data, paste0(resultsdir, "LNAA_data_and_ratios_with_metadata.csv"))
  
  
  # Set up PDF device
  pdf(file = paste0(resultsdir, "Tyrosine_and_Tryptophan_to_LNAA_comparison_plots.pdf"), width = 24, height = 5 * (length(Donors) %/% 4 + 1))
  
  # Stratified plots by group
  plot_list2 <- list()
  for (Group in Group_list) {
    merged_data[[Group]] <- as.factor(merged_data[[Group]])
    
    # Calculate statistics
    summary_data <- merged_data %>%
      group_by(Timepoint, !!sym(Group)) %>%
      summarise(across(c(Tyr_to_LNAA_excl_Tyr, Trp_to_LNAA_excl_Trp, Tyr_Trp_to_LNAA_excl_both), list(Median = ~median(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE))))
    
    # Loop over each ratio and create plots
    for (ratio in ratios) {
      plot <- ggplot(summary_data, aes(x = Timepoint , y = !!sym(paste0(ratio, "_Median")), color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = !!sym(paste0(ratio, "_Median")) - !!sym(paste0(ratio, "_SD")), 
                          ymax = !!sym(paste0(ratio, "_Median")) + !!sym(paste0(ratio, "_SD"))), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "Ratio",
          title = paste(ratio, "over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list2[[length(plot_list2) + 1]] <- plot
    }
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
  }
  
  
  if (length(Group_list) > 1) {
    merged_data[[Group_list[[1]]]] <- as.factor(merged_data[[Group_list[[1]]]])
    merged_data[[Group_list[[2]]]] <- as.factor(merged_data[[Group_list[[2]]]])
    
    summary_data <- merged_data %>%
      group_by(Timepoint, !!sym(Group_list[[1]]), !!sym(Group_list[[2]])) %>%
      summarise(across(ratios, list(Median = ~median(., na.rm = TRUE), SD = ~sd(., na.rm = TRUE))))
    
    for (ratio in ratios) {
      plot <- ggplot(summary_data, aes(x = Timepoint , y = !!sym(paste0(ratio, "_Median")), color = !!sym(Group_list[[1]]), linetype = !!sym(Group_list[[2]]))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = !!sym(paste0(ratio, "_Median")) - !!sym(paste0(ratio, "_SD")), 
                          ymax = !!sym(paste0(ratio, "_Median")) + !!sym(paste0(ratio, "_SD"))), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "Ratio",
          title = paste(ratio, "over Time by", Group_list[[1]], "and", Group_list[[2]])
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list2[[length(plot_list2) + 1]] <- plot
    }
    plot_list2[[length(plot_list2) + 1]] <- ggplot()
  
    # Additional plots for ratios, shows each donor and the mean for groups
    for (ration in ratios) {
      for (Group in unique(merged_data[[Group_list[[1]]]])) {
        for (Subgroup in unique(merged_data[[Group_list[[2]]]])) {
          filtered_data <- merged_data %>%
            filter(!!sym(Group_list[[1]]) == Group) %>%
            filter(!!sym(Group_list[[2]]) == Subgroup)
          
          # Summary statistics for each ratio
          summary_data <- filtered_data %>%
            group_by(Timepoint) %>%
            summarise(Median = median(get(ratio), na.rm = TRUE),
                      SD = sd(get(ratio), na.rm = TRUE))
          
          # Plot each ratio
          plot_ratio <- ggplot() +
            geom_line(data = filtered_data, aes(x = Timepoint , y = get(ratio), group = Donor), color = "gray", alpha = 0.5, linewidth = 0.5) +
            geom_line(data = summary_data, aes(x = Timepoint , y = Median), color = "blue", linewidth = 1) +
            labs(
              x = "Time Point",
              y = "Ratio",
              title = paste(ratio, "over Time by", Group, Subgroup)
            ) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),   # axis labels
              axis.text = element_text(size = 12)    # tick labels
            )
          
          plot_list2[[length(plot_list2) + 1]] <- plot_ratio
        }
      }
    }
   
    do.call(grid.arrange, c(plot_list2, ncol = 4)) 
  }
  
  
  plot_list_aa <- list()
  # Plot each amino acid stratified over groups
  if (length(Group_list) > 1) {
    for (aa in amino_acids) {
      for (Group in unique(merged_data[[Group_list[[1]]]])) {
        for (Subgroup in unique(merged_data[[Group_list[[2]]]])) {
          filtered_data <- merged_data %>%
            filter(!!sym(Group_list[[1]]) == Group) %>%
            filter(!!sym(Group_list[[2]]) == Subgroup)
          
          # Summary statistics for each amino acid
          summary_data <- filtered_data %>%
            group_by(Timepoint) %>%
            summarise(Median_log2 = median(get(paste0(aa, ".log2")), na.rm = TRUE),
                      SD_log2 = sd(get(paste0(aa, ".log2")), na.rm = TRUE))
          
          # Plot each amino acid
          plot_combined_aa <- ggplot() +
            geom_line(data = filtered_data, aes(x = Timepoint , y = !!sym(paste0(aa, ".log2")), group = Donor), color = "gray", alpha = 1, linewidth = 0.5) +
            geom_line(data = summary_data, aes(x = Timepoint , y = Median_log2), color = "blue", linewidth = 1) +
            labs(
              x = "Time Point",
              y = "log2(Area)",
              title = paste(aa, "over Time by", Group, Subgroup)
            ) +
            theme_minimal() +
            theme(
              axis.title = element_text(size = 14),   # axis labels
              axis.text = element_text(size = 12)    # tick labels
            )
          
          plot_list_aa[[length(plot_list_aa) + 1]] <- plot_combined_aa
        }
      }
    }
    do.call(grid.arrange, c(plot_list_aa, ncol = 4))
  }

  
  
  # Fold change plots for each amino acid
  plot_list_fc <- list()
  for (aa in amino_acids) {
    for (Group in Group_list) {
      merged_data[[Group]] <- as.factor(merged_data[[Group]])
      
      summary_data <- merged_data %>%
        group_by(Timepoint, !!sym(Group)) %>%
        summarise(Median_logFC = median(get(paste0("logFC_", aa)), na.rm = TRUE),
                  SD_logFC = sd(get(paste0("logFC_", aa)), na.rm = TRUE))
      
      plot <- ggplot(summary_data, aes(x = Timepoint , y = Median_logFC, color = !!sym(Group))) +
        geom_point(size = 4) +
        geom_errorbar(aes(ymin = Median_logFC - SD_logFC, ymax = Median_logFC + SD_logFC), width = 0.2) +
        geom_line() + 
        labs(
          x = "Time Point",
          y = "log2 Fold Change",
          title = paste(aa, "log2 Fold Change over Time by", Group)
        ) +
        theme_minimal() +
        theme(
          axis.title = element_text(size = 14),   # axis labels
          axis.text = element_text(size = 12)    # tick labels
        )
      
      plot_list_fc[[length(plot_list_fc) + 1]] <- plot
    }
  }
  
  do.call(grid.arrange, c(plot_list_fc, ncol = 2))
  
  # Plot for each donor
  plot_list <- list()
  for (donor in Donors) {
    filtered_data <- merged_data %>% filter(Donor == donor)
    
    plot <- ggplot(filtered_data, aes(x = Timepoint , linetype = !!sym(Group_list[[1]]))) + 
      geom_line(aes(y = Tyr_to_LNAA_excl_Tyr, color = "Tyr to LNAAs excl. Tyr")) +
      geom_line(aes(y = Trp_to_LNAA_excl_Trp, color = "Trp to LNAAs excl. Trp"))+
      geom_line(aes(y = Tyr_Trp_to_LNAA_excl_both, color = "Tyr + Trp to LNAAs excl. both")) +
      geom_point(aes(y = Tyr_to_LNAA_excl_Tyr, color = "Tyr to LNAAs excl. Tyr")) + 
      geom_point(aes(y = Trp_to_LNAA_excl_Trp, color = "Trp to LNAAs excl. Trp")) +
      geom_point(aes(y = Tyr_Trp_to_LNAA_excl_both, color = "Tyr + Trp to LNAAs excl. both")) +
      labs(x = "Time Point", y = "Ratio", title = paste0("Ratios over Time by Intervention  \nDonor ", donor)) +
      theme_minimal() +
      theme(
        axis.title = element_text(size = 14),   # axis labels
        axis.text = element_text(size = 12)    # tick labels
      ) +
      scale_color_manual(values = c("Tyr to LNAAs excl. Tyr" = "blue", "Trp to LNAAs excl. Trp" = "green", "Tyr + Trp to LNAAs excl. both" = "red"))
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  do.call(grid.arrange, c(plot_list, ncol = 3))
  
  # Close PDF device
  dev.off()
}




# Helper Function for Volcano plots
volcano_plot <- function(fit, data, title = "Volcano Plot", subtitel = "End - Start") {
  # Extract results
  results <- topTable(fit, number=Inf)
  
  # Plot
  ggplot(results, aes(x=logFC, y=-log10(adj.P.Val))) +
    geom_point(alpha=0.5, color='blue') +
    theme_minimal() +
    labs(title=paste0(title,"\n",subtitel), x='Log Fold Change', y='-Log10(adj. p-value)') +
    geom_hline(yintercept=-log10(0.05), color='red', linetype='dashed') +
    geom_vline(xintercept=c(-1, 1), color='red', linetype='dashed')
}



# Helper Function for beter Volcano plots
enhanced_volcano_plot <- function(fit, data, title = "Volcano Plot", subtitel = "End - Start", adj_p_val, limma_results, xrange = NULL, yrange = NULL, nr_ann_features = 20, ann_features = NULL) {
  # Extract results and make new column with adjusted p.values
  results <- topTable(fit, number=Inf, adjust.method = "BH") # adjust with Benjamini & Hochberg (is standard option)
                                                             # tested all variants and this is the only one that adjust the data not completely into insignificance
  
  results$id <- rownames(results)
  results <- results %>%
    left_join(data %>% mutate(id = as.character(id)), by = "id") %>%
    mutate(Annotation = coalesce(Annotation, paste0(rt, '@', mz))) # if Annotation is NA, put in rt@mz

  # Clean the title
  clean_title <- gsub("[\n(),]", "", title)  # Remove \n, ',', ( and )
  clean_title <- gsub(" ", "-", clean_title)  # Replace spaces with -
  clean_title <- gsub("[/]", ".", clean_title)  # Replace / with .
  
  # Clean the subtitle
  clean_subtitel <- gsub(" ", "", subtitel)  # Remove spaces
  clean_subtitel <- gsub("[\n(),]", "", clean_subtitel)  # Remove \n, ',', ( and )
  clean_subtitel <- gsub("[/]", ".", clean_subtitel)  # Replace / with .
  
  fwrite(results, paste0(limma_results, clean_title, '_' ,clean_subtitel,'.csv'), row.names = FALSE)
  
  
  if (is.null(ann_features)) {
    if (nr_ann_features > 0) {
      top_hits <- results[order(abs(-log10(results$adj.P.Val) * results$logFC), decreasing = TRUE), ][1:nr_ann_features, "Annotation"]
      labsize = 6
    } else {
      top_hits <- list()
      labsize = 4
    }
  } else {
    top_hits <- results$Annotation[
      results$Annotation %in% ann_features & 
        results$adj.P.Val < 0.05 & 
        (results$logFC < -1 | results$logFC > 1)
    ]
    labsize = 7
  }
  
  
  # Plot enhanced Volcano, with adjusted or unadjusted p-value
  if (is.null(xrange) | is.null(yrange)) {
    if(adj_p_val) {
      EnhancedVolcano(results, 
                      lab = ifelse(results$Annotation %in% top_hits, results$Annotation, ""),
                      x = 'logFC',
                      y = 'adj.P.Val', # pick adjusted p values for y axis
                      ylab = bquote(~-Log[10] ~ "adj. P value"),
                      pCutoff = 0.05,
                      labSize = labsize,
                      drawConnectors = TRUE,
                      widthConnectors = 1,
                      arrowheads = FALSE,
                      max.overlaps = 60,
                      title = title,
                      subtitle = subtitel,
                      axisLabSize = 24)
    } else {
      EnhancedVolcano(results, 
                      lab = ifelse(results$Annotation %in% top_hits, results$Annotation, ""),
                      x = 'logFC',
                      y = 'P.Value', # pick non-adjusted p values for y axis
                      ylab = bquote(~-Log[10] ~ "P value"),
                      pCutoff = 0.05,
                      labSize = labsize,
                      drawConnectors = TRUE,
                      widthConnectors = 1,
                      arrowheads = FALSE,
                      max.overlaps = 60,
                      title = title,
                      subtitle = subtitel,
                      axisLabSize = 24)
    }
    
  } else {
    if(adj_p_val) {
      EnhancedVolcano(results, 
                      lab = ifelse(results$Annotation %in% top_hits, results$Annotation, ""),
                      x = 'logFC',
                      y = 'adj.P.Val', # pick adjusted p values for y axis
                      ylab = bquote(~-Log[10] ~ "adj. P value"),
                      xlim = xrange,
                      ylim = yrange,
                      pCutoff = 0.05,
                      labSize = labsize,
                      drawConnectors = TRUE,
                      widthConnectors = 1,
                      arrowheads = FALSE,
                      max.overlaps = 60,
                      title = title,
                      subtitle = subtitel,
                      axisLabSize = 24)
    } else {
      EnhancedVolcano(results, 
                      lab = ifelse(results$Annotation %in% top_hits, results$Annotation, ""),
                      x = 'logFC',
                      y = 'P.Value', # pick non-adjusted p values for y axis
                      ylab = bquote(~-Log[10] ~ "P value"),
                      xlim = xrange,
                      ylim = yrange,
                      pCutoff = 0.05,
                      labSize = labsize,
                      drawConnectors = TRUE,
                      widthConnectors = 1,
                      arrowheads = FALSE,
                      max.overlaps = 60,
                      title = title,
                      subtitle = subtitel,
                      axisLabSize = 24)
    }
  }
}



# Function: Test data with limma and plot Volcanos to pdf
# Input: resultsdir (directory for result files), all_needed_features, add_info, 
#        list of Groups (first entry is the Group one wants to compare, second entry is the one to stratify by)
# Output: PDF with volcano plots
plot_to_pdf_limma_test_results <- function(resultsdir, all_needed_features, add_info, Group_list, 
                                           only_annotated_features = TRUE, adj_p_val = TRUE, Confidence_cutoff = 0.0, min_groups = 1, nr_ann_features = 20, ann_features = NULL) {
  
  make_timecomp_list <- function(time_points) {
    n_time <- length(time_points)
    comps <- lapply(2:n_time, function(i) c(time_points[i], time_points[1]))
    # Optional: add last vs middle if n_time >= 3, ensure it's unique
    if (n_time >= 3) {
      mid <- time_points[ceiling(n_time/2)]
      last <- time_points[n_time]
      if (!any(sapply(comps, function(x) all(x == c(last, mid))))) {
        comps[[length(comps) + 1]] <- c(last, mid)
      }
    }
    comps
  }
  
  # Function to get next higher 'nice' number (divisible by 2 or 5)
  next_nice_up <- function(x) {
    candidate <- ceiling(x)
    while (!(candidate %% 2 == 0 || candidate %% 5 == 0)) {
      candidate <- candidate + 1
    }
    return(candidate)
  }
  
  # Function to get next lower 'nice' number
  next_nice_down <- function(x) {
    candidate <- floor(x)
    while (!(candidate %% 2 == 0 || candidate %% 5 == 0)) {
      candidate <- candidate - 1
    }
    return(candidate)
  }
  
  all_needed_features <- all_needed_features %>% filter(Confidence >= Confidence_cutoff)
  
  
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
  
  nr_timepoints <- length(unique(full_data$Timepoint))
  
  if(!is.null(ann_features)) {
    nr_ann_features <- length(ann_features)
  }
  
  if(only_annotated_features) {
    all_needed_features <- all_needed_features[!is.na(all_needed_features$Annotation), ]
    
    if (Confidence_cutoff == 0.0) {
      if(adj_p_val) {
        pdf(file = paste0(resultsdir, "limma-test-results_only-annotated-features_adj-p-value_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
        limma_results <- paste0(resultsdir, 'limma_all-resulttables/limma-test-results_only-annotated-features_adj-p-value/')
      } else {
        pdf(file = paste0(resultsdir, "limma-test-results_only-annotated-features_unadj-p-value_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
        limma_results <- paste0(resultsdir, 'limma_all-resulttables/limma-test-results_only-annotated-features_unadj-p-value/')
      }
      
    } else {
      if(adj_p_val) {
        pdf(file = paste0(resultsdir, "limma-test-results_only-annotated-features_adj-p-value_confidence-cutoff", Confidence_cutoff,"_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
        limma_results <- paste0(resultsdir, 'limma_all-resulttables/only-annotated-features_adj-p-value_confidence-cutoff', Confidence_cutoff, '/')
      } else {
        pdf(file = paste0(resultsdir, "limma-test-results_only-annotated-features_unadj-p-value_confidence-cutoff", Confidence_cutoff,"_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
        limma_results <- paste0(resultsdir, 'limma_all-resulttables/only-annotated-features_unadj-p-value_confidence-cutoff', Confidence_cutoff, '/')      
      }
    }
    
  } else {
    
    if(adj_p_val) {
      pdf(file = paste0(resultsdir, "limma-test-results_all-features_adj-p-value_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
      limma_results <- paste0(resultsdir, 'limma_all-resulttables/all-features_adj-p-value/')
    } else {
      pdf(file = paste0(resultsdir, "limma-test-results_all-features_unadj-p-value_nr-ann-feature-", nr_ann_features,".pdf"), width = 8 * nr_timepoints, height = 8 * (length(Group_list) %/% 2 * 6 + 3)) 
      limma_results <- paste0(resultsdir, 'limma_all-resulttables/all-features_unadj-p-value/')
    }
  }
  
  if (!dir.exists(limma_results)) { # Create directory if needed
    dir.create(limma_results, recursive = TRUE)
  }
  
  # Determine unique values for each group
  unique_groups <- lapply(Group_list, function(group) unique(full_data[[group]]))
  names(unique_groups) <- Group_list
  
  # Determine the time points
  time_points <- sort(unique(full_data$Timepoint))
  timecomp_list <- make_timecomp_list(time_points)
  
  plot_list <- list()
  fit_list <- list()
  title_list <- list()
  cond1_list <- list()
  cond2_list <- list()
  
  # Compare the different factors in the first entry of Group_list (Intervention) for certain time points
  for (tp in time_points) {
    data_time <- full_data[full_data$Timepoint == tp, ]

    # Define a blocking factor for subjects to account for repeated measures, pair the samples
    blocking_factor <- factor(data_time$Donor)

    if (nrow(data_time) > 0) {
      data_time[[Group_list[[1]]]] <- as.factor(data_time[[Group_list[[1]]]])
      design <- model.matrix(~ 0 + factor(data_time[[Group_list[[1]]]]), data = data_time)
      colnames(design) <- make.names(colnames(design))

      cond1 <- make.names(sort(unique(factor(data_time[[Group_list[[1]]]])))[1])
      cond2 <- make.names(sort(unique(factor(data_time[[Group_list[[1]]]])))[2])
      
      # Remove the "X" from the condition columns (happens if factor is just a number)
      cond1 <- gsub("^X", "", cond1)
      cond2 <- gsub("^X", "", cond2)
      
      condition_1 <- colnames(design)[grep(cond1, colnames(design))]
      condition_2 <- colnames(design)[grep(cond2, colnames(design))]
      
      assign("condition_1", condition_1, envir = .GlobalEnv) # assigin to global variable, as makeContrasts uses own environment and doesnt find it otherise
      assign("condition_2", condition_2, envir = .GlobalEnv)
      contrast <- makeContrasts(Comparison = paste0(condition_1, " - ", condition_2), levels = design)
      
      nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(data_time)))

      # Prepare the data for analysis
      expression_matrix <- t(data_time[, -c(1:nonnumeric_cols)])

      # Estimate the correlation within subjects using duplicateCorrelation
      corfit <- duplicateCorrelation(expression_matrix, design, block=blocking_factor)

      if (!is.na(corfit$consensus)) {
        fit <- lmFit(expression_matrix, design, block=blocking_factor, correlation=corfit$consensus)
        plot_title <- paste("Intervention Comparison at Time Point ", tp, " (n = ", length(unique(data_time$Donor)), ", paired)")
      } else {
        fit <- lmFit(expression_matrix, design)
        plot_title <- paste("Intervention Comparison at Time Point ", tp, " (n = ", length(unique(data_time$Donor)), ")")
      }
      fit <- eBayes(contrasts.fit(fit, contrast))

      fit_list[[length(fit_list) + 1]] <- fit  # save the fits, need this to have the same axis scale for every plot
      title_list[[length(title_list) + 1]] <- plot_title
      cond1_list[[length(cond1_list) + 1]] <- cond1
      cond2_list[[length(cond2_list) + 1]] <- cond2
    }
  }
  
    logFC_vals <- c()
    adjP_vals <- c()
    
    for (fit in fit_list) {
      tt <- topTable(fit, adjust.method="BH", number=Inf)
      logFC_vals <- c(logFC_vals, tt$logFC)
      adjP_vals <- c(adjP_vals, tt$adj.P.Val)
    }
  
    min_logFC <- next_nice_down(min(logFC_vals, na.rm=TRUE)) # get next "nice" number of the minimal value for axis
    max_logFC <- next_nice_up(max(logFC_vals, na.rm=TRUE))   # get next "nice" number of the maximal value for axis
    max_neglog10_adjP <- next_nice_up(max(-log10(adjP_vals), na.rm=TRUE)) # get next "nice" number of the maximal value for axis
  
    for (i in seq_along(fit_list)) {
      fit <- fit_list[[i]]
      plot_title <- title_list[[i]]
      cond1 <- cond1_list[[i]]
      cond2 <- cond2_list[[i]]
      # If your function doesn't accept axis limits, you must add them
      plot <- enhanced_volcano_plot(
        fit, all_needed_features, plot_title, 
        paste(cond1, "-", cond2), adj_p_val, limma_results,
        xrange = c(min_logFC, max_logFC),
        yrange = c(0, max_neglog10_adjP),
        nr_ann_features = nr_ann_features,
        ann_features = ann_features
      )
      plot_list[[length(plot_list) + 1]] <- plot
    }
    
  
  fit_list <- list()
  title_list <- list()
  cond1_list <- list()
  cond2_list <- list()
  
  # Compare the same factors for certain time points but stratified by the second entry in Group_list (e.g., Sex)
  for (strat_level in unique_groups[[Group_list[[2]]]]) {
    stratified_data <- full_data[full_data[[Group_list[[2]]]] == strat_level, ]

    for (tp in time_points) {
      data_time <- stratified_data[stratified_data$Timepoint == tp, ]

      if (nrow(data_time) > 0) {
        # Define a blocking factor for subjects to account for repeated measures
        blocking_factor <- factor(data_time$Donor)
        
        design <- model.matrix(~ 0 + factor(data_time[[Group_list[[1]]]]), data = data_time)
        colnames(design) <- make.names(colnames(design))
        
        cond1 <- make.names(sort(unique(factor(data_time[[Group_list[[1]]]])))[1])
        cond2 <- make.names(sort(unique(factor(data_time[[Group_list[[1]]]])))[2])
        
        # Remove the "X" from the condition columns (happens if factor is just a number)
        cond1 <- gsub("^X", "", cond1)
        cond2 <- gsub("^X", "", cond2)
        
        condition_1 <- colnames(design)[grep(cond1, colnames(design))]
        condition_2 <- colnames(design)[grep(cond2, colnames(design))]
        
        assign("condition_1", condition_1, envir = .GlobalEnv) # assigin to global variable, as makeContrasts uses own environment and doesnt find it otherise
        assign("condition_2", condition_2, envir = .GlobalEnv)
        contrast <- makeContrasts(Comparison = paste0(condition_1, " - ", condition_2), levels = design)
        
        nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(data_time)))

        # Prepare the data for analysis
        expression_matrix <- t(data_time[, -c(1:nonnumeric_cols)])

        # Estimate the correlation within subjects using duplicateCorrelation
        corfit <- duplicateCorrelation(expression_matrix, design, block=blocking_factor)

        if (!is.na(corfit$consensus)) {
          fit <- lmFit(expression_matrix, design, block=blocking_factor, correlation=corfit$consensus)
          plot_title <- paste0("Intervention Comparison at Time Point ", tp, "\nStratified by ", strat_level, " ", Group_list[[2]], " (n = ", length(unique(data_time$Donor)), ", paired)")
        } else {
          fit <- lmFit(expression_matrix, design)
          plot_title <- paste0("Intervention Comparison at Time Point ", tp, "\nStratified by ", strat_level, " ", Group_list[[2]], " (n = ", length(unique(data_time$Donor)), ")")
        }
        fit <- eBayes(contrasts.fit(fit, contrast))

        fit_list[[length(fit_list) + 1]] <- fit  # save the fits, need this to have the same axis scale for every plot
        title_list[[length(title_list) + 1]] <- plot_title
        cond1_list[[length(cond1_list) + 1]] <- cond1
        cond2_list[[length(cond2_list) + 1]] <- cond2
      }
    }
  }
  
  logFC_vals <- c()
  adjP_vals <- c()
  
  for (fit in fit_list) {
    tt <- topTable(fit, adjust.method="BH", number=Inf)
    logFC_vals <- c(logFC_vals, tt$logFC)
    adjP_vals <- c(adjP_vals, tt$adj.P.Val)
  }
  
  min_logFC <- next_nice_down(min(logFC_vals, na.rm=TRUE)) # get next "nice" number of the minimal value for axis
  max_logFC <- next_nice_up(max(logFC_vals, na.rm=TRUE))   # get next "nice" number of the maximal value for axis
  max_neglog10_adjP <- next_nice_up(max(-log10(adjP_vals), na.rm=TRUE)) # get next "nice" number of the maximal value for axis
  
  for (i in seq_along(fit_list)) {
    fit <- fit_list[[i]]
    plot_title <- title_list[[i]]
    cond1 <- cond1_list[[i]]
    cond2 <- cond2_list[[i]]
    # If your function doesn't accept axis limits, you must add them
    plot <- enhanced_volcano_plot(
      fit, all_needed_features, plot_title, 
      paste(cond1, "-", cond2), adj_p_val, limma_results,
      xrange = c(min_logFC, max_logFC),
      yrange = c(0, max_neglog10_adjP),
      nr_ann_features = nr_ann_features,
      ann_features = ann_features
    )
    plot_list[[length(plot_list) + 1]] <- plot
  }
  

  fit_list <- list()
  title_list <- list()
  cond1_list <- list()
  cond2_list <- list()
  
  # Compare each factor in Intervention between time points (all data)
  for (level in unique_groups[[Group_list[[1]]]]) {
    subset_data_ <- full_data[full_data[[Group_list[[1]]]] == level, ]

    for (time_pair in timecomp_list) {
      time_2 <- time_pair[1] 
      time_1 <- time_pair[2]

      data_t1 <- subset_data_[subset_data_$Timepoint == time_1, ]
      data_t2 <- subset_data_[subset_data_$Timepoint == time_2, ]
      
      subset_data <- rbind(data_t1, data_t2)

      if (nrow(data_t1) > 0 && nrow(data_t2) > 0) {
        # Define a blocking factor for subjects to account for repeated measures
        blocking_factor <- factor(subset_data$Donor)

        subset_data$Timepoint <- as.factor(subset_data$Timepoint)
        design <- model.matrix(~ 0 + subset_data$Timepoint, data = subset_data)
        colnames(design) <- make.names(colnames(design))

        time_point_1_ <- colnames(design)[grep(time_1, colnames(design))]
        time_point_2_ <- colnames(design)[grep(time_2, colnames(design))]
        
        assign("time_point_1_", time_point_1_, envir = .GlobalEnv) # assigin to global variable, as makeContrasts uses own environment and doesnt find it otherise
        assign("time_point_2_", time_point_2_, envir = .GlobalEnv)
        contrast <- makeContrasts(Comparison = paste0(time_point_2_, " - ", time_point_1_), levels = design)

        nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(subset_data)))

        # Prepare the data for analysis
        expression_matrix <- t(subset_data[, -c(1:nonnumeric_cols)])

        # Estimate the correlation within subjects using duplicateCorrelation
        corfit <- duplicateCorrelation(expression_matrix, design, block=blocking_factor)

        if (!is.na(corfit$consensus)) {
          fit <- lmFit(expression_matrix, design, block=blocking_factor, correlation=corfit$consensus)
          plot_title <- paste0(level, " ", Group_list[[1]], " Time Comparison ", "(n = ", length(unique(subset_data$Donor)), ", paired)")
        } else {
          fit <- lmFit(expression_matrix, design)
          plot_title <- paste0(level, " ", Group_list[[1]], " Time Comparison ", "(n = ", length(unique(subset_data$Donor)), ")")
        }
        fit <- eBayes(contrasts.fit(fit, contrast))

        fit_list[[length(fit_list) + 1]] <- fit  # save the fits, need this to have the same axis scale for every plot
        title_list[[length(title_list) + 1]] <- plot_title
        cond1_list[[length(cond1_list) + 1]] <- time_1
        cond2_list[[length(cond2_list) + 1]] <- time_2
      }
    }
  }

  
  logFC_vals <- c()
  adjP_vals <- c()
  
  for (fit in fit_list) {
    tt <- topTable(fit, adjust.method="BH", number=Inf)
    logFC_vals <- c(logFC_vals, tt$logFC)
    adjP_vals <- c(adjP_vals, tt$adj.P.Val)
  }
  
  min_logFC <- next_nice_down(min(logFC_vals, na.rm=TRUE)) # get next "nice" number of the minimal value for axis
  max_logFC <- next_nice_up(max(logFC_vals, na.rm=TRUE))   # get next "nice" number of the maximal value for axis
  max_neglog10_adjP <- next_nice_up(max(-log10(adjP_vals), na.rm=TRUE)) # get next "nice" number of the maximal value for axis
  
  for (i in seq_along(fit_list)) {
    fit <- fit_list[[i]]
    plot_title <- title_list[[i]]
    time1 <- cond1_list[[i]]
    time2 <- cond2_list[[i]]
    # If your function doesn't accept axis limits, you must add them
    plot <- enhanced_volcano_plot(
      fit, all_needed_features, plot_title, 
      paste("Time point", time2, "vs", time1), adj_p_val, limma_results,
      xrange = c(min_logFC, max_logFC),
      yrange = c(0, max_neglog10_adjP),
      nr_ann_features = nr_ann_features,
      ann_features = ann_features
    )
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  
  fit_list <- list()
  title_list <- list()
  cond1_list <- list()
  cond2_list <- list()
  
  # Compare each factor in Intervention between time points (stratified)
  for (strat_level in unique_groups[[Group_list[[2]]]]) {
    stratified_data <- full_data[full_data[[Group_list[[2]]]] == strat_level, ]

    for (level in unique_groups[[Group_list[[1]]]]) {
      subset_data_ <- stratified_data[stratified_data[[Group_list[[1]]]] == level, ]

      for (time_pair in timecomp_list) {
        time_2 <- time_pair[1] 
        time_1 <- time_pair[2]

        data_t1 <- subset_data_[subset_data_$Timepoint == time_1, ]
        data_t2 <- subset_data_[subset_data_$Timepoint == time_2, ]
        
        subset_data <- rbind(data_t1, data_t2)

        if (nrow(data_t1) > 0 && nrow(data_t2) > 0) {
          # Define a blocking factor for subjects to account for repeated measures
          blocking_factor <- factor(subset_data$Donor)
          
          subset_data$Timepoint <- factor(subset_data$Timepoint)
          design <- model.matrix(~ 0 + subset_data$Timepoint, data = subset_data)
          colnames(design) <- make.names(colnames(design))

          time_point_1_ <- colnames(design)[grep(time_1, colnames(design))]
          time_point_2_ <- colnames(design)[grep(time_2, colnames(design))]
          
          assign("time_point_1_", time_point_1_, envir = .GlobalEnv) # assigin to global variable, as makeContrasts uses own environment and doesnt find it otherise
          assign("time_point_2_", time_point_2_, envir = .GlobalEnv)
          contrast <- makeContrasts(Comparison = paste0(time_point_2_, " - ", time_point_1_), levels = design)
          nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(subset_data)))

          # Prepare the data for analysis
          expression_matrix <- t(subset_data[, -c(1:nonnumeric_cols)])

          # Estimate the correlation within subjects using duplicateCorrelation
          corfit <- duplicateCorrelation(expression_matrix, design, block=blocking_factor)

          if (!is.na(corfit$consensus)) {
            fit <- lmFit(expression_matrix, design, block=blocking_factor, correlation=corfit$consensus)
            plot_title <- paste0(level, " ", Group_list[[1]], " Time Comparison\nStratified by ", strat_level, " ", Group_list[[2]], " (n = ", length(unique(subset_data$Donor)), ", paired)")
          } else {
            fit <- lmFit(expression_matrix, design)
            plot_title <- paste0(level, " ", Group_list[[1]], " Time Comparison\nStratified by ", strat_level, " ", Group_list[[2]], " (n = ", length(unique(subset_data$Donor)), ")")
          }
          fit <- eBayes(contrasts.fit(fit, contrast))

          fit_list[[length(fit_list) + 1]] <- fit  # save the fits, need this to have the same axis scale for every plot
          title_list[[length(title_list) + 1]] <- plot_title
          cond1_list[[length(cond1_list) + 1]] <- time_1
          cond2_list[[length(cond2_list) + 1]] <- time_2
        }
      }
    }
  }
  
  
  logFC_vals <- c()
  adjP_vals <- c()
  
  for (fit in fit_list) {
    tt <- topTable(fit, adjust.method="BH", number=Inf)
    logFC_vals <- c(logFC_vals, tt$logFC)
    adjP_vals <- c(adjP_vals, tt$adj.P.Val)
  }
  
  min_logFC <- next_nice_down(min(logFC_vals, na.rm=TRUE)) # get next "nice" number of the minimal value for axis
  max_logFC <- next_nice_up(max(logFC_vals, na.rm=TRUE))   # get next "nice" number of the maximal value for axis
  max_neglog10_adjP <- next_nice_up(max(-log10(adjP_vals), na.rm=TRUE)) # get next "nice" number of the maximal value for axis
  
  for (i in seq_along(fit_list)) {
    fit <- fit_list[[i]]
    plot_title <- title_list[[i]]
    time1 <- cond1_list[[i]]
    time2 <- cond2_list[[i]]
    # If your function doesn't accept axis limits, you must add them
    plot <- enhanced_volcano_plot(
      fit, all_needed_features, plot_title, 
      paste("Time point", time2, "vs", time1), adj_p_val, limma_results,
      xrange = c(min_logFC, max_logFC),
      yrange = c(0, max_neglog10_adjP),
      nr_ann_features = nr_ann_features,
      ann_features = ann_features
    )
    plot_list[[length(plot_list) + 1]] <- plot
  }
  
  
  do.call(grid.arrange, c(plot_list, ncol = nr_timepoints))
  
  # Close PDF device
  dev.off()
}




# Function: Compare time points using limma, calculate changes, and perform clustering
# Input: resultsdir (directory for result files), all_needed_features, add_info, 
#        list of Groups (first entry is the Group one wants to compare, second entry is the one to stratify by)
# Output: Clustering results and plots
timepoint_comparison_and_clustering <- function(resultsdir, all_needed_features, add_info, Group_list, only_annotated_features = TRUE, adj_p_val = TRUE, num_clusters = 3) {
  
  if(only_annotated_features) {
    all_needed_features <- all_needed_features[!is.na(all_needed_features$Annotation), ]
  }
  
  # Reshape all_needed_features to have samples as rows and features as columns
  if ("rtime_group" %in% colnames(all_needed_features)) {
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
  
  reshaped_data <- reshape2::dcast(reshaped_data, variable ~ id, value.var = "value")
  colnames(reshaped_data)[1] <- "Sample"
  rownames(reshaped_data) <- reshaped_data[,1]
  pure_data <- reshaped_data[,-1]
  
  # Merge reshaped data with add_info
  merged_data <- merge(reshaped_data, add_info, by.x = "Sample", by.y = "Sample")
  
  merged_data[, 2:(length(merged_data) - (ncol(add_info) -1))] <- as.data.frame(apply(merged_data[, 2:(length(merged_data) - 6)], 2, replace_na_and_Inf))
  pure_data <- as.data.frame(apply(pure_data, 2, replace_na_and_Inf))
  pure_data <- log2(pure_data)
  row.names(pure_data) <- merged_data$Sample 
  
  # Merge pure_data with metadata
  full_data <- cbind(merged_data[, c("Sample", "Donor", "Timepoint", unlist(Group_list))], pure_data)
  
  # Determine unique values for each group
  unique_groups <- lapply(Group_list, function(group) unique(full_data[[group]]))
  names(unique_groups) <- Group_list
  
  # Determine the time points
  time_points <- sort(unique(full_data$Timepoint))
  
  # Matrix to store logFC values
  logFC_matrix <- matrix(NA, nrow = ncol(pure_data), ncol = length(time_points) - 1)
  rownames(logFC_matrix) <- colnames(pure_data)
  colnames(logFC_matrix) <- paste0("LogFC_", time_points[-1], "_vs_", time_points[-length(time_points)])
  
  for (i in 2:length(time_points)) {
    data_t1 <- full_data[full_data$Timepoint == time_points[i-1], ]
    data_t2 <- full_data[full_data$Timepoint == time_points[i], ]
    
    # Ensure that both time points have the same samples
    common_donors <- intersect(data_t1$Donor, data_t2$Donor)
    data_t1 <- data_t1[data_t1$Donor %in% common_donors, ]
    data_t2 <- data_t2[data_t2$Donor %in% common_donors, ]
    
    if (nrow(data_t1) > 0 && nrow(data_t2) > 0) {
      combined_data <- rbind(data_t1, data_t2)
      combined_data$Timepoint <- factor(combined_data$Timepoint)
      
      design <- model.matrix(~ 0 + Timepoint, combined_data)
      
      # Dynamically set up the contrast using colnames(design)
      time_point_1 <- colnames(design)[grep(time_points[i-1], colnames(design))]
      time_point_2 <- colnames(design)[grep(time_points[i], colnames(design))]
      print(paste0(time_point_2, " - ", time_point_1))
      contrast <- makeContrasts(Comparison = paste0(time_point_2, " - ", time_point_1), levels = design)
      
      nonnumeric_cols <- sum(!grepl("^[0-9]+$", colnames(data_t1)))
      
      expression_matrix <- t(combined_data[, -(1:nonnumeric_cols)])
      fit <- lmFit(expression_matrix, design)
      fit <- eBayes(contrasts.fit(fit, contrast))
      
      top_table <- topTable(fit, adjust.method = if(adj_p_val) "BH" else "none", number = Inf)
      
      # Store logFC values in the matrix
      logFC_matrix[, i-1] <- top_table$logFC[match(rownames(logFC_matrix), rownames(top_table))]
    }
  }
  
  # Perform clustering analysis on logFC matrix
  distance_matrix <- dist(logFC_matrix, method = "euclidean")
  clustering <- hclust(distance_matrix, method = "ward.D2")
  
  # Cut the tree into clusters
  cluster_assignments <- cutree(clustering, k = num_clusters)
  
  # Plot the dendrogram
  pdf(file = paste0(resultsdir, "clustering_dendrogram_limma.pdf"), width = 15, height = 10)
  plot(clustering, labels = FALSE, main = "Clustering Dendrogram (Based on LogFC)", sub = "", xlab = "", ylab = "Height")
  rect.hclust(clustering, k = num_clusters, border = 2:6)
  dev.off()
  
  feature_annotations <- all_needed_features$Annotation[match(rownames(logFC_matrix), all_needed_features$id)]
  # Combine cluster assignments with annotations
  clustered_features <- data.frame(
    id = rownames(logFC_matrix),
    Annotation = feature_annotations,
    Cluster = cluster_assignments
  )
  
  logFC_matrix <- as.data.frame(logFC_matrix)
  logFC_matrix$id <- rownames(logFC_matrix)
  
  clustered_feature_logFC <- merge(clustered_features, logFC_matrix, by = "id")
  clustered_feature_logFC <- clustered_feature_logFC[order(clustered_feature_logFC$Cluster), ]
  
  # Return the clustering results
  return(list(clustered_features = clustered_features, logFC_matrix = logFC_matrix))
}




