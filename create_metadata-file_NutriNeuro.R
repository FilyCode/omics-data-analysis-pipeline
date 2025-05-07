source("//prot-data/MS Storage 3/tmp_PT/R Project/MS-data-analysis_functions.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/feature-visualization-script.R")

library(stringr)
library(openxlsx)
library(dplyr)

wd = "//prot-data/MS Storage 3/Exploris480/"
setwd(wd)

exp = "NutriNeuro"
datadir = paste0(wd, exp, "/data/")
filedir = paste0(wd, exp, "/mzML/")
more_meta_data_file <- paste0(wd, exp, "/participant_file.csv")
delivery_file <- paste0(wd, exp, "/Bloodsample-delivery.xlsx")
storage_file <- paste0(wd, exp, "/Sample-Storage-List.xlsx")

exp_data <- process_merged_files(datadir)
all_data <- exp_data$data

# Make Metadata dataframe
metadata_df <- data.frame(Sample = colnames(all_data)) %>% filter(!Sample %in% c('QC', 'ProcBlank'))
metadata_df$Experiment <- exp

# Extract Information from sample names
# Function to handle extracting parts of the sample name
  extract_sample_parts <- function(sample_name) {
      parts <- str_split(sample_name, "_")[[1]]
      
      if ("Plasma" %in% parts) {
        Sample.Type <- 'Plasma'
        donor <- as.integer(parts[3])
        sampling_day <- as.integer(parts[4])
        timepoint <- as.integer(parts[5])
      } else if  ("FiS" %in% parts) {
        Sample.Type <- 'FiS'
        donor <- as.integer(parts[3])
        sampling_day <- as.integer(parts[4])
        timepoint <- as.integer(parts[5])
      } else {
        Sample.Type <- 'FiS'
        donor <- as.integer(parts[2])
        sampling_day <- as.integer(parts[3])
        timepoint <- as.integer(parts[4])
      }
    
    return(c(Sample.Type, donor, sampling_day, timepoint))
  }
  
  # Apply the function to extract parts and create new columns
  extracted_parts <- t(sapply(metadata_df$Sample, extract_sample_parts))
  
  # Assign the extracted parts to new columns
  metadata_df$Sample.Type <- as.character(extracted_parts[, 1])
  metadata_df$Donor <- as.integer(extracted_parts[, 2])
  metadata_df$Sampling.Day <- as.integer(extracted_parts[, 3])
  metadata_df$Timepoint <- as.integer(extracted_parts[, 4])
  metadata_df$Sample.ID <- paste0(metadata_df$Sample.Type, '.', metadata_df$Donor, '.', metadata_df$Sampling.Day, '.', metadata_df$Timepoint)
  metadata_df$Donor.Sampling.Day <- paste0(metadata_df$Donor, '.', metadata_df$Sampling.Day)
  
  # adjust timepoints to make it more intuitive, shift first time point (fasted state) to T0 instead of T1
  metadata_df$Timepoint <- metadata_df$Timepoint - 1
  
  # align plasma sampling time points to finger sweat time points
  metadata_df <- metadata_df %>%
    mutate(Timepoint = ifelse(Sample.Type == "Plasma" & Timepoint != 0, Timepoint + 2, Timepoint))
  
# Read in additional data
delivery_data <- read.xlsx(delivery_file, cols = c(1,2,4))
delivery_data[[1]] <- as.Date(delivery_data[[1]], origin = "1899-12-30")
all_storage_data <- rbind(read.xlsx(storage_file, sheet="FingerSweat", cols = c(3,6,7,8)), read.xlsx(storage_file, sheet="Plasma (EDTA)", cols = c(3,8,9,10)))
all_storage_data$Sample.ID <- paste0(all_storage_data$Sample.Type, '.', all_storage_data$Sample.ID) 
all_storage_data <- all_storage_data[!is.na(all_storage_data$Batch.Nr), ]
all_storage_data$Sample.Type <- NULL
more_meta_data <- read.csv(more_meta_data_file)

# Merge more data into metadata df
metadata_df <- merge(metadata_df, all_storage_data, by = "Sample.ID")
metadata_df <- merge(metadata_df, delivery_data, by.x = 'Donor.Sampling.Day', by.y = 'Sample.(Proband.to.Food)')
metadata_df <- metadata_df %>% mutate(Sample.Temp.at.delivery = as.numeric(str_remove(Sample.Temp.at.delivery, "°C")))

# Extract donor ID from participant_id (removing 'sub-' and converting to numeric)
more_meta_data <- more_meta_data %>% mutate(Donor = as.numeric(str_remove(participant_id, "sub-")))

# Merge gender, age, and BMI based on Donor
metadata_df <- metadata_df %>% 
  left_join(more_meta_data, by = "Donor")

# Function to assign meal type based on Donor and Sampling Day
assign_meal_type <- function(donor, sampling_day, meal_data) {
  # Get the first meal type for the donor
  first_meal <- meal_data %>% filter(Donor == donor) %>% pull(first_meal)
  
  # Assign meal type based on Sampling Day
  if (length(first_meal) > 0) {
    return(ifelse(sampling_day == 1, first_meal, ifelse(first_meal == "high c/p", "low c/p", "high c/p"))) # assign second meal type
  } else {
    return(NA)  # Return NA if donor not found in meal data
  }
}

# Apply function to assign meal type
metadata_df$Intervention <- mapply(assign_meal_type, metadata_df$Donor, metadata_df$Sampling.Day, MoreArgs = list(meal_data = more_meta_data))

write.csv(metadata_df, paste0(wd, exp, '/Metadata.csv'), row.names = FALSE)

