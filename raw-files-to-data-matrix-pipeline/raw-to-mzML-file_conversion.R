# File to perform the raw to mzML file conversion automatically

library(stringr)
library(parallel)
library(IDBacApp)
library(fs)

input_directory <- function() {
  check_existence <- FALSE
  
  cat("Enter the Path of directory where your raw files are located (use '/' to separate folders).\nExample: //prot-data/MS Storage 3/Exploris480/MyData ")
  
  while(!check_existence) {
    directory <- readline(prompt = "Directory: ")
    
    if (!dir.exists(directory)) {
      cat("This directory doesn't exist. Please input an existing one.")
    } else {
      check_existence <- TRUE
    }
  }
  
  return(directory)
}

convert_raw_to_mzML_file <- function(directory) {
  # Define variables for the output directory and the path of the software 
  output_dir <- paste0(directory, "/mzML")
  software_dir <- "C:/Users/proteomics/AppData/Local/Apps/ProteoWizard 3.0.24123.38027b6 64-bit/"
  msconvert_path <- paste0(software_dir, "msconvert.exe") # Define the msconvert executable path
  datadir <- paste0(directory, "/data/SIRIUS")
  resultsdir <- paste0(directory, "/results")
  
  # Check if directory for mzML files already exists and already has mzML files in them
  # create output directory if it doesnt exist
  if (dir.exists(output_dir)) {
    all_files_ <- list.files(output_dir)
    mzML_files <- all_files_[grepl("\\.mzML$", all_files_, ignore.case = TRUE)]

    partly_conversion <- FALSE
    
    if (length(mzML_files) > 0) {
      cat("\nThe directory already contains .mzML files.\nDo you want to skip the conversion of raw files to mzML files completly?")
      input <- readline(prompt = "(y/n): ")
      
      if (any(input == c("y","Y", "yes", "YES", "j", "J", "ja", "JA"))) {
        return(cat("\nFile conversion got skipped!\n"))
      }
      
      cat("\nThe directory already contains some .mzML files.\nDo you want to skip the conversion of of already converted raw files to mzML files?")
      input <- readline(prompt = "(y/n): ")
      
      if (any(input == c("y","Y", "yes", "YES", "j", "J", "ja", "JA"))) {
        partly_conversion <- TRUE
      }
    }
    
  } else {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (!dir.exists(datadir)) {
    dir.create(datadir, recursive = TRUE)
  }
  
  if (!dir.exists(resultsdir)) {
    dir.create(resultsdir, recursive = TRUE)
  }
  
  cat("\nThe mzML file conversion has started! This can take up to several minutes.\n")
  
  # List all .raw files in the import directory
  # all_files <- list.files(directory, pattern = "\\.raw$", full.names = TRUE, recursive = TRUE)
  all_files <- as.character(dir_ls(directory, recurse = TRUE, glob = "*.raw")) # should be faster
  
  # Remove already converted files if wanted so
  if (partly_conversion) {
    all_files <- all_files[!sub("\\.raw$", ".mzML", basename(all_files)) %in% mzML_files]
  }
  
    # Filter files to exclude those with "Blank" but include those with "ProcBlank"
  filtered_files <- all_files[(!str_detect(all_files, "Blank") | str_detect(all_files, "ProcBlank")) & !str_detect(all_files, "Standby")]
  
  
  # Define the command line options
  options <- '--zlib --filter "peakPicking vendor msLevel=1-" --filter "titleMaker <RunId>.<ScanNumber>.<ScanNumber>.<ChargeState> File:"""^<SourcePath^>""", NativeID:"""^<Id^>""""'
  
  # Construct the full command with all filtered files
  command_for_parallel <- lapply(filtered_files,
                                 function(file){
                                   paste0(shQuote(msconvert_path), ' ', 
                                          shQuote(file), ' ',
                                          options, 
                                          ' -o ', shQuote(output_dir))
                                 })
  
  # Execute the command and run it parallel on cluster
    numCores <- detectCores() - 1
    
    # reduce number of used cores for smaller assignments to save time at cluster creation
    if (length(filtered_files) < numCores) {
      numCores <- length(filtered_files)
    }
    
    cl <- makeCluster(numCores)

    run_MSConverter_parallel <- function(x){
      system(command = x,
             wait = TRUE)
    }
    
    parLapply(cl,
              command_for_parallel,
              run_MSConverter_parallel)
    
    stopCluster(cl)
  
  cat("\nThe mzML file conversion has finished!\n")
}
