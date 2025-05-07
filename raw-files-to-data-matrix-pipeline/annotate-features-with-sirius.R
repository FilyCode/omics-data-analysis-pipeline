# File to perform the feature annotation with SIRIUS automatically

library(future.apply)

annotate_features_with_sirius <- function(directory) {
  # Define the function to run SIRIUS login
  sirius_login <- function(email, password, sirius_path) {
    # Set environment variables for username and password
    Sys.setenv(SIRIUS_USER = email, SIRIUS_PASSWORD = password)
    
    # Construct the SIRIUS login command using environment variables
    sirius_login_cmd <- sprintf('"%s" login --user-env=SIRIUS_USER --password-env=SIRIUS_PASSWORD', sirius_path)
    
    # Execute the SIRIUS login command
    login_output <- system(sirius_login_cmd, intern = TRUE, wait = TRUE)
    
    # Check for successful login
    if (any(grepl("Login successful!", login_output))) {
      cat("Login successful!")
    } else {
      cat("Login failed!")
      print(login_output)
    }
  }
  
  # Define the function to run SIRIUS annotation
  run_sirius <- function(input_mgf, output_dir, config_options, sirius_path) {
    # Construct the SIRIUS command
    sirius_cmd <- sprintf(
      '%s -i %s -o %s %s write-summaries',
      shQuote(sirius_path),
      shQuote(input_mgf),
      shQuote(output_dir),
      config_options
    )
    
    # sirius_cmd <- sprintf(
    #   '%s -i %s -o %s %s -o %s',
    #   shQuote(sirius_path),
    #   shQuote(input_mgf),
    #   shQuote(paste0(input_dir,"/sirius-data.sirius")),
    #   config_options,
    #   shQuote(output_dir)
    # )
    
    # sirius_cmd <- sprintf(
    #   '%s -i %s -o %s %s',
    #   shQuote(sirius_path),
    #   shQuote(input_mgf),
    #   shQuote(paste0(input_dir,"/sirius-data.sirius")),
    #   config_options
    # )
    
    # sirius_cmd <- sprintf(
    #   '%s --help',
    #   shQuote(sirius_path)
    # )
    
    # Execute the SIRIUS command
    system(sirius_cmd)
  }
  
  # Define the function to calculate the tanimoto similaritie score
  get_tanimoto_score <- function(output_dir, sirius_path) {
    sim_dir <- paste0(directory, "/data/similarity")
    dir.create(sim_dir, recursive = TRUE)
    # Construct the SIRIUS command
    sirius_cmd <- sprintf(
      '%s -i %s similarity --numpy --tanimoto --tanimoto-canopus -d %s',
      shQuote(sirius_path),
      shQuote(output_dir),
      shQuote(sim_dir)
    )
    
    # Execute the SIRIUS command
    system(sirius_cmd)
   }
  
  # Function to delete all files and subdirectories in a directory
  delete_directory_contents <- function(dir_path) {
    cat("\nDeleting old annotation data from: ", dir_path, "\nThis can take several minutes.\n")
    
    # Function to delete file or directory
    delete_file_or_directory <- function(path) {
      if (file.info(path)$isdir) {
        # Recursively delete directories
        unlink(path, recursive = TRUE)
      } else {
        # Delete files
        file.remove(path)
      }
    }
    
    # List all files and directories in the specified directory
    files <- list.files(dir_path, full.names = TRUE)
    
    # Parallel processing setup
    plan(multicore)  # Use multicore backend for parallel processing
    
    # Parallel foreach loop to delete files and directories
    future_lapply(files, function(file) {
      delete_file_or_directory(file)
    })
    
    # Verify if deletion was successful
    if (length(list.files(dir_path, full.names = TRUE)) == 0) {
      cat("\nFiles deleted successfully.\n")
    } else {
      cat("\nError deleting files.\n")
    }
  }
  
  
  
  # Define the path to the SIRIUS executable
  sirius_path <- "C:/Program Files/sirius/sirius.exe"
  # Define the input .mgf file and output directory
  input_dir <- paste0(directory, "/data")
  input_mgf <- paste0(directory, "/data/feature-data.mgf")
  output_dir <- paste0(directory, "/data/SIRIUS")
  
  # Check if directory with input mgf file exists, if not stop calculation
  if (dir.exists(input_dir)) {
    if (!file.exists(input_mgf)) {
      return(cat("\nThe directory doesn't contain the feature-data.mgf file! Please check it and then start again."))
    }
  } else {
    return(cat("\nThere is no data directory! Please check it and then start again."))
  }
    
  # Check if directory for output files already exists and already has SIRUS files in it
  # create output directory if it doesnt exist
  if (dir.exists(output_dir)) {
    if (file.exists(file.path(output_dir, "compound_identifications.tsv"))) {
      cat("\nThe directory already contains the result files of the SIRIUS annotation.\nDo you want to skip the SIRUS annotation?")
      cat("\nIf you want to make a new annotation run with SIRIUS the old files will be deleted (can take some minutes).")
      input <- readline(prompt = "(y/n): ")
      
      if (any(input == c("y", "Y", "yes", "YES", "j", "J", "ja", "JA"))) {
        return(cat("\nSIRIUS annotation got skipped!\n"))
      }
      
      # if it already exists and we calculate it new, than it is needed to delete the old data
      delete_directory_contents(output_dir)
      
    } else if (length(list.files(output_dir, full.names = TRUE)) != 0) {
      cat("There are files in the output directory without the correct summary file. Files will be deleted before the new annotation.\n")
      delete_directory_contents(output_dir)    
    }
    
  } else {
    dir.create(output_dir, recursive = TRUE)
  }
  
  
  cat("\n\nFeature Annotation with SIRIUS is starting! This will probably take around 30 minutes (but can take up to several hours).\n\n\n")
  
  # Prompt the user for email and password
  # email <- readline(prompt = "Enter your SIRIUS email: ")
  # password <- readline(prompt = "Enter your SIRIUS password: ")
  
  # Later use a account only created for this
  email <- "a11915661@unet.univie.ac.at"
  password <- "#16Krashka26"
  
  # Login to SIRIUS
  sirius_login(email, password, sirius_path)
  
  # Define the configuration options
  config_options <- "config --IsotopeSettings.filter=true --FormulaSearchDB= --Timeout.secondsPerTree=0 --FormulaSettings.enforced=HCNOP --Timeout.secondsPerInstance=0 --AdductSettings.detectable=[[M-H2O+H]+,[M+K]+,[M-H]-,[M+Cl]-,[M+Na]+,[M+H3N+H]+,[M+H]+,[M+Br]-,[M-H2O-H]-,[M-H4O2+H]+] --UseHeuristic.mzToUseHeuristicOnly=650 --AlgorithmProfile=orbitrap --IsotopeMs2Settings=IGNORE --MS2MassDeviation.allowedMassDeviation=5.0ppm --NumberOfCandidatesPerIon=1 --UseHeuristic.mzToUseHeuristic=300 --FormulaSettings.detectable=B,Cl,Br,Se,S --NumberOfCandidates=10 --AdductSettings.enforced=, --AdductSettings.fallback=[[M+K]+,[M+Cl]-,[M-H]-,[M+Na]+,[M+H]+,[M+Br]-] --FormulaResultThreshold=true --InjectElGordoCompounds=true --StructureSearchDB=BIO --RecomputeResults=false formula fingerprint structure canopus"
  #config_options <- "config --AlgorithmProfile=orbitrap formulas fingerprints classes structures denovo-structures"
  
  # Run the SIRIUS function
  run_sirius(input_mgf, output_dir, config_options, sirius_path)
  #get_tanimoto_score(output_dir, sirius_path)
    
  cat("\n\nFeature Annotation with SIRIUS has finished!\n\n")
}
