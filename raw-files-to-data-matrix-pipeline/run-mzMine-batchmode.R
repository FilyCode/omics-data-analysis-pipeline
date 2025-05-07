# File to perform the feature detection with mzMine and save the results

run_mzMine_batchmode <- function(directory) {
  mzmine_path <- "C:/Program Files/mzmine/mzmine.exe"
  #user_file <- "C:/Users/student/.mzmine/users/trollmann.mzuser"
  batch_template <- "//prot-data2/Group/Studenten_Schueler/Philipp_Trollmann/mzMine_Batch-methods/full_method_for_automatisation.mzbatch"
  input_dir <- paste0(directory, "/mzML/")
  output_dir <- paste0(directory, "/data")

  # Check if directory with input mzML files exists, if not stop calculation
  if (dir.exists(input_dir)) {
    all_files <- list.files(input_dir)
    mzML_files <- all_files[grepl("\\.mzML$", all_files, ignore.case = TRUE)]
    
    if (length(mzML_files) == 0) {
      return(cat("\nThe directory doesn't contain any mzML files! Please check it and then start again."))
    }
    
    if (length(mzML_files) < 6) {
      return(cat("\nThe directory contains to few files. You need at least 1 file for the experimental data and ProcBlank of each measuring mode (hpos, npos, neg)! So at least 6 files. Please check it and then start again."))
    }
  } else {
    return(cat("\nThere is no mzML directory! Please check it and then start again."))
  }
    
  # Check if directory for output files already exists and already has area.csv and feature-data.mgf files in them
  # create output directory if it doesnt exist
  if (dir.exists(output_dir)) {
    csv_file <- file.exists(file.path(output_dir, "areas.csv"))
    mgf_file <- file.exists(file.path(output_dir, "feature-data.mgf"))
    
    if (csv_file & mgf_file) {
      cat("\nThe directory already contains the result files of the mzMine Batchmode.\nDo you want to skip the mzMine Batchmode?")
      input <- readline(prompt = "(y/n): ")
      
      if (any(input == c("y","Y", "yes", "YES", "j", "J", "ja", "JA"))) {
        return(cat("\nmzMine Batchmode got skipped!\n"))
      }
    }
    
  } else {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Ask if the MS/MS data shall be merged
  cat("\nDo you want to merge the MS/MS data when writing the data to the 'feature-data.mgf' file?\nThe feature annotation by SIRIUS takes longer for unmerged data.\nKeeping it unmerged is especially recommended for not so well known samples.")
  input_merge <- readline(prompt = "Merge MS/MS data? (y/n): ")
  merge <- "false"
  
  if (any(input_merge == c("y","Y", "yes", "YES", "j", "J", "ja", "JA"))) {
    merge <- "true"
  }
  
  
  # Construct input files section
  input_files <- list.files(input_dir, pattern = "\\.mzML$", full.names = TRUE)
  input_section <- paste(
    paste(sapply(input_files, function(file) paste0("<file>", file, "</file>")), collapse = "\n"),
    sep = "\n"
  )
  
  # Construct output files section
  mgf_output <- paste0(
    "<current_file>", paste0(output_dir, "\\feature-data.mgf"), "</current_file>"
  )
  
  csv_output <- paste0(
    "<current_file>", paste0(output_dir, "\\areas.csv"), "</current_file>"
  )
  
  ms_ms_merge <- paste0(
    '<parameter name="Merge MS/MS" selected="', merge, '">'
  )
  
  # Read the batch template content
  batch_template_content <- readLines(batch_template, warn = FALSE)
  
  # Find where to insert input_section, mgf_output, and csv_output in batch_template_content
  idx_input <- grep("<insert Filepaths here>", batch_template_content)
  idx_mgf <- grep("feature-data.mgf</current_file>", batch_template_content)
  idx_csv <- grep("areas.csv</current_file>", batch_template_content)
  idx_merge <- grep('<parameter name="Merge MS/MS" selected', batch_template_content)
  
  # Replace the placeholder sections with actual content
  batch_template_content[idx_input] <- input_section
  batch_template_content[idx_mgf] <- mgf_output
  batch_template_content[idx_csv] <- csv_output
  batch_template_content[idx_merge] <- ms_ms_merge
  
  # Write modified content to a temporary .mzbatch file
  temp_batch_file <- tempfile(fileext = ".mzbatch")
  writeLines(batch_template_content, temp_batch_file)
  
  # Creates command for system call to command line (CLI)
  command <- paste(
    shQuote(mzmine_path),
    "-batch", shQuote(temp_batch_file)
  )
  
  # Print the command to verify
  #cat("\nCommand to be executed:\n", command, "\n")
  
  # Execute the command
  system(command)
  
  # Clean up
  file.remove(temp_batch_file)
  #cat("\nTemporary files deleted.")
  
  cat("\nThe mzMine Batchmode has finished!\n")
}
