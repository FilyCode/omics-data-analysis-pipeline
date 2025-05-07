source("//prot-data/MS Storage 3/tmp_PT/R Project/raw-to-mzML-file_conversion.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/run-mzMine-batchmode.R")
source("//prot-data/MS Storage 3/tmp_PT/R Project/annotate-features-with-sirius.R")

# Get directory where data is located
directory <- input_directory()

# Convert raw files to mzML files and save to according directory
convert_raw_to_mzML_file(directory)

# Run mzMine Batchmode, extract features from mzML files and save data to according directory
run_mzMine_batchmode(directory)

# Annotate features with SIRIUS and save data to according directory
annotate_features_with_sirius(directory)








# can be deleted later, just needed for checking data

area_automatic1 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_CLI1.csv")
area_automatic2 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_CLI2.csv")
area_automatic3 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_CLI3.csv")
area_manual1 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_autom-batchfile-done-manually.csv")
area_manual2 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_autom-batchfile-done-manually2.csv")
area_manual3 <- read.csv("W:/Exploris480/2024/test_directory/data/areas_autom-batchfile-done-manually3.csv")


identical(area_manual1, area_automatic1)
identical(area_manual2, area_automatic1)
identical(area_manual3, area_automatic1)
identical(area_manual1, area_automatic2)
identical(area_manual2, area_automatic2)
identical(area_manual3, area_automatic2)
identical(area_automatic1, area_automatic2)
identical(area_automatic1, area_automatic3)
identical(area_automatic2, area_automatic3)
identical(area_manual1, area_manual2)
identical(area_manual1, area_manual3)
identical(area_manual2, area_manual3)


# Function to identify differing columns
find_differences <- function(df1, df2) {
  differing_cols <- character()  # Vector to store differing column names
  
  for (col in names(df1)) {
    if (!identical(df1[[col]], df2[[col]])) {
      differing_cols <- c(differing_cols, col)  # Add differing column name
    }
  }
  
  return(differing_cols)
}

# Get the differing columns
differing_columns <- find_differences(area_automatic2, area_automatic3)
differing_columns

length(area_automatic2$mz_range.min[area_automatic3$mz_range.min != area_automatic2$mz_range.min])
length(area_automatic2$mz_range.max[area_automatic3$mz_range.max != area_automatic2$mz_range.max])
length(area_manual2$mz_range.min[area_manual3$mz_range.min != area_manual2$mz_range.min])
length(area_manual2$mz_range.max[area_manual3$mz_range.max != area_manual2$mz_range.max])
length(area_automatic2$mz_range.min[area_automatic2$mz_range.min != area_manual2$mz_range.min])
length(area_automatic2$mz_range.max[area_automatic2$mz_range.max != area_manual2$mz_range.max])

identical(area_manual2$mz, area_manual3$mz)
identical(area_manual2$mz, area_automatic2$mz)
identical(area_automatic2$mz, area_automatic3$mz)
