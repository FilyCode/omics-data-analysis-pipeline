# Omics Data Analysis Pipeline

This repository contains an R-based pipeline for the analysis of omics datasets. Originally developed for LC-MS metabolomics data from the NutriNeuro study (including finger sweat and blood plasma), it is designed to be adaptable to various omics data types.

## Pipeline Components

The core data analysis, statistical tests, and visualization are orchestrated by:

*   `fingersweat_analysis_script.R`: The main script that calls all necessary functions.

Supporting functions are organized into:

*   `feature-visualization-script.R`: Contains functions for data visualization.
*   `MS-data-analysis_functions.R`: Houses core functions for MS data analysis.

## Data Preparation

Before running the analysis pipeline, raw omics data files require pre-processing into a feature matrix (Samples as columns, Features as rows) and associated metadata.

### Required Input Files:

1.  **Feature Matrix:** A `.csv` file containing the summarized features (e.g., from LC-MS raw data), accompanied by a SIRIUS summary file for annotation.
2.  **Metadata File:** A `.csv` file containing sample metadata.

### Pre-processing Workflow:

The initial `.raw` files can be processed into the required feature matrix and SIRIUS annotation files through the following steps:

1.  **Raw to mzML Conversion:** Convert `.raw` files to `.mzML` format.
2.  **Feature Detection & Quantification:** Run a mzMine Batch mode for feature detection and quantification.
3.  **SIRIUS Annotation:** Perform SIRIUS annotation on the processed features.

Alternatively, the `get_untargeted_annotated_features_from_raw_files.R` script can automate these steps.

Metadata can be generated using the the `create_metadata-file_NutriNeuro.R` script.

## Authorship
This pipeline and its associated scripts were solely developed by Philipp Trollmann as part of his Master Thesis in Biological Chemistry in the group of Dr. Samuel Meier-Menches and Dr. Christopher Gerner at University of Vienna.
Script used for Master Thesis: "Investigating dietary interventions on the human metabolome using LC/MS for matched plasma and finger sweat samples" (DOI: 10.25365/thesis.78505)

