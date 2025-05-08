# omics-data-analysis-pipeline
This code was used to analyze he NuriNeuro study data including finger sweat and blood plasma LC-MS metabolomics data but can be adapted to any omics dataset.

The 'fingersweat_analysis_script.R' file calls all the necessary functions to run the data analysis, statistical tests and visualization.

The called functions are in the files 'feature-visualization-script.R' and 'MS-data-analysis_functions.R'.

First the .raw files need to prepared and summarized into a matrix (Samples as columns, features as rows).
This has been done by converting them into .mzML files, running a mzMine Batch mode, with a subsequent SIRIUS annotation. All this can be either done by hand or by running the 'get_untargeted_annotated_features_from_raw_files.R' script.
The resulting .csv file (+SIRIUS summary file) is used to run the data analysis.
Metadata are also need, which can be create with 'create_metadata-file_NutriNeuro.R'.
