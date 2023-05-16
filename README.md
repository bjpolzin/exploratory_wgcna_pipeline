# Table of Contents
1. [Exploratory WGCNA Pipeline](#exploratory-wgcna-pipeline)
2. [Prerequisites](#prerequisites)
3. [Pipeline Overview](#pipeline-overview)
4. [Usage](#usage)
5. [Output](#output)
6. [Acknowledgements](#acknowledgements)

# Exploratory WGCNA Pipeline <a name="exploratory-wgcna-pipeline"></a>
This repository contains an R-based exploratory pipeline for Weighted Gene Co-expression Network Analysis (WGCNA), designed for initial exploratory analyses to compare the results of different soft threshold powers. The pipeline takes gene expression and behavioral data as input and generates co-expression networks, soft-threshold plots, module-trait relationship plots, and correlation/p-value tables. The pipeline is modular, allowing users to modify or extend the steps as needed.

This pipeline is built upon the foundational WGCNA codebase, originally developed by Peter Langfelder and Steve Horvath. For citation information, please refer to the acknowledgements section at the bottom of [this page](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/). If you are interested in further exploring the WGCNA method, additional resources can be found on this page as well.

## Prerequisites <a name="prerequisites"></a>
Before running the pipeline, please make sure you have the following R packages installed:

* `tidyverse`
* `dplyr`
* `WGCNA`
* `data.table`

## Pipeline Overview <a name="pipeline-overview"></a>
The main pipeline script is `exploratory_wgcna_pipeline.R`. This script sets up the necessary parameters and data directories and sources the following three modules:

1. `01_read_organize_data.R` - Reads and organizes gene expression and behavioral data.
2. `02_investigate_sft.R` - Creates soft-threshold plots to determine appropriate powers to explore.
3. `03_create_networks_export_results.R` - Performs WGCNA and exports the results.

## Usage <a name="usage"></a>
1. Open `exploratory_wgcna_pipeline.R` in your working directory. Note that multiple subdirectories will be created relative to the home directory. To avoid potential issues, it is recommended to use an R project file for easier management of directories and dependencies.
2. Configure the parameters at the beginning of `exploratory_wgcna_pipeline.R`. Set the static variables, such as the study identifier, data directories, etc.
3. Execute Module 1 (`01_read_organize_data.R`) and Module 2 (`02_investigate_sft.R`).
4. Examine the generated `sft_plot` to determine the powers you want to investigate. Create a vector of all the soft-threshold powers you wish to explore in your results.
5. Run Module 3 (`03_create_networks_export_results.R`).
6. Enjoy reviewing your results! I hope you discover something interesting and insightful!

The following is a concise summary of each module's function. For a comprehensive understanding and in-depth explanation of the code, please refer to the individual R script modules, which contain detailed instructions and descriptions.

### Module 1: Read in data and organize data for input to WGCNA
The first module reads and processes the gene expression and behavioral data. It takes care of transposing the data and filtering out low-expression genes. The module ensures that the sample IDs for both datasets are identical and prints a message if they are consistent or a warning if they do not match.

### Module 2: Create soft-threshold plots to determine powers to explore
The second module creates soft-threshold plots to determine the powers to explore. It performs network topology analysis and plots the scale-free topology fit index and mean connectivity as functions of the soft-thresholding power. The generated plots are saved as PDF files in the `sft_results/plots` directory, and the scale-free topology data is saved as an RData file in the `sft_results/working_rdata` directory.

### Module 3: Perform WGCNA and export results
The third module iterates through the specified soft-threshold powers, creating networks, calculating module eigengenes and correlations, and exporting the results. The module generates dendrogram plots, module-trait relationship plots, correlation and p-value tables, and hub genes for each power. The results are saved in separate directories for each power under the `results` directory.

## Output <a name="output"></a>
The pipeline generates the following output files in automatically generated directories:

* Soft-threshold plots (PDF)
* Dendrogram plots (PNG)
* Module-trait relationship plots (PDF)
* Correlation and p-value tables (CSV)
* Hub genes (TXT)

The output files are organized into separate directories for each soft-threshold power under the results directory.

## Acknowledgements <a name = "acknowledgements"></a>
This pipeline utilizes the WGCNA package, originally developed by Peter Langfelder and Steve Horvath. To cite the original creators of WGCNA, please refer to the following publication:

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9(1), 559.
