# Exploratory WGCNA Pipeline

This repository contains an R-based exploratory pipeline for Weighted Gene Co-expression Network Analysis (WGCNA). The pipeline takes gene expression and behavioral data as input and generates co-expression networks, soft-threshold plots, module-trait relationship plots, and correlation/p-value tables. The pipeline is modular, allowing users to modify or extend the steps as needed.

## Prerequisites

Before running the pipeline, please make sure you have the following R packages installed:

* `tidyverse`
* `dplyr`
* `WGCNA`
* `data.table`

## Pipeline Overview
The main pipeline script is Exploratory_WGCNA_Pipeline.R. This script sets up the necessary parameters and data directories and sources the following three modules:

1. `01_read_organize_data.R` - Reads and organizes gene expression and behavioral data.
2. `02_investigate_sft.R` - Creates soft-threshold plots to determine appropriate powers to explore.
3. `03_create_networks_export_results.R` - Performs WGCNA and exports the results.

## Usage
1. Configure the parameters at the beginning of `exploratory_wgcna_pipeline.R`. Set the study identifier, data directories, and variables for filtering out low data.
2. Indicate the group selection based on the sample ID and create a vector of soft-threshold powers to investigate.
3. Run `exploratory_wgcna_pipeline.R`.

### Module 1: Read in data and organize data for input to WGCNA
The first module reads and processes the gene expression and behavioral data. It takes care of transposing the data and filtering out low-expression genes. The module ensures that the sample IDs for both datasets are identical and prints a message if they are consistent or a warning if they do not match.

### Module 2: Create soft-threshold plots to determine powers to explore
The second module creates soft-threshold plots to determine the powers to explore. It performs network topology analysis and plots the scale-free topology fit index and mean connectivity as functions of the soft-thresholding power. The generated plots are saved as PDF files in the sft_results/plots directory, and the scale-free topology data is saved as an RData file in the sft_results/working_rdata directory.

### Module 3: Perform WGCNA and export results
The third module iterates through the specified soft-threshold powers, creating networks, calculating module eigengenes and correlations, and exporting the results. The module generates dendrogram plots, module-trait relationship plots, correlation and p-value tables, and hub genes for each power. The results are saved in separate directories for each power under the results directory.

## Output

The pipeline generates the following output files:

* Soft-threshold plots (PDF)
* Dendrogram plots (PNG)
* Module-trait relationship plots (PDF)
* Correlation and p-value tables (CSV)
* Hub genes (TXT)

The output files are organized into separate directories for each soft-threshold power under the results directory.

## License

This project is licensed under the MIT License. See the LICENSE file for details.
