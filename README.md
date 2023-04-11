# Table of Contents
1. [Exploratory WGCNA Pipeline](#exploratory-wgcna-pipeline)
2. [Prerequisites](#prerequisites)
3. [Pipeline Overview](#pipeline-overview)
4. [Usage](#usage)
5. [Output](#output)
6. [License](#license)
7. [Acknowledgements](#acknowledgements)

# Exploratory WGCNA Pipeline <a name="exploratory-wgcna-pipeline"></a>
This repository contains an R-based exploratory pipeline for Weighted Gene Co-expression Network Analysis (WGCNA), designed for initial exploratory analyses to compare the results of different soft threshold powers. The pipeline takes gene expression and behavioral data as input and generates co-expression networks, soft-threshold plots, module-trait relationship plots, and correlation/p-value tables. The pipeline is modular, allowing users to modify or extend the steps as needed.

This is all based on the WGCNA codebase originally created by Peter Langfelder and Steve Horvath. Please see acknowledgements at the bottom of this page for citation and refer to [this page](https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/) if you are interested in learning more about the WGCNA method.

## Prerequisites <a name="prerequisites"></a>
Before running the pipeline, please make sure you have the following R packages installed:

* `tidyverse`
* `dplyr`
* `WGCNA`
* `data.table`

## Pipeline Overview <a name="pipeline-overview"></a>
The main pipeline script is `exploratory_wgcna_pipeline.`R. This script sets up the necessary parameters and data directories and sources the following three modules:

1. `01_read_organize_data.R` - Reads and organizes gene expression and behavioral data.
2. `02_investigate_sft.R` - Creates soft-threshold plots to determine appropriate powers to explore.
3. `03_create_networks_export_results.R` - Performs WGCNA and exports the results.

## Usage <a name="usage"></a>
1. Configure the parameters at the beginning of `exploratory_wgcna_pipeline.R`. Set the study identifier, data directories, and variables for filtering out low data.
2. Indicate the group selection based on the sample ID and create a vector of soft-threshold powers to investigate.
3. Run `exploratory_wgcna_pipeline.R`.

### Module 1: Read in data and organize data for input to WGCNA
The first module reads and processes the gene expression and behavioral data. It takes care of transposing the data and filtering out low-expression genes. The module ensures that the sample IDs for both datasets are identical and prints a message if they are consistent or a warning if they do not match.

### Module 2: Create soft-threshold plots to determine powers to explore
The second module creates soft-threshold plots to determine the powers to explore. It performs network topology analysis and plots the scale-free topology fit index and mean connectivity as functions of the soft-thresholding power. The generated plots are saved as PDF files in the `sft_results/plots directory`, and the scale-free topology data is saved as an RData file in the `sft_results/working_rdata` directory.

### Module 3: Perform WGCNA and export results
The third module iterates through the specified soft-threshold powers, creating networks, calculating module eigengenes and correlations, and exporting the results. The module generates dendrogram plots, module-trait relationship plots, correlation and p-value tables, and hub genes for each power. The results are saved in separate directories for each power under the `results` directory.

## Output <a name="output"></a>
The pipeline generates the following output files in automatically generated folders:

* Soft-threshold plots (PDF)
* Dendrogram plots (PNG)
* Module-trait relationship plots (PDF)
* Correlation and p-value tables (CSV)
* Hub genes (TXT)

The output files are organized into separate directories for each soft-threshold power under the results directory.

## License <a name="license"></a>
This project is licensed under the MIT License. See the LICENSE file for details.

## Acknowledgements <a name = "acknowledgements"></a
This pipeline utilizes the WGCNA package, originally developed by Peter Langfelder and Steve Horvath. To cite the original creators of WGCNA, please refer to the following publication:

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics, 9(1), 559.
