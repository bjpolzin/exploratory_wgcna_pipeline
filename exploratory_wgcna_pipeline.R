### Exploratory WGCNA Pipeline ###

# Load packages
library(tidyverse)
library(WGCNA)
# Install the cleanRNAseq package from GitHub if necessary
# install.packages("remotes")
# remotes::install_github("bjpolzin/cleanRNAseq")

# Establish study identifier
study_id <- ""  # Identifier for the study, tagged onto some output results

# Set data directories
expr_dir <-  ""  # Expression data directory
# (Columns = sample IDs, rows = gene IDs)
behav_dir <- ""   # Behavior data directory
# (Columns = behaviors, rows = sample IDs)

# Set variables to filter out low data
# (refer to ?cleanRNAseq::filter_low_genes() for details)
min_expr <- 10
percent_cutoff <- 90
metric <- "units"

# Indicate group selection based on sample ID
qual <- ""  # Qualifier: string containing sample qualifier found in the
# sample ID names (e.g. tissue type, brain region, etc.)

# Create vector of soft-threshold powers to investigate
sft_powers <-  # Integer vector of soft-threshold powers to investigate
# (6-12 is a good range, but computationally expensive;
# review scale independence graph and change if necessary)

# MODULE 1: Read in data and organize data for input to WGCNA
source("./pipeline_modules/01_read_organize_data.R") 

# MODULE 2: Create soft-threshold plots to determine powers to explore
source("./pipeline_modules/02_investigate_sft.R")

# MODULE 3: Perform WGCNA and export results
source("./pipeline_modules/03_create_networks_export_results.R")
