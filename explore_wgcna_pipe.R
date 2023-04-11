### Exploratory WGCNA Pipeline ###

## Load packages
library(tidyverse)
library(WGCNA)
# install.packages("remotes")
# remotes::install_github("bjpolzin/cleanRNAseq")

## Study signifier
study_sig <- "mpoa_play"

## Set your data directories
expr_dir <-  "./data/rsem_play_data.txt" # Directory where expression data is located
                                         # Columns = sample IDs, rows = gene IDs
behav_dir <- "./data/play_behavior.txt"  # Directory where behavior data is located
                                         # Columns = behaviors, rows = sample IDs

## Set variables to filter out low data
min_expr <- 10
percent_cutoff <- 90
metric <- "FPKM"

## Indicate group selection based on sample ID
qual <- "mpoa" # Qualifier: string containing sample qualifier found in the 
               # sample ID names (e.g. tissue type, brain region, etc.)

## Create vector of all soft-threshold powers you want to investigate ##
sft_power <- 6:12 # Integer vector of all soft-threshold powers you want to 
                  # investigate (6-12 good range, but computationally expensive, 
                  # review scale independence graph and change if necessary)

#### MODULE 1: Read in data and organize data for input to WGCNA 
source("./pipeline_modules/01_read_organize_data.R") 
# Review message for how many genes were filtered out, in addition to ensuring
# that sample IDs are consistent between both datasets

#### MODULE 2: Create soft-threshold plots to determine powers you want to explore 
source("./pipeline_modules/02_investigate_sft.R")

#### MODULE 3: Perform WGCNA for selected powers
source("")














### DIRECTORIES WILL NEED TO BE CREATED ##
orig_wd <- getwd()
### Create networks and correlate behaviors for all powers desired ##
for (beta_threshold in threshold_powers) {
  # Create appropriate directories
  power_wd <- file.path(".", paste("power_", beta_threshold, sep = ""))
  dir.create(power_wd)
  setwd(power_wd)
  # Construct the gene network and identifying modules
  net = blockwiseModules(datExpr, power = beta_threshold, TOMType = "signed", 
                         minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.3, 
                         numericLabels = TRUE, pamRespectsDendro = FALSE, 
                         saveTOMs = FALSE, saveTOMFileBase = "mydataTOM1", 
                         verbose = 3, deepSplit=2)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(net$colors)
  # Write out dendro plot to appropriate folder
  png(file = paste("power_", beta_threshold, "_dendrogram.png"), 
      width = 900, height = 500)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], 
                      "Module colors", dendroLabels = FALSE, hang = 0.03, 
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs
  geneTree = net$dendrograms[[1]]
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)
  
  MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3))
  # Write module trait plot
  pdf(file = paste("power_", beta_threshold, "_module_trait_plot.pdf", sep = ""),
      height = 10, width = 13)
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), 
                 ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), 
                 textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.6, zlim = c(-1,1), 
                 main = paste("Module-trait relationships"), cex.lab.y = 0.5, cex.lab.x = 0.6)
  dev.off()
  # Write out csv of cor values
  cor_df <- as.data.frame(moduleTraitCor)
  write.csv(cor_df, file = paste("power_", beta_threshold, "_cor_value_table.csv",
                                 sep = ""))
  # Write out csv of p values
  p_vals_df <- as.data.frame(moduleTraitPvalue)
  write.csv(p_vals_df, file = paste("power_", beta_threshold, "_p_vals_table.csv",
                                    sep = ""))
  # Grab hub genes
  hubs <- chooseTopHubInEachModule(datExpr, omitColors = "grey", colorh = moduleColors,
                                   power = beta_threshold, type = "signed") %>%
    as.data.frame()
  write.csv(hubs, file = paste("power_", beta_threshold, "_hub_genes.csv",
                               sep = ""))
  
  # Grab hub genes
  setwd(orig_wd)
}


# Could be worth it to grab number of significant modules for all of the given behaviors
# Could be worth it to investigate number of modules that are broadly created in each
# Could be worth it to export the hub genes here too for the summary of results





