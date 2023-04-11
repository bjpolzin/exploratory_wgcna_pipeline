# Create networks and export results
for (beta_threshold in sft_powers) {
  
  # Create unique directory name
  # power_dir <- paste0("./results/power_", beta_threshold)
  # Create result output markers
  power_mark <- paste0(study_sig, "_power_", beta_threshold, "_")
  power_results_dir <- paste0("./results/", "power_", beta_threshold, "_results/")
  power_plots_dir <- paste0(power_results_dir, "plots/")
  
  
  
  # Create the directory if it doesn't exist
  
  if (!file.exists(power_results_dir)) {
    dir.create(power_results_dir)
  }
  
  if (!file.exists(power_plots_dir)) {
    dir.create(power_plots_dir)
  }
  
  # setwd(power_wd)
  
  # Construct the gene network and identifying modules
  net <- blockwiseModules(
    datExpr,
    power = beta_threshold,
    TOMType = "signed",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.3,
    numericLabels = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs = FALSE,
    saveTOMFileBase = "mydataTOM1",
    verbose = 3,
    deepSplit = 2,
    maxBlockSize = 20000
  )
  
  # Convert labels to colors for plotting
  mergedColors <- labels2colors(net$colors)
  
  # Save dendrogram plot
  png(file = paste0(power_plots_dir, "power_", beta_threshold, "_dendrogram.png"),
      width = 900, height = 500)
  plotDendroAndColors(
    net$dendrograms[[1]],
    mergedColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05
  )
  dev.off()
  
  # Calculate module eigengenes and correlations
  moduleLabels <- net$colors
  moduleColors <- labels2colors(net$colors)
  MEs <- net$MEs
  geneTree <- net$dendrograms[[1]]
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  # Save module trait plot
  pdf(file = paste0(power_plots_dir, "power_", beta_threshold, "_module_trait_plot.pdf"),
      height = 10, width = 13)
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = names(datTraits),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.6,
    zlim = c(-1, 1),
    main = paste("Module-trait relationships"),
    cex.lab.y = 0.5,
    cex.lab.x = 0.6
  )
  dev.off()
  
  # Save correlation and p-value tables
  cor_df <- as.data.frame(moduleTraitCor)
  write.csv(cor_df, file = paste0(power_results_dir, "power_", beta_threshold, "_cor_value_table.csv"))
  p_vals_df <- as.data.frame(moduleTraitPvalue)
  write.csv(p_vals_df, file = paste0(power_results_dir, "power_", beta_threshold, "_p_vals_table.csv"))
  
  # Save hub genes
  hubs <- chooseTopHubInEachModule(datExpr, omitColors = "grey", colorh = moduleColors,
                                   power = beta_threshold, type = "signed")
  write.table(hubs, file = paste0(power_plots_dir, "power_", beta_threshold, 
                                  "_hub_genes.txt"), quote = F, col.names = F)
}

rm(power_plots_dir, power_results_dir, power_mark)
