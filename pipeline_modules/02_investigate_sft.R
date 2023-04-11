# Allow WGCNA threads (might not work on Apple Silicon Macs)
allowWGCNAThreads()

# Set options for strings
options(stringsAsFactors = FALSE)

# Define soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Perform network topology analysis
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Set the plotting window size
sizeGrWindow(9, 5)

# Set plotting parameters
par(mfrow = c(1, 2))
cex1 <- 0.9

# Save the plots to a PDF file
results_dir <- "./results/"
sft_results_dir <- "./results/sft_results/"
sft_plots_dir <- "./results/sft_results/plots/"

dir.create(results_dir)
dir.create(sft_results_dir)
dir.create(sft_plots_dir)

pdf(file = paste0(sft_plots_dir, study_sig, "_soft_thresholding_power_plots.pdf"))

# Plot scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,
     cex = cex1,
     col = "red")
abline(h = 0.80, col = "red")  # R^2 cut-off of h (should be above 0.80)

# Plot mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5],
     labels = powers,
     cex = cex1,
     col = "red")

# Close PDF file
dev.off()

# Save scale-free topology data
working_rdata_dir <- "./results/sft_results/working_rdata/"
dir.create(working_rdata_dir)
saveRDS(sft, paste0(working_rdata_dir, study_sig, "_sft_data.RData"))

# Clean up garbage
rm(sft_plots_dir, sft_results_dir)