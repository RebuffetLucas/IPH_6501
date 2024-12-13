# A few additional figures

# **UMAP Visualization of Treatment Conditions**
# This section visualizes the UMAP of PBMCs colored by treatment conditions (T5 vs US).

# Load the Seurat object for PBMC data
PBMC = readRDS(file.path(PATH_EXPERIMENT_REFERENCE_Data, SAMPLE_ANK))

# Generate a UMAP plot colored by "Condition" metadata and remove axes for a cleaner plot
p0 = DimPlot(PBMC, label = TRUE, group.by = "Condition") & NoAxes()

# Save the UMAP plot as a PNG file
png(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "UMAP_Treatment.png"), 
    width = 12, height = 10, units = "cm", res = 600)
print(p0)
dev.off()

# Save the UMAP plot as a PDF file
pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "UMAP_Treatment.pdf"), 
    width = 12 / 2.54, height = 10 / 2.54)
print(p0)
dev.off()

# **CIML Score Visualization**
# This section visualizes the CIML score as a feature plot to show its distribution across cells.

# Generate a feature plot for the "CIML" score with a reversed "RdBu" color gradient and no axes
p1 = FeaturePlot(PBMC, label = FALSE, features = "CIML") & 
  NoAxes() & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Display the feature plot
p1

# Save the CIML score plot as a PNG file
png(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "CIML_score.png"), 
    width = 12, height = 10, units = "cm", res = 600)
print(p1)
dev.off()

# Save the CIML score plot as a PDF file
pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "CIML_score.pdf"), 
    width = 12 / 2.54, height = 10 / 2.54)
print(p1)
dev.off()

# Print the UMAP treatment plot again (if needed for immediate display in the script)
p0
