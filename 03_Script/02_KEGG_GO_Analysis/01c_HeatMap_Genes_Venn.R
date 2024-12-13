## @knitr Multiple_Enrichment

# **Objective:** Perform multiple types of enrichment analysis and generate a heatmap
# to visualize genes that are upregulated at 4H, 24H, and those shared across both time points.

###########################################
#######   Part 1: Extracting the genes ####
###########################################

# Read gene lists from external sources for 4H and 24H intersections
Results_4H_Intersection = read.xlsx("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/IPH_6501/05_Output/01_GlobalHeterogeneity/DEG/Paired_DEG/Ven_Diagram_Gene_Lists/4H/gene_lists_4H_results.xlsx", colNames = FALSE)
Results_4H_Intersection = Results_4H_Intersection$X1

Results_24H_Intersection = read.xlsx("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/IPH_6501/05_Output/01_GlobalHeterogeneity/DEG/Paired_DEG/Ven_Diagram_Gene_Lists/24H/gene_lists_24H_results.xlsx", colNames = FALSE)
Results_24H_Intersection = Results_24H_Intersection$X1

###########################################
######   Part 2: Drawing the heatmap  #####
###########################################

# Load the Seurat Object
PBMC = readRDS(file.path(PATH_EXPERIMENT_REFERENCE_Data, SAMPLE_ANK))

# Visualize UMAP
DimPlot(PBMC)

# Subset PBMC to exclude specific conditions
PBMC = subset(PBMC, idents = c("CIML_T5_24H", "CIML_T5_4H"), invert = TRUE)
PBMC$final_annotation3 = droplevels(PBMC$final_annotation3)

# Normalize and scale the data
PBMC_Norm = NormalizeData(PBMC)
All_genes = rownames(PBMC_Norm)
PBMC = ScaleData(PBMC_Norm, features = All_genes, do.scale = DATA_SCALE, do.center = DATA_CENTER)

# Compare gene sets:
# Genes unique to 4H
List_Genes1 = setdiff(Results_4H_Intersection, Results_24H_Intersection)

# Genes shared between 4H and 24H
List_Genes2 = intersect(Results_4H_Intersection, Results_24H_Intersection)

# Genes unique to 24H
List_Genes3 = setdiff(Results_24H_Intersection, Results_4H_Intersection)

# Combine the gene lists
List_Markers = list.append(List_Genes1, List_Genes2, List_Genes3)
List_of_Lists = list(List_Genes1, List_Genes2, List_Genes3)
List_total = intersect(List_Markers, rownames(PBMC))

# Calculate average expression for each cluster
PBMC$Condition_t = droplevels(PBMC$Condition_t)
PBMC = SetIdent(PBMC, value = "Condition_t")

cluster.averages = AverageExpression(PBMC, group.by = "Condition_t", features = List_total, slot = "data")
mat = cluster.averages$RNA

# Scale the expression matrix
mat.scaled = t(apply(mat, MARGIN = 1, FUN = function(X) (X - mean(X)) / sd(X)))

# Define color function for heatmap
col_fun = circlize::colorRamp2(c(-2, -1, 0, 1, 2), rev(brewer.pal(n = 5, name = "RdBu")))

# Create heatmaps for each gene group
SIZE_ROW = 6
p1 = Heatmap(mat.scaled[intersect(List_of_Lists[[1]], rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = col_fun, border = TRUE, 
             row_names_gp = gpar(fontsize = SIZE_ROW),
             heatmap_legend_param = list(title = "Z score", direction = "horizontal", legend_width = unit(5, "cm")))
p2 = Heatmap(mat.scaled[intersect(List_of_Lists[[2]], rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = col_fun, border = TRUE, 
             show_heatmap_legend = FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))
p3 = Heatmap(mat.scaled[intersect(List_of_Lists[[3]], rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = col_fun, border = TRUE, 
             show_heatmap_legend = FALSE, row_names_gp = gpar(fontsize = SIZE_ROW))

# Combine heatmaps into a single visualization
ht_list = p1 %v% p2 %v% p3
p14 = draw(ht_list, heatmap_legend_side = "top")

# Save heatmap as PNG and PDF
png(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "HeatmapVennDG.png"), width = 10, height = 25, units = "cm", res = 600)
print(p14)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot, "HeatmapVennDG.pdf"), width = 10 / 2.54, height = 25 / 2.54)
print(p14)
dev.off()

# Additional plots for metadata visualization
VlnPlot(PBMC, features = c("GZMB"))
DimPlot(PBMC, group.by = "nkg2c")

# Heatmaps for individual gene groups (examples)
p1 = Heatmap(mat.scaled[intersect(List_Genes1, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")))
p2 = Heatmap(mat.scaled[intersect(List_Genes2, rownames(cluster.averages$RNA)),], cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)
# Additional gene lists (List_Genes3, List_Genes4, etc.) can be added similarly.

# Calculate and visualize average expression for specified gene lists
for (i in 1:length(List_of_Lists)) {
  cluster.averages = AverageExpression(PBMC, group.by = "Condition_t", features = List_of_Lists[[i]], slot = "data")
  print(pheatmap(cluster.averages$RNA, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, fontsize_col = 10, fontsize_row = 15))
}



# Additional plots for metadata visualization

# **Visualizing expression of specific genes and metadata:**
# - Generate violin and dimension plots for specific genes and metadata.

# Plot the expression of "GZMB" using a violin plot
VlnPlot(PBMC, features = c("GZMB"))

# Create a UMAP visualization grouped by the "nkg2c" metadata
DimPlot(PBMC_for_Metadata, group.by = "nkg2c")

# **Generating heatmaps for gene expression:**
# - Visualize scaled expression levels of different gene lists using heatmaps.

# Heatmap for genes unique to 4H (List_Genes1)
p1 = Heatmap(mat.scaled[intersect(List_Genes1, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")))

# Heatmap for genes shared between 4H and 24H (List_Genes2)
p2 = Heatmap(mat.scaled[intersect(List_Genes2, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)

# Additional heatmaps for other gene lists (e.g., List_Genes3, List_Genes4, etc.)
p3 = Heatmap(mat.scaled[intersect(List_Genes3, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)
p4 = Heatmap(mat.scaled[intersect(List_Genes4, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)
p5 = Heatmap(mat.scaled[intersect(List_Genes5, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)
p6 = Heatmap(mat.scaled[intersect(List_Genes6, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)
p7 = Heatmap(mat.scaled[intersect(List_Genes7, rownames(cluster.averages$RNA)),], 
             cluster_columns = FALSE, col = rev(brewer.pal(n = 4, name = "RdYlBu")), border = TRUE)

# **Calculate average expression for specific gene sets grouped by condition:**
# - Compute and visualize cluster averages for different conditions and gene lists.

cluster.averages <- AverageExpression(PBMC, group.by = "Condition_t", features = List_Demaria, slot = "data")

# Iterate through a list of gene sets and generate heatmaps
for (i in 1:6) {
  cluster.averages <- AverageExpression(PBMC, group.by = "Condition_t", features = List_of_Lists[[i]], slot = "data")
  print(pheatmap(cluster.averages$RNA, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, 
                 fontsize_col = 10, fontsize_row = 15))
}

# Heatmaps for individual gene sets
cluster.averages <- AverageExpression(PBMC, group.by = "Condition_t", features = List_Genes1, slot = "data")
p1 = pheatmap(cluster.averages$RNA, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, 
              fontsize_col = 10, fontsize_row = 15)

cluster.averages <- AverageExpression(PBMC, group.by = "Condition_t", features = List_Genes2, slot = "data")
p2 = pheatmap(cluster.averages$RNA, scale = "row", cluster_rows = TRUE, cluster_cols = FALSE, 
              fontsize_col = 10, fontsize_row = 15)

# Adjust clustering and gap rows for improved visualization
pheatmap(cluster.averages$RNA, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE, 
         fontsize_col = 15, fontsize_row = 10, cutree_rows = 6, 
         gaps_row = c(10, 16, 25, 34, 39))

# Plot specific features in violin plots
VlnPlot(Merged_Seurat_Rescaled, feature = "KLRC1")

# Display the first few rows of the RNA matrix from cluster averages
head(cluster.averages[["RNA"]][, 1:5])

# **Convert average expression into a Seurat object for visualization:**

# Clean up and reformat the metadata identifiers
orig.levels = levels(PBMC$Condition_t)
Idents(PBMC) <- gsub(pattern = " ", replacement = "_", x = Idents(PBMC))
orig.levels <- gsub(pattern = " ", replacement = "_", x = orig.levels)
levels(PBMC) <- orig.levels

# Convert average expression into a Seurat object
cluster.averages <- AverageExpression(PBMC, return.seurat = TRUE)

# Visualize average expression using a heatmap
DoHeatmap(cluster.averages, features = rownames(cluster.averages), size = 3, 
          draw.lines = FALSE, label = TRUE, group.by = "Condition_t")

# Directly visualize scaled data for selected features in the original object
DoHeatmap(PBMC, features = List_total, group.by = "Condition_t", slot = "scale.data")
