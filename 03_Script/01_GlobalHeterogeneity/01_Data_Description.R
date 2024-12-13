## @knitr First_Look_All_Data

# Define a helper function to calculate the median and ± standard deviation (for diagnostic purposes)
data_summary <- function(x) {
  m <- median(x)    # Median
  ymin <- m - sd(x) # Lower bound
  ymax <- m + sd(x) # Upper bound
  return(c(y = m, ymin = ymin, ymax = ymax))
}

# Add a section header for the analysis
cat(" \n \n")
cat("# IPH 6501 analysis {.tabset .tabset-fade} \n\n")
cat(" \n \n")

cat(" \n \n")
cat("## Global look at the data")
cat(" \n \n")

# Load the PBMC dataset provided by collaborators
PBMC = readRDS((paste0(PATH_EXPERIMENT_REFERENCE_Data, "/", SAMPLE_ANK)))

# Perform dimensionality reduction and visualize clusters
p0 = DimPlot(PBMC, label = TRUE) # UMAP of all clusters
cat(" \n \n")
print(p0)
cat(" \n \n")

# Plot the data grouped by resolution of clustering
p1 = DimPlot(PBMC, group.by = "RNA_snn_res.0.5", label = TRUE) 
cat(" \n \n")
print(p1)
cat(" \n \n")

# Perform marker discovery for all clusters to identify the top genes per cluster
All_Markers = FindAllMarkers(
  PBMC, 
  only.pos = FINDMARKERS_ONLYPOS, 
  method = FINDMARKERS_METHOD, 
  min.pct = FINDMARKERS_MINPCT, 
  logfc.threshold = FINDMARKERS_LOGFC_THR, 
  verbose = TRUE
)

# Extract and order the top genes by log2 fold change for each cluster
All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

# Create a dot plot to visualize the top markers for each cluster
p1 = DotPlot(PBMC, features = unique(top10$gene), cols = "RdBu") + 
  NoLegend() + 
  theme(axis.text.x = element_text(angle = 90)) & 
  ggtitle(paste0("TOP ", FINDMARKERS_SHOWTOP, " Genes Signature"))
cat(" \n \n")
p1
cat(" \n \n")

# Scoring clusters with NK1, NK2, and NK3 signatures
# Load predefined markers for NK cell subsets
Markers_Seurat = readRDS(paste0(PATH_EXPERIMENT_REFERENCE, "/external_signatures/All_Markers_NK123.rds"))

# Filter markers based on significance, log fold change, and remove irrelevant genes
Markers_Seurat %>%
  filter(avg_log2FC > 0) %>%
  filter(p_val_adj < 5e-2) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>%
  group_by(cluster) %>%
  arrange(-avg_log2FC) %>%
  top_n(n = NUMBER_TOP_SCORING, wt = avg_log2FC) -> top_All

# Generate a list of top genes per cluster
list_Top_Genes = list()
for (i in levels(top_All$cluster)) {
  top_All %>%
    filter(cluster == i) -> top_clust
  list_Top_Genes = append(list_Top_Genes, list(top_clust$gene))
}
names(list_Top_Genes) = levels(top_All$cluster)

# Convert the list of top genes to a format usable for scoring
List_To_Use = lapply(list_Top_Genes, function(x) as.data.frame(x))
MONITORED_Markers = List_To_Use

# Add module scores for each NK subset signature
for (i in names(MONITORED_Markers)) {
  PBMC = AddModuleScore(PBMC, features = as.list(MONITORED_Markers[[i]]), pool = NULL, name = i, seed = 19)
}

# Visualize scoring with violin plots for NK1, NK2, and NK3 signatures
p1 = VlnPlot(PBMC, features = "NK11", pt.size = 0, sort = "descending") & 
  stat_summary(fun.data = data_summary, color = "black") & ggtitle("NK1")
p2 = VlnPlot(PBMC, features = "NK21", pt.size = 0, sort = "descending") & 
  stat_summary(fun.data = data_summary, color = "black") & ggtitle("NK2")
p3 = VlnPlot(PBMC, features = "NK31", pt.size = 0, sort = "descending") & 
  stat_summary(fun.data = data_summary, color = "black") & ggtitle("NK3")

cat(" \n \n")
print(p1)
cat(" \n \n")
cat(" \n \n")
print(p2)
cat(" \n \n")
cat(" \n \n")
print(p3)
cat(" \n \n")

# Generate feature plots to spatially visualize NK scores across the UMAP
p1 = FeaturePlot(PBMC, features = "NK11") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & 
  ggtitle("NK1")
p2 = FeaturePlot(PBMC, features = "NK21") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & 
  ggtitle("NK2")
p3 = FeaturePlot(PBMC, features = "NK31") & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & 
  ggtitle("NK3")

cat(" \n \n")
print(p1)
cat(" \n \n")
cat(" \n \n")
print(p2)
cat(" \n \n")
cat(" \n \n")
print(p3)
cat(" \n \n")

