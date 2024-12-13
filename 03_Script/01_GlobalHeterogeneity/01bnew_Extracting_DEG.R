## @knitr Extracting_DEG

# Extracting Differentially Expressed Genes (DEG) and preparing data for functional analysis.

# Section to set the purpose and context of the analysis
# The goal is to retrieve functional annotation datasets and format DEG lists
# and background lists for Over-Representation Analysis (ORA).

cat(" \n \n")
cat("## Look at Genes Lists data {.tabset .tabset-fade} \n\n")
cat(" \n \n")

# Convert the `final_annotation3` column to a factor for categorical data handling
PBMC$final_annotation3 = as.factor(PBMC$final_annotation3)

# Save all identified markers globally for later comparison or analysis
saveRDS(All_Markers, paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "/All_Markers_Allpops.rds"))

cat(" \n \n")
cat("### Genes Lists data T5 vs US in each subgroup")
cat(" \n \n")

# Define populations and timepoints for pairwise comparisons
list_pop = list("NK1", "NK2", "NK3")
list_timepoint = list("4H", "24H")

# Loop through populations and timepoints to perform DEG analysis
for (pop in list_pop) {
  for (timepoint in list_timepoint) {
    # Define groups for comparison
    group1 = paste0(pop, "_T5_", timepoint)
    group2 = paste0(pop, "_US_", timepoint)
    title = paste0(group1, "vs", group2)
    
    # Perform pairwise DEG analysis using Seurat's FindMarkers function
    comparison = FindMarkers(PBMC, ident.1 = group1, ident.2 = group2, 
                             min.pct = MIN_PCT, logfc.threshold = LOGFC_THR, verbose = FALSE)
    
    # Add gene names and cluster direction (UP or DOWN)
    comparison$gene = rownames(comparison)
    comparison <- comparison %>%
      relocate(gene, .before = everything()) %>%
      mutate(cluster = ifelse(avg_log2FC > 0, "UP", "DOWN"))
    
    # Save the DEG results
    saveRDS(comparison, paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "/Paired_DEG/", title, ".rds"))
    
    # Filter significant DEGs for visualization
    DEG_sample_filtered = comparison %>%
      filter(p_val_adj < THR_adj_pval) %>%
      arrange(desc(avg_log2FC))
    
    # Print an interactive table of the filtered DEGs
    cat(" \n \n")
    cat("#####", "Markers Table", title, "\n")
    cat(" \n \n")
    
    Nb_markers = length(unique(PBMC@active.ident))
    mypalette = hue_pal()(Nb_markers)  # Generate color palette
    print(htmltools::tagList(
      DT::datatable(DEG_sample_filtered, rownames = FALSE, extensions = 'Buttons',
                    options = list(dom = 'Blfrtip', buttons = c('excel', "csv"), fixedHeader = TRUE)) %>% 
        DT::formatStyle('cluster', backgroundColor = DT::styleEqual(c("DOWN", "UP"), PALETTE_UP_DOWN))
    ))
  }
}

# Perform additional DEG analysis for aggregated NK populations
# for specific timepoints (4H and 24H).

cat(" \n \n")
cat("### Genes Lists data all NK US_4H vs T5_4H ")
cat(" \n \n")

# Subset data for 4H and set up group comparison
data.query = subset(PBMC, subset = Timepoint == "4H")
data.query = SetIdent(data.query, value = "Condition")
group1 = "T5"
group2 = "US"
title = paste0("All_NK_", group1, "_vs_", group2, "_at_4H")

# Perform DEG analysis and annotate results
comparison = FindMarkers(data.query, ident.1 = group1, ident.2 = group2, 
                         min.pct = MIN_PCT, logfc.threshold = LOGFC_THR, verbose = FALSE)
comparison$gene = rownames(comparison)
comparison <- comparison %>%
  relocate(gene, .before = everything()) %>%
  mutate(cluster = ifelse(avg_log2FC > 0, "UP", "DOWN"))

# Save results
saveRDS(comparison, file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, paste0("Paired_DEG_", title, ".rds")))

# Visualize the filtered DEGs interactively
DEG_sample_filtered = comparison %>%
  filter(p_val_adj < THR_adj_pval) %>%
  arrange(desc(avg_log2FC))

cat(" \n \n")
cat("#####", "Markers Table", title, "\n")
cat(" \n \n")

Nb_markers = length(unique(data.query@active.ident))
mypalette = hue_pal()(Nb_markers)
print(htmltools::tagList(
  DT::datatable(DEG_sample_filtered, rownames = FALSE, extensions = 'Buttons',
                options = list(dom = 'Blfrtip', buttons = c('excel', "csv"), fixedHeader = TRUE)) %>% 
    DT::formatStyle('cluster', backgroundColor = DT::styleEqual(c("DOWN", "UP"), PALETTE_UP_DOWN))
))

cat(" \n \n")
cat("### Genes Lists data all NK US_24H vs T5_24H ")
cat(" \n \n")

# Repeat the process for 24H timepoint
data.query = subset(PBMC, subset = Timepoint == "24H")
data.query = SetIdent(data.query, value = "Condition")
group1 = "T5"
group2 = "US"
title = paste0("All_NK_", group1, "_vs_", group2, "_at_24H")

comparison = FindMarkers(data.query, ident.1 = group1, ident.2 = group2, 
                         min.pct = MIN_PCT, logfc.threshold = LOGFC_THR, verbose = FALSE)
comparison$gene = rownames(comparison)
comparison <- comparison %>%
  relocate(gene, .before = everything()) %>%
  mutate(cluster = ifelse(avg_log2FC > 0, "UP", "DOWN"))

# Save results
saveRDS(comparison, file.path(PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, paste0("Paired_DEG_", title, ".rds")))

# Visualize the filtered DEGs interactively
DEG_sample_filtered = comparison %>%
  filter(p_val_adj < THR_adj_pval) %>%
  arrange(desc(avg_log2FC))

cat(" \n \n")
cat("#####", "Markers Table", title, "\n")
cat(" \n \n")

Nb_markers = length(unique(data.query@active.ident))
mypalette = hue_pal()(Nb_markers)
print(htmltools::tagList(
  DT::datatable(DEG_sample_filtered, rownames = FALSE, extensions = 'Buttons',
                options = list(dom = 'Blfrtip', buttons = c('excel', "csv"), fixedHeader = TRUE)) %>% 
    DT::formatStyle('cluster', backgroundColor = DT::styleEqual(c("DOWN", "UP"), PALETTE_UP_DOWN))
))

