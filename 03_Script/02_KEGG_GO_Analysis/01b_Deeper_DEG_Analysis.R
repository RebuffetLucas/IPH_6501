## @knitr Multiple_Enrichment
#Perform Multiple kind of enrichment analysis and store them

#Part 1:
#  The aim of this script is to retrieve Functionnal Annotation Datasets
#  & Format both these datasets, DEG lists and Background Lists for ORA

#install.packages("msigdbr")
#BiocManager::install("ReactomePA")
#install.packages("rlist")


#Build the list of comparison
list_pop= list( "NK1", "NK2", "NK3")
list_timepoint = list("4H", "24H")

list_paired_rds = list()
list_titles= list()
for (pop in list_pop){
  for (timepoint in list_timepoint){
    group1= paste0(pop , "_T5_", timepoint)
    group2= paste0(pop , "_US_", timepoint)
    title = paste0(group1, "vs" , group2)
    
    comparison =  readRDS(file.path( PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "Paired_DEG" , paste0(title, ".rds" )))
    list_paired_rds = list.append(list_paired_rds, comparison)
    list_titles = list.append(list_titles, title)
  }
}

names(list_paired_rds) = list_titles


#Sort the genes and store them in a table
my_list = list_paired_rds

# Function to extract genes based on the conditions
extract_genes <- function(df) {
  df %>%
    filter(cluster == "UP" & avg_log2FC > TRESHOLD_LOGFC_SELECTION) %>%
    select(gene)
}


#Extract genes
gene_lists <- lapply(my_list, extract_genes)


# Convert each data frame in the result list to a character vector of gene names
gene_lists <- lapply(gene_lists, function(df) df$gene)



############################
###### Ven Diagram #########
############################

# Separate the gene_lists into 4H and 24H conditions
gene_lists_4H <- gene_lists[c("NK1_T5_4HvsNK1_US_4H", "NK2_T5_4HvsNK2_US_4H", "NK3_T5_4HvsNK3_US_4H")]
gene_lists_24H <- gene_lists[c("NK1_T5_24HvsNK1_US_24H", "NK2_T5_24HvsNK2_US_24H", "NK3_T5_24HvsNK3_US_24H")]

# Create labels for the Venn diagrams
labels_4H <- c("NK1", "NK2", "NK3")
labels_24H <- c("NK1", "NK2", "NK3")

# Draw the Venn diagram for 4H conditions
venn_diagram_4H <- venn.diagram(
  x = gene_lists_4H,
  category.names = labels_4H,
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = palette3,
  main = "Venn Diagram for 4H Conditions"
)


# Display the Venn diagram for 4H conditions
grid.newpage()
grid.draw(venn_diagram_4H)


#Save Figure
png(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG , "Venn_Diagram_4H_HD.png"), width = 10, height = 10,  units = "cm", res=600 )
grid.draw(venn_diagram_4H)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG, "Venn_Diagram_4H_HD.pdf"), width = 10/2.54, height = 10/2.54)
grid.draw(venn_diagram_4H)
dev.off()




# Draw the Venn diagram for 24H conditions
venn_diagram_24H <- venn.diagram(
  x = gene_lists_24H,
  category.names = labels_24H,
  filename = NULL,
  output = TRUE,
  imagetype = "png",
  height = 480,
  width = 480,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  lty = 'blank',
  fill = palette3,
  main = "Venn Diagram for 24H Conditions"
) 

# Display the Venn diagram for 24H conditions
grid.newpage()
grid.draw(venn_diagram_24H)



#Save Figure
png(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG , "Venn_Diagram_24H_HD.png"), width = 10, height = 10,  units = "cm", res=600 )
grid.draw(venn_diagram_24H)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG, "Venn_Diagram_24H_HD.pdf"), width = 10/2.54, height = 10/2.54)
grid.draw(venn_diagram_24H)
dev.off()

###############################################
####  Store the intersections of the lists ####
###############################################


# For 4H
gene_lists_4H <- gene_lists[c("NK1_T5_4HvsNK1_US_4H", "NK2_T5_4HvsNK2_US_4H", "NK3_T5_4HvsNK3_US_4H")]

# Extract the genes from each sublist
genes_NK1 <- gene_lists_4H$NK1_T5_4HvsNK1_US_4H
genes_NK2 <- gene_lists_4H$NK2_T5_4HvsNK2_US_4H
genes_NK3 <- gene_lists_4H$NK3_T5_4HvsNK3_US_4H

# Convert to character vectors
genes_NK1 <- as.character(genes_NK1)
genes_NK2 <- as.character(genes_NK2)
genes_NK3 <- as.character(genes_NK3)

# Find the intersection of the three sublists
intersect_genes <- Reduce(intersect, list(genes_NK1, genes_NK2, genes_NK3))

# Find the unique elements in each sublist
unique_NK1 <- setdiff(genes_NK1, union(genes_NK2, genes_NK3))
unique_NK2 <- setdiff(genes_NK2, union(genes_NK1, genes_NK3))
unique_NK3 <- setdiff(genes_NK3, union(genes_NK1, genes_NK2))

# Create a data frame to store the results
results <- list(
  "Intersection" = intersect_genes,
  "Unique_NK1" = unique_NK1,
  "Unique_NK2" = unique_NK2,
  "Unique_NK3" = unique_NK3
)

# Write the results to an Excel file
write.xlsx(results, file = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/IPH_6501/05_Output/01_GlobalHeterogeneity/DEG/Paired_DEG/Ven_Diagram_Gene_Lists/4H/gene_lists_4H_results.xlsx")


# For 24H
gene_lists_24H <- gene_lists[c("NK1_T5_24HvsNK1_US_24H", "NK2_T5_24HvsNK2_US_24H", "NK3_T5_24HvsNK3_US_24H")]

# Extract the genes from each sublist
genes_NK1 <- gene_lists_24H$NK1_T5_24HvsNK1_US_24H
genes_NK2 <- gene_lists_24H$NK2_T5_24HvsNK2_US_24H
genes_NK3 <- gene_lists_24H$NK3_T5_24HvsNK3_US_24H

# Convert to character vectors
genes_NK1 <- as.character(genes_NK1)
genes_NK2 <- as.character(genes_NK2)
genes_NK3 <- as.character(genes_NK3)

# Find the intersection of the three sublists
intersect_genes <- Reduce(intersect, list(genes_NK1, genes_NK2, genes_NK3))

# Find the unique elements in each sublist
unique_NK1 <- setdiff(genes_NK1, union(genes_NK2, genes_NK3))
unique_NK2 <- setdiff(genes_NK2, union(genes_NK1, genes_NK3))
unique_NK3 <- setdiff(genes_NK3, union(genes_NK1, genes_NK2))

# Create a data frame to store the results
results <- list(
  "Intersection" = intersect_genes,
  "Unique_NK1" = unique_NK1,
  "Unique_NK2" = unique_NK2,
  "Unique_NK3" = unique_NK3
)

write.xlsx(results, file = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/IPH_6501/05_Output/01_GlobalHeterogeneity/DEG/Paired_DEG/Ven_Diagram_Gene_Lists/24H/gene_lists_24H_results.xlsx")



############################
####  Heatmap Jacquard ####
############################

#Calculate Jaccard Index
#Function to calculate Jaccard Index
jaccard_index <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

#Calculate Jaccard Index for each combination of gene lists
jaccard_results <- matrix(0, nrow = length(gene_lists), ncol = length(gene_lists))
rownames(jaccard_results) <- names(gene_lists)
colnames(jaccard_results) <- names(gene_lists)

for (i in seq_along(gene_lists)) {
  for (j in seq_along(gene_lists)) {
    jaccard_results[i, j] <- jaccard_index(gene_lists[[i]], gene_lists[[j]])
  }
}

# Perform hierarchical clustering
dist_matrix <- as.dist(1 - jaccard_results)  # Convert similarity to distance
hc <- hclust(dist_matrix, method = "average")  # Hierarchical clustering

# Reorder the Jaccard index matrix according to clustering
jaccard_results <- jaccard_results[hc$order, hc$order]

# Print the Jaccard Index matrix
print(jaccard_results) #Have a look at the results

# Convert the matrix to a long format suitable for ggplot
jaccard_data_long <- as.data.frame(as.table(jaccard_results))

# Renaming the columns for clarity
colnames(jaccard_data_long) <- c("Row", "Column", "JaccardIndex")

# Create the heatmap
heat_jaccard <- ggplot(jaccard_data_long, aes(x = Row, y = Column, fill = JaccardIndex)) +
  geom_tile(color = "white") +  # Use geom_tile() for heatmap squares
  scale_fill_viridis() +
  theme_minimal() +  # Use a minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        axis.title = element_blank()) +  # Remove default axis titles
  labs(fill = "Jaccard\nIndex", title = "Jaccard Index Heatmap") +  # Set fill label and title
  coord_fixed()  # Ensure the aspect ratio is preserved

print(heat_jaccard)


#Save Figure
png(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG , "Jaccard_HD.png"), width = 15, height = 15,  units = "cm", res=600 )
print(heat_jaccard)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG, "Jaccard_HD.pdf"), width = 15/2.54, height = 15/2.54)
print(heat_jaccard)
dev.off()



############################
## Heatmap Overlap Index  ##
############################


# Function to calculate Overlap Index
overlap_index <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection_size <- length(intersect(set1, set2))
  min_size <- min(length(set1), length(set2))
  return(intersection_size / min_size)
}

# Calculate Overlap Index for each combination of gene lists
overlap_results <- matrix(0, nrow = length(gene_lists), ncol = length(gene_lists))
rownames(overlap_results) <- names(gene_lists)
colnames(overlap_results) <- names(gene_lists)

for (i in seq_along(gene_lists)) {
  for (j in seq_along(gene_lists)) {
    overlap_results[i, j] <- overlap_index(gene_lists[[i]], gene_lists[[j]])
  }
}

# Perform hierarchical clustering
dist_matrix <- as.dist(1 - overlap_results)  # Convert similarity to distance
hc <- hclust(dist_matrix, method = "average")  # Hierarchical clustering

# Reorder the Overlap index matrix according to clustering
overlap_results <- overlap_results[hc$order, hc$order]

# Print the Overlap Index matrix
print(overlap_results) # Have a look at the results

# Convert the matrix to a long format suitable for ggplot
overlap_data_long <- as.data.frame(as.table(overlap_results))

# Renaming the columns for clarity
colnames(overlap_data_long) <- c("Row", "Column", "OverlapIndex")

# Create the heatmap
heat_overlap <- ggplot(overlap_data_long, aes(x = Row, y = Column, fill = OverlapIndex)) +
  geom_tile(color = "white") +  # Use geom_tile() for heatmap squares
  scale_fill_viridis() +
  theme_minimal() +  # Use a minimal theme for cleaner appearance
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for better readability
        axis.title = element_blank()) +  # Remove default axis titles
  labs(fill = "Overlap\nIndex", title = "Overlap Index Heatmap") +  # Set fill label and title
  coord_fixed()  # Ensure the aspect ratio is preserved

print(heat_overlap)


#Save Figure
png(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG , "Overlap_HD.png"), width = 15, height = 15,  units = "cm", res=600 )
print(heat_overlap)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG, "Overlap_HD.pdf"), width = 15/2.54, height = 15/2.54)
print(heat_overlap)
dev.off()


