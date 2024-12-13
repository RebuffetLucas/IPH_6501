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




#Improved version
# Create a function to generate the plot for each data frame in the list
plot_function <- function(data, title) {
  # Count the number of genes that pass the thresholds
  n_0_5 <- nrow(data[data$avg_log2FC >= 0.5, ])
  n_0_75 <- nrow(data[data$avg_log2FC >= 0.75, ])
  
  # Determine the maximum y value for placing the labels
  max_y <- min(300 , max(-log10(data$p_val_adj), na.rm = TRUE))
  
  ggplot(data, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(size = 1) +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "red") +
    geom_vline(xintercept = 0.75, linetype = "dotted", color = "blue") +
    annotate("text", x = 0.5, y = max_y * 0.9, label = paste("n =", n_0_5), color = "red", vjust = -0.5) +
    annotate("text", x = 0.75, y = max_y * 0.8, label = paste("n =", n_0_75), color = "blue", vjust = -0.5) +
    labs(title = title, x = "avg_log2FC", y = "-log10(p_val_adj)") +
    theme_minimal()
}

# Generate a list of plots
plots <- lapply(names(list_paired_rds), function(name) {
  plot_function(list_paired_rds[[name]], name)
})

# Arrange the plots in a grid
grid.arrange(grobs = plots, ncol = 2)



#Save Figure
png(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG , "Select_LOGFC_THR.png"), width = 40, height = 40,  units = "cm", res=600 )
grid.arrange(grobs = plots, ncol = 2)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG, "Select_LOGFC_THR.pdf"), width = 40/2.54, height = 40/2.54)
grid.arrange(grobs = plots, ncol = 2)
dev.off()


