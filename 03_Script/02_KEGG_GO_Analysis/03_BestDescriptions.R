## @knitr Best_Enrichment

#############################################################
## Find descriptions of interest####
###############################################################

#Read all the enrichment results:

#Build the list of enrichment with the 6 subpops
list_pop= list( "NK1", "NK2", "NK3", "AllNK")
list_timepoint = list("4H", "24H")

list_paired_rds = list()
list_titles= list()
for (pop in list_pop){
  for (timepoint in list_timepoint){
    group1= paste0(pop , "_T5_", timepoint)
    group2= paste0(pop , "_US_", timepoint)
    title = paste0(group1, "vs" , group2)
    
    
    enrichment =  readRDS(file.path( PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "Paired_Enriched_Pathway" , paste0(title, ".rds" )))
    
    list_paired_rds = list.append(list_paired_rds, enrichment)
    list_titles = list.append(list_titles, title)
  }
}

names(list_paired_rds) = list_titles


CompareCLUST_list = list_paired_rds

for (i in 1:length(CompareCLUST_list)){
  for (j in 1:length(CompareCLUST_list[[1]])){
  
  TITLE_comparison  = names(CompareCLUST_list[[i]])[j] 
  
  print(TITLE_comparison)
  
  object_plot = CompareCLUST_list[[i]][[j]]
  
  if (dim(object_plot)[1]!=0 ){
  
  DotPlot_Kegg=  enrichplot::dotplot(object_plot, font.size = 7, label_format = 90, showCategory=30 , title= TITLE_comparison) + 
    theme(axis.text.x=element_text(size = 7), axis.text.y=element_text(size = 7), 
        legend.text =  element_text(size = 7), legend.title =  element_text(size = 7))
  
  #png(file = file.path( PATH_EXPERIMENT_OUTPUT_FIGURE_KEGG_DotPlot , paste0("DotPlot_KEGG_", TITLE_comparison,".png")), width = 30, height = 20,  units = "cm", res=600 )
  #print(DotPlot_Kegg)
  #dev.off()
  
  pdf(file = file.path( PATH_EXPERIMENT_OUTPUT_FIGURE_KEGG_DotPlot , paste0("DotPlot_KEGG_", TITLE_comparison,".pdf")), width = 30/2.54, height = 20/2.54)
  print(DotPlot_Kegg)
  dev.off()
  
  }
  }
}




GO_results = CompareCLUST[[COMP_OF_INTEREST_NAME]]

Go_results_CompareClust = GO_results@compareClusterResult

#Have a look at Shared Descriptions

extract_duplicates <- function(lst) {
  table_lst <- table(lst)
  duplicates <- names(table_lst[table_lst > 1])
  return(duplicates)
}


Shared_Descriptions = extract_duplicates(Go_results_CompareClust$Description)


#all_Descriptions[ Go_results_CompareClust$Description=="oxidative phosphorylation" ] 
#Remove unwanted terms
all_Descriptions = unique(Go_results_CompareClust$Description)

all_Descriptions

black_List2 <- data.frame(all_Descriptions) %>%
  filter( grepl(TERMS_TO_REMOVE, all_Descriptions))

black_List2= black_List2$all_Descriptions

black_List = black_List2

Go_results_CompareClust2 = filter(Go_results_CompareClust,  !(Description %in% black_List))


Go_results_CompareClust2$GeneRatio = parse_ratio(Go_results_CompareClust2$GeneRatio)

top_descriptions <- Go_results_CompareClust2 %>%
  group_by(Cluster) %>%
  slice_min(p.adjust, n = NUMBER_OF_DESCRIPTIONS_TO_KEEP)



#Extract in Top10 the most interesting 
unique(top_descriptions$Description)


