## @knitr Multiple_Enrichment
#Perform Multiple kind of enrichment analysis and store them

#Part 1:
#  The aim of this script is to retrieve Functionnal Annotation Datasets
#  & Format both these datasets, DEG lists and Background Lists for ORA

#install.packages("msigdbr")
#BiocManager::install("ReactomePA")
#install.packages("rlist")


#Build the list of comparison with the 6 subpops
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

#Add the two other comparisons (all NK at 4H and All NK at 24H)

#4H
comparison =  readRDS(file.path( PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "Paired_DEG_All_NK_T5_vs_US_at_4H.rds" ))
list_paired_rds = list.append(list_paired_rds, comparison)
names(list_paired_rds)[7] = "AllNK_T5_4HvsAllNK_US_4H"

#24H
comparison =  readRDS(file.path( PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "Paired_DEG_All_NK_T5_vs_US_at_24H.rds" ))
list_paired_rds = list.append(list_paired_rds, comparison)
names(list_paired_rds)[8] = "AllNK_T5_24HvsAllNK_US_24H"


#First, Filter all the data based on the logFC Treshold determine before
filtered_list_paired_rds <- lapply(list_paired_rds, function(df) {
  df[abs(df$avg_log2FC) > TRESHOLD_LOGFC_SELECTION, ]
})



#Define universe
PBMC = readRDS((paste0(PATH_EXPERIMENT_REFERENCE_Data,"/",SAMPLE_ANK)))
UNIVERSE_SymbolID = rownames(PBMC@assays[["RNA"]])



#UNIVERSE_SymbolID = VariableFeatures(PBMC)


#THR_adj_pval = 0.05
#THR_log2FC = 0

FX_TO_RUN = "enrichPathway" # One of "groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway" .



for (indice in 1:length(filtered_list_paired_rds)){

Markers= filtered_list_paired_rds[[indice]]
title = names(filtered_list_paired_rds)[indice]

# Search for ENTREZID equivalent ID
##### ADD ENTREZID TO ORA DF #####
# Rename Marker column to SYMBOL
ORA_DF <- Markers %>% dplyr::select(cluster,gene) %>% dplyr::rename(SYMBOL=gene)


# Conversion Gene Symbol to ENTREID for some functional enrichment analysis (such enrichKEGG)
genes_EntrezID<-AnnotationDbi::select(org.Hs.eg.db, ORA_DF$SYMBOL, 'ENTREZID','SYMBOL')
genes_EntrezID<-genes_EntrezID[which(genes_EntrezID$ENTREZID != is.na(genes_EntrezID$ENTREZID)),]
ORA_DF<-ORA_DF %>% filter(SYMBOL %in% genes_EntrezID$SYMBOL)
ORA_DF<-inner_join(ORA_DF,genes_EntrezID, by="SYMBOL") %>% dplyr::distinct()




#####################################################################################
## Prepare Background list for SYMBOL ID functionnal datasets ORA  ####
#####################################################################################
## @knitr format_background_data

# Background list conversion for ENTREID functional enrichment analysis (such enrichKEGG)
UNIVERSE_EntrezID <- AnnotationDbi::select(org.Hs.eg.db, UNIVERSE_SymbolID, 'ENTREZID','SYMBOL')
UNIVERSE_EntrezID <- UNIVERSE_EntrezID[which(UNIVERSE_EntrezID$ENTREZID != is.na(UNIVERSE_EntrezID$ENTREZID)),]
UNIVERSE_SymbolID <- UNIVERSE_EntrezID$SYMBOL
UNIVERSE_EntrezID <- UNIVERSE_EntrezID$ENTREZID

#######################
## Compare Cluster ####
#######################

#Define list up and list down:
# Create list_entrez_UP containing ENTREZID for cluster "UP"
list_entrez_UP <- ORA_DF$ENTREZID[ORA_DF$cluster == "UP"]

# Create list_entrez_DOWN containing ENTREZID for cluster "DOWN"
list_entrez_DOWN <- ORA_DF$ENTREZID[ORA_DF$cluster == "DOWN"]


# Create an empty list to store the 2 pathway analysis
CompareCLUST <- vector(mode = "list", length = 2)
names(CompareCLUST)<-c(paste0( title, "UP" ) , paste0(title, "DOWN" ))



#Run comparison:

CompareCLUST[[1]] <- enrichPathway(
  gene= list_entrez_UP,
  organism = "human",
  pvalueCutoff = ENRICHPATHWAY_P_VALUE_CUTOFF,
  pAdjustMethod = "BH",
  qvalueCutoff = ENRICHPATHWAY_Q_VALUE_CUTOFF,
  universe= UNIVERSE_EntrezID ,
  minGSSize = ENRICHPATHWAY_minGS_SIZE,
  maxGSSize = ENRICHPATHWAY_maxGS_SIZE,
  readable = TRUE
)

CompareCLUST[[2]] <- enrichPathway(
  gene= list_entrez_DOWN,
  organism = "human",
  pvalueCutoff = ENRICHPATHWAY_P_VALUE_CUTOFF,
  pAdjustMethod = "BH",
  qvalueCutoff = ENRICHPATHWAY_Q_VALUE_CUTOFF,
  universe= UNIVERSE_EntrezID ,
  minGSSize = ENRICHPATHWAY_minGS_SIZE,
  maxGSSize = ENRICHPATHWAY_maxGS_SIZE,
  readable = TRUE
)


#Save 
saveRDS(CompareCLUST,paste0( PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG, "/Paired_Enriched_Pathway/", title, ".rds" ))
gc()

print(paste0("list" , indice, " finished"))
}
