

for(line in lines){
  cat(paste0("## ", line, " {.tabset .tabset-fade}  \n \n"))
  cat("### List of genes {.tabset .tabset-fade}  \n \n")
  print(tagList(DT::datatable(as.matrix(setdiff(intersect(up_drug2[,line], up_drug1[,line]), list_pan)))))
  cat("\n \n")
  
  cat("### Pathway Analysis (KEGG) {.tabset .tabset-fade} \n \n")cat("### Pathway Analysis (KEGG) {.tabset .tabset-fade} \n \n")

############################################################### EnrichKEGG Analysis

liste = as.character(as.matrix(setdiff(intersect(up_drug2[,line], up_drug1[,line]), list_pan)))
#library(organism, character.only = TRUE)

entrez_id = as.character(na.omit(mapIds(org.Hs.eg.db, na.omit(liste), 'ENTREZID', 'SYMBOL')))
gse = enrichKEGG(entrez_id,organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE) #universe=universe_blue_entrez, 

if(length(gse)>2){
  # Plot
  cat("#### DotPlot  {.tabset .tabset-fade} \n \n")
  print(enrichplot::dotplot(gse, showCategory=10,font.size=18)) #, split=".sign" + facet_grid(.~.sign)
  cat("\n \n")
  
  cat("#### cNET Plot  {.tabset .tabset-fade} \n \n")
  edox <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
  print(cnetplot(edox, showCategory=10,font.size=18))
  cat("\n \n")
  
  cat("#### Table  {.tabset .tabset-fade} \n \n")
  gene_names=rep("NA", nrow(summary(gse)))
  for(i in 1:nrow(summary(gse))){gene_names[i] = paste(as.character(mapIds(org.Hs.eg.db, strsplit(summary(gse)[,8],"/")[[i]], 'SYMBOL', 'ENTREZID')),collapse="/")}
  print(tagList(DT::datatable(cbind(summary(gse), gene_names), extensions = 'Buttons',options = list(dom = 'Blfrtip',buttons = c('excel', "csv"), fixedHeader = TRUE))))
  #write.csv(cbind(summary(gse), gene_names), "/mnt/DOSI/BNSRLAB/BIOINFO/Project/CRISPR_Melanie/03_StoryTelling/GO/KEGG/34_KEGG_Table.csv")
}
cat("\n \n")

cat("### Pathway Analysis (GO) {.tabset .tabset-fade} \n \n")
gse = enrichGO(entrez_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont="BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.2, readable = T) #universe=universe_blue_entrez, 

# Plot
cat("#### DotPlot  {.tabset .tabset-fade} \n \n")
print(enrichplot::dotplot(gse, showCategory=10,font.size=18)) #, split=".sign" + facet_grid(.~.sign)
cat("\n \n")

cat("#### cNET Plot  {.tabset .tabset-fade} \n \n")
edox <- setReadable(gse, 'org.Hs.eg.db', 'ENTREZID')
print(cnetplot(edox, showCategory=10,font.size=18))
cat("\n \n")

# Add length of pathways and GeneRatio
PathwayLength = rep(0, nrow(summary(gse)))
PercentPathwayMatch = rep(0, nrow(summary(gse)))
for(i in 1:nrow(gse)){
  PathwayLength[i] = length(gse@geneSets[[rownames(as.data.frame(summary(gse)))[i]]])
  PercentPathwayMatch[i] = paste0(round((summary(gse)[i,"Count"] / PathwayLength[i])*100,1), " %")}

cat("#### Table  {.tabset .tabset-fade} \n \n")
print(tagList(DT::datatable(cbind(summary(gse), PathwayLength, PercentPathwayMatch), extensions = 'Buttons',options = list(dom = 'Blfrtip',buttons = c('excel', "csv"), fixedHeader = TRUE))))
cat("\n \n")

}