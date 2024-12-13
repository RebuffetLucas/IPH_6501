#Try solution for the enrichKEGG function that does not work

#1- Update the clusterprofilier to the latest Github version ( the lastest version is 4.7.1.3)
remotes::install_github("YuLab-SMU/clusterProfiler") 
#BiocManager::install("DOSE", version = "3.3")
BiocManager::install("clusterProfiler")
# 2- Establish a local KEGG database
# install the packages
remotes::install_github("YuLab-SMU/createKEGGdb")
# import the library and create a KEGG database locally 
library(createKEGGdb)
species <-c("ath","hsa","mmu", "rno","dre","dme","cel")
createKEGGdb::create_kegg_db(species = species)
# You will get KEGG.db_1.0.tar.gz file in your working directory 
# 3- install the KEGG.db and import it

install.packages("KEGG.db_1.0.tar.gz", repos=NULL,type="source")
library(KEGG.db)

#add use_internal_data=T in your enrichKEGG function

data(gcSample)
yy = enrichKEGG(gcSample[[5]], pvalueCutoff=0.01, use_internal_data=T)
head(summary(yy))