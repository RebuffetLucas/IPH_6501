## @knitr Visualize_Enrichment

#Simple vizualisation
#Not used in the pipeline

#Table of extraction
#For TableS3 page 1

lanes_of_interest= c(1,2,8,11,13,19,26,33,42,48,66,76,82) #Lanes of interest of unique(top_descriptions$Description)
Description_of_interest  = unique(top_descriptions$Description)[lanes_of_interest]





#  !!!!!!!    CODE à adapter pour eviter le soucis avec les descrptions redondantes !!!!!!!! #



dfm  = top_descriptions[top_descriptions$Description %in%Description_of_interest,]


#Modify the legend and values for visualization
dfm <- dfm %>%
  dplyr::mutate(cluster = ifelse(grepl("T5", cluster), "UP", ifelse(grepl("US", cluster), "DOWN", cluster)))

dfm <- dfm %>%
  mutate(Generatio_bis = ifelse(cluster == "UP", GeneRatio, -GeneRatio))

dfm <- dfm %>%
  mutate(neg_log_pvalue = -log(pvalue, base = 10))



#Vizualisation
ggbarplot(dfm, x = "Description", y = "GeneRatio",
          fill = "Cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          #palette = PALETTE_UP_DOWN,            # jco journal color palett. see ?ggpar
          sort.val = "desc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15))



ggbarplot(dfm, x = "Description", y = "GeneRatio",
          fill = "Cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          #palette = "jco",            # jco journal color palett. see ?ggpar
          sort.val = "desc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15)) + coord_flip()


ggbarplot(dfm, x = "Description", y = "Generatio_bis",
          fill = "cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          palette = PALETTE_UP_DOWN,            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = FALSE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15)) + coord_flip()


ggbarplot(dfm, x = "Description", y = "neg_log_pvalue",
          fill = "cluster",               # change fill color by cyl
          color = "Black",            # Set bar border colors to black
          palette = PALETTE_UP_DOWN,            # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # Sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal(),
          width= 0.7
) +  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + theme(text = element_text(size = 15)) + coord_flip()




