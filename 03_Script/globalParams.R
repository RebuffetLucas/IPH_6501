###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples and several analysis steps
###############################################################################


#### General
#### Variables that contains general information on the project and experiment

# A global description of the projet
# It will appear at the beginning of all analysis reports
# Example : "Study of Mouse Memory B Cell infected by Flu"
GLOBAL_DESCRIPTION = "This project aims at providing additional analysis for Demaria's paper to provide insights into how the different subsets of NK cells respond to IPH_6501"

# The name of you scientific group as it appear in the DOSI folder
# Example : "MGLAB"
SCIENTIFIC_GROUP = "EVLAB"

# The name of the scientific projet. 
# It must be the same name as the name of the project folder you defined
# Example : "moFluMemB"
SCIENTIFIC_PROJECT_NAME = "IPH_6501"

# The experiment name this file is in
# It must be the same name as the name of the experiment folder you defined
# Example : "10x_190712_m_moFluMemB"
EXPERIMENT_NAME = "Description_GO_KEGG_analysis"


#### Automatic definition of path
#### You may not modify the following section where standard path are defined in
#### constants so to be used in your code. Have a look at the declared variable
#### and use it in your code.

#### Input / Output

# This is the path of the global project folder
PATH_PROJECT = file.path( "/mnt/DOSI", 
                                        SCIENTIFIC_GROUP,
                                        "BIOINFO", 
                                        "BIOINFO_PROJECT",
                                        SCIENTIFIC_PROJECT_NAME)

# This is the path of the current experiment this file is in
PATH_EXPERIMENT = file.path( PATH_PROJECT)
#PATH_EXPERIMENT = file.path( PATH_PROJECT,EXPERIMENT_NAME)

# Those are the path to the main folder : RAW data, REFERENCE data and Output of analysis
PATH_EXPERIMENT_RAWDATA   = file.path( PATH_EXPERIMENT, "00_RawData")
PATH_EXPERIMENT_REFERENCE = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_REFERENCE_Data = file.path(PATH_EXPERIMENT_REFERENCE , "Data")
PATH_EXPERIMENT_OUTPUT    = file.path( PATH_EXPERIMENT, "05_Output")

PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity")

#Some more specific path:
#PATH_EXPERIMENT_OUTPUT_RegNet_BM_Figs = file.path(PATH_EXPERIMENT_OUTPUT ,  "/02_RegulatoryNetworkAnalysis/BM/Final_Figures")
#PATH_EXPERIMENT_OUTPUT_RegNet_FL_Figs  = file.path(PATH_EXPERIMENT_OUTPUT ,  "/02_RegulatoryNetworkAnalysis/FL/Final_Figures")
#PATH_EXPERIMENT_OUTPUT_RDSFiles = file.path(PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity/Seurat_Objects")
PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures = file.path(PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity/Figures")
PATH_EXPERIMENT_OUTPUT_GlobalHeteroDEG= file.path(PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity/DEG")

PATH_EXPERIMENT_OUTPUT_FIGURE_GlobalHeteroDEG = file.path(PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity/Figures/Genes_Correlations")
PATH_EXPERIMENT_OUTPUT_FIGURE_KEGG_DotPlot = file.path(PATH_EXPERIMENT_OUTPUT , "02_KEGG_GO/Figures/DotPlot")


PATH_EXPERIMENT_OUTPUT_GlobalHetero_FeaturePlot= file.path(PATH_EXPERIMENT_OUTPUT , "01_GlobalHeterogeneity/FeaturePlot")


#PATH_EXPERIMENT_OUTPUT_DESTINY_FIGURES = file.path(PATH_EXPERIMENT_OUTPUT , "03_Trancriptionnal_Trajectories/Destiny/Figures")
#PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES = file.path(PATH_EXPERIMENT_OUTPUT , "03_Trancriptionnal_Trajectories/Monocle/Figures")


# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME   )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



