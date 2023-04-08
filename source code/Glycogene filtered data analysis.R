############################################
##### GLYCOGENE FILTERED DATA ANALYSIS #####
############################################


'''
In this script, we are going to make seurat objects of the matrices with the cells filtered
for those that are T-cells, and genes filtered for glycogenes only.
We will perform a FindAllMarkers of DESeq analysis to look at differential expression of
the glycogenes in these sets.
'''

#----------------------------------------------------------------------------------------
# 0) Install necessary packages and set file path/working directory to raw data directory
#----------------------------------------------------------------------------------------

# Run set up script to install key packages and set working directory
path_to_package <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/source code/setup script.R"
source(path_to_package)

# set the working directory to the file path
path_to_TIL <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/TIL output/TIL glycofiltered"
path_to_LN <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/LN output/LN glycofiltered"

#set working directory to output folder
path_to_output <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/"
setwd(path_to_output)

#-----------------------------------------------------------------
# 1) Initialize data & create Seurat object of glycofiltered cells
#-----------------------------------------------------------------
path_to_raw <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/raw data/TIL data/"

glycoTIL.data <- Read10X(path_to_raw)
glycoTIL <- CreateSeuratObject(counts = glycoTIL.data[['Gene Expression']], min.cells = 0) #default assay is RNA

#-------------------
##current attempt###
#-------------------
# 2) DESeq on raw non-normalized data
library(DESeq2)
glycoTIL <- DESeq(glycoTIL)
TIL_res <- results(glycoTIL)


##do find all markers instead , compare to glycogenes 
glycoTIL.markers <- FindAllMarkers(glycoTIL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'DESeq2')
glycoTIL.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

