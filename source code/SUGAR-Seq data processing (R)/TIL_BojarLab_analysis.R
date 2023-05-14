#####################################
### TIL clustering & data analysis ###
#####################################

#----------------------------------------------------------------------------------------
# 0) Install necessary packages and set file path/working directory to raw data directory
#----------------------------------------------------------------------------------------

# Run set up script to install key packages and set working directory
path_to_package <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/source code/setup script.R"
source(path_to_package)

# set the working directory to the file path
path_to_raw <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/raw data/TIL data"

#set working directory to output folder
path_to_output <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/TIL ouput/"
setwd(path_to_output)

# list files in the directory to check if the file is present
list.files()

### SET PARAMETERS ### 
wfile <- TRUE
# The number of decimal places to save in the final .csv table (may help reduce file size)
n.decimal <- 3 # n.decimal <- FALSE, if rounding up is not desired
# Name of the column containing the lectin data
lectin.col.name <- "Biotin"
# Lectin name, for output
lectin.name <- "PHA-L"

#sample set name
sample.set.name <- "TIL"

# Create a dataframe for data summary
ss_seuratTIL.data.summary <- as.data.frame(matrix(nrow = 1, ncol = 1))

#-----------------------------------------------------------------
# 1) Initialize data and create Seurat objects
#-----------------------------------------------------------------
# Use Read10X() to extract data
ss.data <- Read10X(path_to_raw)

# Make Seurat object with info from the three matrices with assays: "RNA", "ADT", "CUSTOM"
TILss_seurat <- CreateSeuratObject(counts = ss.data[['Gene Expression']], project = 'sugarseq_seurat', min.cells = 0) #default assay is RNA
TILss_seurat[['ADT']] <- CreateAssayObject(counts = ss.data[['Antibody Capture']])
TILss_seurat[['HTO']] <- CreateAssayObject(counts = ss.data[['Custom']])

DefaultAssay(TILss_seurat) <- 'RNA'

#export each assay as csv file
TIL_RNA <- GetAssayData(object = TILss_seurat, assay = "RNA", slot = "data")
write.csv(TIL_RNA, "TIL RNA matrix.csv")

TIL_ADT <- GetAssayData(object = TILss_seurat, assay = "ADT", slot = "data")
write.csv(TIL_ADT, "TIL ADT matrix.csv")

#add cell information from SUGAR-Seq seurat object to dataframe
ss_seuratTIL.data.summary[["Number of all genes"]] <- dim(TILss_seurat@assays$RNA)[1]
ss_seuratTIL.data.summary[["Number of all cells"]] <- dim(TILss_seurat@assays$RNA)[2]
ss_seuratTIL.data.summary[, 1] <- NULL

#-----------------------------------------------------------------
# 2) Run SCTransform and PCA
#-----------------------------------------------------------------
var_genes <- VariableFeatures(TILss_seurat)
TILss_seurat <- SCTransform(TILss_seurat, assay = 'RNA', new.assay.name = 'SCT',
                               variable.features.n = NULL, variable.features.rv.th = 1.3, 
                               method = "glmGamPoi", verbose = TRUE)

#save csv of normalized data for
TIL_SCT <- GetAssayData(object = TILss_seurat, assay = "SCT", slot = "data")
write.csv(TIL_SCT, "TIL SCT matrix.csv")

TILss_seurat <-RunPCA(TILss_seurat, verbose = FALSE)
TILss_seurat <- FindNeighbors(TILss_seurat, dims = 1:10,   k.param = 50, annoy.metric = 'cosine', verbose = FALSE)
TILss_seurat <-FindClusters(TILss_seurat, resolution=0.6)
TILss_seurat <-RunUMAP(TILss_seurat, dims = 1:10, n.neighbors = 50, metric = "cosine", verbose = FALSE)
DimPlot(TILss_seurat)


#-----------------------------------------------------------------
# 3) Normalize lectin data and write to file, keeping only variable genes
#-----------------------------------------------------------------

ss_seuratTIL.data.summary[["Number of variable genes"]] <- length(var_genes)

TILss_seurat <- NormalizeData(TILss_seurat, assay = 'ADT', normalization.method = "CLR", margin = 1)
final.matrix.rna <- as.data.frame(t(as.data.frame(GetAssayData(TILss_seurat), 
                                                  assay = "SCT", slot = "data")[var_genes,]))
final.matrix.adt <- as.data.frame(t(TILss_seurat@assays$ADT@data))
final.matrix.rna[lectin.name] <- as.data.frame(final.matrix.adt[lectin.col.name][[1]])

if (n.decimal != 0){
  final.matrix.rna <- round(final.matrix.rna, digits = n.decimal)}


#-----------------------------------------------------------------
# 4) Identify T cell subtype using ProjecTIL  
#-----------------------------------------------------------------

# Load reference map
ref.map <- readRDS("/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/TIL data/ref_TILAtlas_mouse_v1.rds")

# Run projection
query.projected <- make.projection(TILss_seurat, ref=ref.map, query.assay = "RNA")
plot.projection(ref.map, query.projected, linesize = 0, ref.size = 0, pointsize = 0)+ ggtitle("LN Make.Projection ver.")
DimPlot(query.projected)


# Classify cells
query.projected <- cellstate.predict(ref=ref.map, query=query.projected)
plot.statepred.composition(ref.map, query.projected, metric = "Percent")

# Export cell states
fun.cluster <- data.frame(names(Idents(query.projected)))
colnames(fun.cluster) <- c("Cell")
fun.cluster["Functional Cluster"] <- query.projected@meta.data[["functional.cluster"]]
fun.cluster[[lectin.name]] <- as.data.frame(t(query.projected@assays$ADT@data))[[lectin.col.name]]


## Remove non-T cells from the final matrix identified by scGate
filter.1 <- rownames(final.matrix.rna)%in%(fun.cluster$Cell)
sum(as.integer(filter.1))

final.matrix.rna <- final.matrix.rna[filter.1, ]
final.matrix.adt <- final.matrix.adt[filter.1, ]
final.matrix.adt["Type"] <- fun.cluster$`Functional Cluster`
final.matrix.adt <- final.matrix.adt[, c(3,4)]
colnames(final.matrix.adt)[1] <- lectin.name
ss_seuratTIL.data.summary[["Number of non-T cells"]] <- sum(as.integer(!filter.1))
ss_seuratTIL.data.summary[["Final number of cells"]] <- sum(as.integer(filter.1))


## Write to file
if (wfile){
  write.table(final.matrix.rna,
              file = paste0(path_to_output, sample.set.name, '_transformed_data.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(final.matrix.adt,
              file = paste0(path_to_output, sample.set.name, '_transformed_identity.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(ss_seuratTIL.data.summary, file = paste0(path_to_output, sample.set.name, '_data_summary.csv'),
              sep = ',', row.names = F, col.names = T, quote = F)
}

## # Replace cell Idents
new.ident <- query.projected@meta.data[["functional.cluster"]]
names(new.ident) <- names(Idents(query.projected))
Idents(query.projected) <- new.ident


#-------------------
##current attempt###
#-------------------

# filter TILss_seurat based on glycogenes and do findallmarkers
glycogenes <- read.csv("Mouse glycogenes.csv")
glycoTILss <- query.projected[,glycogenes]

FindAllMarkers(TILss_seurat)

