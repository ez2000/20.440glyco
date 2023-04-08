#####################################
### LN clustering & data analysis ###
#####################################

#----------------------------------------------------------------------------------------
# 0) Install necessary packages and set file path/working directory to raw data directory
#----------------------------------------------------------------------------------------

# Run set up script to install key packages and set working directory
path_to_package <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/source code/setup script.R"
source(path_to_package)

path_to_raw <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/raw data/LN data"

# set the working directory to the output
path_to_output <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/LN output/"
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
sample.set.name <- "LN"

# Create a dataframe for data summary
ss_seuratLN.data.summary <- as.data.frame(matrix(nrow = 1, ncol = 1))

#-----------------------------------------------------------------
# 1) Initialize data and create Seurat objects
#-----------------------------------------------------------------
# Use Read10X() to extract data
LNss_seurat.data <- Read10X(path_to_raw)

# Make Seurat object with info from the three matrices with assays: "RNA", "ADT", "CUSTOM"
LNss_seurat <- CreateSeuratObject(counts = LNss_seurat.data[['Gene Expression']], project = 'sugarseq_seurat', min.cells = 0) #default assay is RNA
LNss_seurat[['ADT']] <- CreateAssayObject(counts = LNss_seurat.data[['Antibody Capture']])
LNss_seurat[['HTO']] <- CreateAssayObject(counts = LNss_seurat.data[['Custom']])

DefaultAssay(LNss_seurat) <- 'RNA'

#export each assay as csv file
LN_RNA <- GetAssayData(object = LNss_seurat, assay = "RNA", slot = "data")
write.csv(TIL_RNA, "LN RNA matrix.csv")

LN_ADT <- GetAssayData(object = LNss_seurat, assay = "ADT", slot = "data")
write.csv(TIL_ADT, "LN ADT matrix.csv")


#add cell information from SUGAR-Seq seurat object to dataframe
ss_seuratLN.data.summary[["Number of all genes"]] <- dim(ss_seurat@assays$RNA)[1]
ss_seuratLN.data.summary[["Number of all cells"]] <- dim(ss_seurat@assays$RNA)[2]
ss_seuratLN.data.summary[, 1] <- NULL

#-----------------------------------------------------------------
# 2) Run SCTransform and PCA
#-----------------------------------------------------------------
var_genes <- VariableFeatures(ss_seurat)
LNss_seurat <- SCTransform(LNss_seurat, assay = 'RNA', new.assay.name = 'SCT',
                               variable.features.n = NULL, variable.features.rv.th = 1.3, 
                               method = "glmGamPoi", verbose = TRUE)

#save csv of normalized data for
LN_SCT <- GetAssayData(object = LNss_seurat, assay = "SCT", slot = "data")
write.csv(LN_SCT, "LN SCT matrix.csv")


LNss_seurat <-RunPCA(LNss_seurat, verbose = FALSE)
LNss_seurat <- FindNeighbors(LNss_seurat, dims = 1:10,   k.param = 50, annoy.metric = 'cosine', verbose = FALSE)
LNss_seurat <-FindClusters(LNss_seurat, resolution=0.6)
LNss_seurat <-RunUMAP(LNss_seurat, dims = 1:10, n.neighbors = 50, metric = "cosine", verbose = FALSE)
DimPlot(LNss_seurat)

#-----------------------------------------------------------------
# 3) Normalize lectin data and write to file, keeping only variable genes
#-----------------------------------------------------------------

ss_seuratLN.data.summary[["Number of variable genes"]] <- length(var_genes)

LNss_seurat <- NormalizeData(LNss_seurat, assay = 'ADT', normalization.method = "CLR", margin = 1)
final.matrix.rna <- as.data.frame(t(as.data.frame(GetAssayData(LNss_seurat), 
                                                  assay = "SCT", slot = "data")[var_genes,]))
final.matrix.adt <- as.data.frame(t(LNss_seurat@assays$ADT@data))
final.matrix.rna[lectin.name] <- as.data.frame(final.matrix.adt[lectin.col.name][[1]])

if (n.decimal != 0){
  final.matrix.rna <- round(final.matrix.rna, digits = n.decimal)}


#-----------------------------------------------------------------
# 4) Identify T cell subtype using ProjecTIL  
#-----------------------------------------------------------------

# Load reference map
ref.map <- readRDS("/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/TIL data/ref_TILAtlas_mouse_v1.rds")

# Run projection
query.projected <- make.projection(ss_seurat, ref=ref.map) 
plot.projection(ref.map, query.projected, linesize = 0, ref.size = 0, pointsize = 0)+ ggtitle("LN Make.Projection ver.")
DimPlot(query.projected)

# Classify cells
query.projected <- cellstate.predict(ref=ref.map, query=query.projected)
plot.statepred.composition(ref.map, query.projected, metric = "Percent") + ggtitle("LN Cell state composition")

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
ss_seuratLN.data.summary[["Number of non-T cells"]] <- sum(as.integer(!filter.1))
ss_seuratLN.data.summary[["Final number of cells"]] <- sum(as.integer(filter.1))

## Write to file, export to desktop

if (wfile){
  write.table(final.matrix.rna,
              file = paste0(path_to_output, sample.set.name, '_transformed_data.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(final.matrix.adt,
              file = paste0(path_to_output, sample.set.name, '_transformed_identity.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(ss_seuratLN.data.summary, file = paste0(path_to_output, sample.set.name, '_data_summary.csv'),
              sep = ',', row.names = F, col.names = T, quote = F)
}



## # Replace cell Idents
new.ident <- query.projected@meta.data[["functional.cluster"]]
names(new.ident) <- names(Idents(query.projected))
Idents(query.projected) <- new.ident

