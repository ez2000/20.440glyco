#####################
# CLUSTER ANALAYSIS #
#####################

#using Deep Learning paper github by Bojar et. al. 
# https://github.com/BojarLab/scGlycomics_b16_branching
#-----------------------------------------------------------------
# 0) Set up
#-----------------------------------------------------------------

# set the file path to folder containing baseline source script to load data and packages
path_to_init <- "/Users/erikazhang/Dropbox (MIT)/20.440 Biological Networks/project/source code"
# set the working directory to the file path for
setwd(path_to_init)
# source 
source("baseline source script.R") #can't run while open 

##-----------------------------------------------------------------
# 00) Remove mitochondrial and RNA counts below 200 and over 4000 --> ALREADY DONE! code left here just to show how 
#-----------------------------------------------------------------
# Make new column in gene expression matrix with percent mitochondrial, hemoglobin RNA 
#geneexp_seurat <- PercentageFeatureSet(geneexp_seurat, pattern = "^MT-", col.name = 'pct_mito')

# Filter out samples with QC characteristics. numbers taken from nbisweden.github.io doc
# RNA counts more than 200 is standard, less than 4000 is parameter given in paper
#selected_cells <- WhichCells(geneexp_seurat, expression = nFeature_RNA > 200 & nFeature_RNA < 4000)
#selected_feats <- rownames(geneexp_seurat)[Matrix::rowSums(geneexp_seurat) > 3]

# Select subset of data to filter out cells
#geneexp_seurat.filt <- subset(geneexp_seurat, features = selected_feats, cells = selected_cells)


##-----------------------------------------------------------------
# 1) QC: Remove Doublets & Demultiplex HTO (hashtag object) data (based on Bojat Lab GitHub)
#-----------------------------------------------------------------

# Normalize HTO data
sugarseq_seurat <- NormalizeData(sugarseq_seurat, assay = "HTO", normalization.method = "CLR", margin = 1)

# Demultiplex by hashtags
cluster_func <- "clara"
positive.quantile <- 0.8
sugarseq_seurat <- HTODemux(sugarseq_seurat, assay = "HTO", kfunc = cluster_func, positive.quantile = positive.quantile, verbose = TRUE)

# Plot results
HTOHeatmap(sugarseq_seurat, assay = "HTO", ncells = 5000)
RidgePlot(sugarseq_seurat, assay = "HTO", features = rownames(sugarseq_seurat[["HTO"]])[1:2], ncol = 2)

# Numbers of singlet/multiplet/negative cells
sugarseq_seurat.class <- as.data.frame(table(sugarseq_seurat$HTO_classification.global))

# Remove negatives and doublets
sugarseq_seurat <- subset(x = sugarseq_seurat, subset = HTO_classification.global == 'Singlet')
sugarseq_seurat.data.summary[["Number of all genes"]] <- dim(sugarseq_seurat@assays$RNA)[1]
sugarseq_seurat.data.summary[["Number of all cells"]] <- dim(sugarseq_seurat@assays$RNA)[2]
sugarseq_seurat.data.summary[c("HTODemux doublet", "HTODemux negative", "HTODemux singlet")] <- sc.class[, 2]
sugarseq_seurat.data.summary[, 1] <- NULL


#-----------------------------------------------------------------
# 2) SCTransform RNA Data
#-----------------------------------------------------------------
var_genes <- VariableFeatures(sugarseq_seurat)
sugarseq_seurat.data.summary[["Number of variable genes"]] <- length(var_genes)

# Compute percentage contribution of mitochondrial and ribosomal genes
sugarseq_seurat <- PercentageFeatureSet(sugarseq_seurat, pattern = "^mt-", col.name = "percent.mt")
sugarseq_seurat <- PercentageFeatureSet(sugarseq_seurat, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa", col.name = "percent.ribo")

# Compute cell cycle scores
cell.cycle.genes <- read.csv(file = paste0(datadir, cellcycle.file), header = TRUE, 
                             as.is = TRUE, stringsAsFactors = FALSE)
s.genes <- cell.cycle.genes[cell.cycle.genes$phase == "S", match("genename", colnames(cell.cycle.genes))]
g2m.genes <- cell.cycle.genes[cell.cycle.genes$phase == "G2/M", match("genename", colnames(cell.cycle.genes))]
sc <- CellCycleScoring(sc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# SCTransform data
regress.out.factors <- c("percent.mt", "percent.ribo", "S.Score", "G2M.Score")
sc <- SCTransform(sc, assay = 'RNA', new.assay.name = 'SCT',
                  vars.to.regress = regress.out.factors,
                  variable.features.n = NULL, variable.features.rv.th = 1.3, 
                  method = "glmGamPoi", verbose = TRUE)

sugarseq_seurat <- SCTransform(sugarseq_seurat, assay = 'RNA', new.assay.name = 'SCT',
                  vars.to.regress = regress.out.factors,
                  variable.features.n = NULL, variable.features.rv.th = 1.3, 
                  method = "glmGamPoi", verbose = TRUE)


#-----------------------------------------------------------------
# ALT 2) SCTransform RNA Data (my way, manually)
#-----------------------------------------------------------------
# Set default assay to RNA to work with gene expression
DefaultAssay(sugarseq_seurat) <- 'RNA'

#SCTransform replaces Normalize, find variable features, and ScaleData
sugarseq_seurat <- SCTransform(sugarseq_seurat, verbose = FALSE)

#select cells with residual variance >1.3 from SCTransform
DefaultAssay(sugarseq_seurat) <- 'SCT'
residual_variance <- rowVars(GetAssayData(object = sugarseq_seurat, slot= "scale.data")) 
selected_cells <- names(residual_variance[residual_variance > 1.3])
sugarseq_seurat <- subset(sugarseq_seurat, cells = selected_cells)

sugarseq_seurat <-RunPCA(sugarseq_seurat, verbose = FALSE)
sugarseq_seurat <- FindNeighbors(sugarseq_seurat, dims = 1:10,   k.param = 50, annoy.metric = 'cosine', verbose = FALSE)
sugarseq_seurat <-FindClusters(sugarseq_seurat, resolution=0.6)
sugarseq_seurat <-RunUMAP(sugarseq_seurat, dims = 1:10, n.neighbors = 50, metric = "cosine", verbose = FALSE)

DimPlot(sugarseq_seurat)


final.matrix.rna <- as.data.frame(t(as.data.frame(GetAssayData(sugarseq_seurat), 
                                                  assay = "SCT", slot = "data")[var_genes,]))

##----------------------------------------------------------------------------------------
# 3) Identify differentially expressed genes (features) per cluster and make feature plots
#------------------------------------------------------------------------------------------

ss_seurat <- FindAllMarkers(seq_seurat)
head(ss_seurat)


#-----------------------------------------------------------------
# 4) Identify T cell subtype using ProjecTIL  
#-----------------------------------------------------------------

# Load reference map
ref.map <- load.reference.map()
# Run projection
query.projected <- make.projection(ss_seurat, ref=ref.map, query.assay = "RNA")
plot.projection(ref.map, query.projected)
# Classify cells
query.projected <- cellstate.predict(ref=ref.map, query=query.projected)
plot.statepred.composition(ref.map, query.projected, metric = "Percent")
# Export cell states
fun.cluster <- data.frame(names(Idents(query.projected)))
colnames(fun.cluster) <- c("Cell")
fun.cluster["Functional Cluster"] <- query.projected@meta.data[["functional.cluster"]]
fun.cluster[[lectin.name]] <- as.data.frame(t(query.projected@assays$ADT@data))[[lectin.col.name]]
```


## Remove non-T cells from the final matrix identified by scGate

```{r}
filter.1 <- rownames(final.matrix.rna)%in%(fun.cluster$Cell)
sum(as.integer(filter.1))
final.matrix.rna <- final.matrix.rna[filter.1, ]
final.matrix.adt <- final.matrix.adt[filter.1, ]
final.matrix.adt["Type"] <- fun.cluster$`Functional Cluster`
final.matrix.adt <- final.matrix.adt[, c(3,4)]
colnames(final.matrix.adt)[1] <- lectin.name
sc.data.summary[["Number of non-T cells"]] <- sum(as.integer(!filter.1))
sc.data.summary[["Final number of cells"]] <- sum(as.integer(filter.1))
```


## Write to file

```{r}
if (wfile){
  write.table(final.matrix.rna,
              file = paste0(datadir, sample.set.name, '_transformed_data.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(final.matrix.adt,
              file = paste0(datadir, sample.set.name, '_transformed_identity.csv'), 
              sep = ',', row.names = T, col.names = NA, quote = F)
  
  write.table(sc.data.summary, file = paste0(datadir, sample.set.name, '_data_summary.csv'),
              sep = ',', row.names = F, col.names = T, quote = F)
}
```


## # Replace cell Idents

```{r}
new.ident <- query.projected@meta.data[["functional.cluster"]]
names(new.ident) <- names(Idents(query.projected))
Idents(query.projected) <- new.ident
