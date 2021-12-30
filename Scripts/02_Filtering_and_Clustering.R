#R script created on November, 5th, 2021
# This script takes the seurat object created on the previous step, reduce it for the desired stages and run the basic seurat pipeline. 
# As an output you will have a new Seurat object as well as a single cell object.
# Work from your 'Script' dir

#Load libraries
library(Seurat)
library(tidyverse)

suppressMessages(require(useful))
suppressMessages(require(scater))
suppressMessages(require(dplyr))
suppressMessages(require(cowplot))

options(stringsAsFactors = FALSE)

## 1. Load Seurat object -------------------------------------------------------
seurat_object <- readRDS("../E-MTAB-6946_10x/Raw/prep_10x/SeuratObj_red")
#It is already filtered as in the original paper
summary(seurat_object@meta.data$nFeature_RNA)

#reduce for meiosis stages P10, 15, 20
stages = c("do17821", "do18195", "do17824")

seurat_object <- subset(x = seurat_object, subset = Library %in%  stages)
rm(stages)

# 3. Normalizing the data ------------------------------------------------------
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

# 4. Higly variable genes ------------------------------------------------------
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

# 5. Scaling the data and removing unwanted sources of variation----------------
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

#Run PCA
seurat_object <- RunPCA(seurat_object, pc.genes = seurat_object@var.genes, npcs = 40, verbose = F) 
rm(all.genes)

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = seurat_object, reduction = "pca", pt.size = .1, group.by = "Library"
              )
p2 <- VlnPlot(object = seurat_object, features = "PC_1", group.by = "Library", 
              pt.size = .1)
p3 <- VlnPlot(object = seurat_object, features = "PC_2", group.by = "Library", 
              pt.size = .1)
p4 <- VlnPlot(object = seurat_object, features = "PC_3", group.by = "Library", 
              pt.size = .1)
plot_grid(p1,p2,p3,p4)


E1 <- ElbowPlot(object = seurat_object,  ndims = 40)

# Determine percent of variation associated with each PC
# https://hbctraining.github.io/scRNA-seq_online/lessons/elbow_plot_metric.html
pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
co2

# Minimum of the two calculation
pcs <- min(co1, co2)

pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))

# Elbow plot to visualize 
E2 <-  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
 
plot_grid(E1, E2)  



  ######################

seurat_object <- seurat_object %>% 
    RunUMAP(reduction = "pca", dims = 1:16) %>% 
    RunTSNE(reduction = "pca", dims = 1:16) %>% 
    FindNeighbors(reduction = "pca", dims = 1:16) %>% 
  # Determine the clusters for various resolutions                                
    FindClusters(resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)) %>% 
    identity()

#Save seurat object
seurat_object <- readRDS("../E-MTAB-6946_10x/Raw/prep_10x/SeuratObj_red")
dir.create("../E-MTAB-6946_10x/Processed")

saveRDS(seurat_object, file = "../E-MTAB-6946_10x/Processed/Meiosis_SeuratObj.rds")

#Save expression as a single cell object
my_order = c("Spermatogonia", "eP1", "eP2", "mP", "lP1", "lP2", "D", "MI", "MII", "S1", "S2", "Fetal_Leydig_1", "Fetal_Leydig_2", "Leydig_1", "Leydig_2", "Sertoli", "Endothelial_cells", "PTM", "Interstitial_tMg")
#remove outliers, S3 and S4

seurat_object <- subset(x = seurat_object, subset = AnnotatedClusters %in%  my_order)

sce_object <- as.SingleCellExperiment(seurat_object)
saveRDS(sce_object, file = "../E-MTAB-6946_10x/Processed/Meiosis_sce.rds")


