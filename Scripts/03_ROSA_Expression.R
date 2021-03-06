#R script created on November, 5th, 2021
# This script takes the single cell and seurat object created on the previous step and produces plots for some marker genes defined in the original paper and the ROSA 26 locus
# Work from your 'Script' dir

#Load libraries
library(tidyverse)
library(useful)
library(cowplot)
library(ggpubr)
library(scater)

options(stringsAsFactors = FALSE)

## 1. Load Single cell experiment object ---------------------------------------
sce_object <- readRDS("../E-MTAB-6946_10x/Processed/Meiosis_sce.rds")

## 2. Gene expression level per cell type --------------------------------------
## 2.1. Organize
my_order = c("Spermatogonia", "eP1", "eP2", "mP", "lP1", "lP2", "D", "MI", "MII", "S1", "S2", "Fetal_Leydig_1", "Fetal_Leydig_2", "Leydig_1", "Leydig_2", "Sertoli", "Endothelial_cells", "PTM", "Interstitial_tMg")

sce_object$AnnotatedClusters <- factor(sce_object$AnnotatedClusters, 
                          levels = my_order)

##2.2 choose your genes
Gene.df = as.data.frame(rownames(sce_object))
View(Gene.df)

##Fig5 gene expression during male meiosis
genes = c("ENSMUSG00000054717 Hmgb2", "ENSMUSG00000027855 Sycp1", "ENSMUSG00000062248 Cks2", "ENSMUSG00000016559 H3f3b", "ENSMUSG00000060445 Sycp2", "ENSMUSG00000029423 Piwil1", "ENSMUSG00000018554 Ybx2", "ENSMUSG00000093668 Pou5f2", "ENSMUSG00000019942 Cdk1", "ENSMUSG00000027496 Aurka", "ENSMUSG00000005233 Spc25", "ENSMUSG00000029580 Actb")


plotExpression(sce_object, 
               features = genes , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 0.6, #tip: if you want to hide the dots you can set point_alpha = 0 
               colour_by = "Sample",
               ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 #
genes = c(
  "ENSMUSG00000030264 Thumpd3",
  "ENSMUSG00000034269 Setd5",
  "ENSMUSG00000086429 Gt(ROSA)26Sor",
  "ENSMUSG00000029580 Actb"
)

plotExpression(sce_object, 
               features = genes , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 0.6, #tip: if you want to hide the dots you can set point_alpha = 0 
               colour_by = "BroadClusters",
               ncol = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


genes = c(
  "ENSMUSG00000027323 Rad51",
  "ENSMUSG00000022429 Dmc1",
  "ENSMUSG00000002076 Hsf2bp",
  "ENSMUSG00000000751 Rpa1",
  "ENSMUSG00000020059 Sycp3",
  "ENSMUSG00000086429 Gt(ROSA)26Sor"
)

plotExpression(sce_object, 
               features = genes , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 0.6, #tip: if you want to hide the dots you can set point_alpha = 0 
               colour_by = "BroadClusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plotExpression(sce_object, 
               features = "ENSMUSG00000086429 Gt(ROSA)26Sor" , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 0.6, #tip: if you want to hide the dots you can set point_alpha = 0 
               colour_by = "BroadClusters") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#
library(Seurat)
ROSA = sce_object[("ENSMUSG00000086429 Gt(ROSA)26Sor"), ]

ROSA.seurat <- CreateSeuratObject(counts = counts(ROSA), meta.data = as.data.frame(colData(ROSA)))
ROSA.seurat <- SetAssayData(object = ROSA.seurat, slot = "data", new.data = logcounts(ROSA))

ROSA.P10.seurat <- subset(x = ROSA.seurat, subset = BroadClusters == "Germ")
ROSA.P10.seurat <- subset(x = ROSA.P10.seurat, subset = Sample == "P10")

ROSA.P15.seurat <- subset(x = ROSA.seurat, subset = BroadClusters == "Germ")
ROSA.P15.seurat <- subset(x = ROSA.P15.seurat, subset = Sample == "P15")

ROSA.P20.seurat <- subset(x = ROSA.seurat, subset = BroadClusters == "Germ")
ROSA.P20.seurat <- subset(x = ROSA.P20.seurat, subset = Sample == "P20")

ROSA.P10 = as.SingleCellExperiment(ROSA.P10.seurat)
ROSA.P15 = as.SingleCellExperiment(ROSA.P15.seurat)
ROSA.P20 = as.SingleCellExperiment(ROSA.P20.seurat)

p10 <- plotExpression(ROSA.P10, 
               features = "ENSMUSG00000086429 Gt(ROSA)26Sor" , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 1, 
               colour_by = "AnnotatedClusters"
               ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")

p15 <- plotExpression(ROSA.P15, 
               features = "ENSMUSG00000086429 Gt(ROSA)26Sor" , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 1, 
               colour_by = "AnnotatedClusters"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")


p20 <- plotExpression(ROSA.P20, 
               features = "ENSMUSG00000086429 Gt(ROSA)26Sor" , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 1, 
               colour_by = "AnnotatedClusters"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")

ggarrange(p10,p15,p20, ncol = 3)


genes = c(
  "ENSMUSG00000030264 Thumpd3",
  "ENSMUSG00000034269 Setd5",
  "ENSMUSG00000086429 Gt(ROSA)26Sor")

ROSA.locus = sce_object[genes, ]

ROSA.locus.seurat <- CreateSeuratObject(counts = counts(ROSA.locus), meta.data = as.data.frame(colData(ROSA.locus)))
ROSA.locus.seurat <- SetAssayData(object = ROSA.locus.seurat, slot = "data", new.data = logcounts(ROSA.locus))


ROSA.locus.P10.seurat <- subset(x = ROSA.locus.seurat, subset = BroadClusters == "Germ")
ROSA.locus.P10.seurat <- subset(x = ROSA.locus.P10.seurat, subset = Sample == "P10")

ROSA.locus.P15.seurat <- subset(x = ROSA.locus.seurat, subset = BroadClusters == "Germ")
ROSA.locus.P15.seurat <- subset(x = ROSA.locus.P15.seurat, subset = Sample == "P15")

ROSA.locus.P20.seurat <- subset(x = ROSA.locus.seurat, subset = BroadClusters == "Germ")
ROSA.locus.P20.seurat <- subset(x = ROSA.locus.P20.seurat, subset = Sample == "P20")

ROSA.locus.P10 = as.SingleCellExperiment(ROSA.locus.P10.seurat)
ROSA.locus.P15 = as.SingleCellExperiment(ROSA.locus.P15.seurat)
ROSA.locus.P20 = as.SingleCellExperiment(ROSA.locus.P20.seurat)


p10 <- plotExpression(ROSA.locus.P10, 
                      features = genes , #define your genes
                      x = "AnnotatedClusters", #group by cluster identity
                      point_alpha = 1, 
                      colour_by = "AnnotatedClusters",
               ncol = 1
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")+ ggtitle("P10")

p15 <- plotExpression(ROSA.locus.P15, 
                      features = genes , #define your genes
                      x = "AnnotatedClusters", #group by cluster identity
                      point_alpha = 1, 
                      colour_by = "AnnotatedClusters",
               ncol = 1
               ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")+ ggtitle("P15")


p20 <- plotExpression(ROSA.locus.P20, 
                      features = genes,
                      x = "AnnotatedClusters", #group by cluster identity
                      point_alpha = 1, 
                      colour_by = "AnnotatedClusters",
                      ncol = 1
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")+ ggtitle("P20")

ggarrange(p10,p15,p20, ncol = 3)


#(...)
