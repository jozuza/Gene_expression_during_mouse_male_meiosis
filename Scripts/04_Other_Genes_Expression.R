library(tidyverse)
library(Seurat)
library(useful)
library(cowplot)
library(ggpubr)
library(scater)

setwd("~/Documents/Joost_Lab_2021/ROSA26_meiosis_expression/E-MTAB-6946_10x")

options(stringsAsFactors = FALSE)

## 1. Load Single cell experiment object ---------------------------------------
sce_object <- readRDS("Seurat_processed/Meiosis_sce.rds")

## 2. Gene expression level per cell type --------------------------------------
## 2.1. Organize
my_order = c("Spermatogonia", "eP1", "eP2", "mP", "lP1", "lP2", "D", "MI", "MII", "S1", "S2", "Fetal_Leydig_1", "Fetal_Leydig_2", "Leydig_1", "Leydig_2", "Sertoli", "Endothelial_cells", "PTM", "Interstitial_tMg")

sce_object$AnnotatedClusters <- factor(sce_object$AnnotatedClusters, 
                          levels = my_order)

##2.2 choose your genes
Gene.df = as.data.frame(rownames(sce_object))
# View(Gene.df)

genes = c("ENSMUSG00000025105 Bnc1",
          "ENSMUSG00000041147 Brca2",
          "ENSMUSG00000022429 Dmc1",
          # "ENSMUSG00000059970 Hspa2",
          # "ENSMUSG00000022556 Hsf1",
          # "ENSMUSG00000019878 Hsf2",
          # "ENSMUSG00000045802 Hsf3",
          # "ENSMUSG00000033249 Hsf4",
          # "ENSMUSG00000002076 Hsf2bp",
          # "ENSMUSG00000044702 Palb2",
          "ENSMUSG00000027323 Rad51")


plotExpression(sce_object, 
               features = genes , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 0.6, #tip: if you want to hide the dots you can set point_alpha = 0 
               colour_by = "Sample",
               ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#
library(Seurat)
Gen_Set = sce_object[genes, ]

Gen_Set.seurat <- CreateSeuratObject(counts = counts(Gen_Set), meta.data = as.data.frame(colData(Gen_Set)))
Gen_Set.seurat <- SetAssayData(object = Gen_Set.seurat, slot = "data", new.data = logcounts(Gen_Set))

Gen_Set.P10.seurat <- subset(x = Gen_Set.seurat, subset = BroadClusters == "Germ")
Gen_Set.P10.seurat <- subset(x = Gen_Set.P10.seurat, subset = Sample == "P10")

Gen_Set.P15.seurat <- subset(x = Gen_Set.seurat, subset = BroadClusters == "Germ")
Gen_Set.P15.seurat <- subset(x = Gen_Set.P15.seurat, subset = Sample == "P15")

Gen_Set.P20.seurat <- subset(x = Gen_Set.seurat, subset = BroadClusters == "Germ")
Gen_Set.P20.seurat <- subset(x = Gen_Set.P20.seurat, subset = Sample == "P20")

Gen_Set.P10 = as.SingleCellExperiment(Gen_Set.P10.seurat)
Gen_Set.P15 = as.SingleCellExperiment(Gen_Set.P15.seurat)
Gen_Set.P20 = as.SingleCellExperiment(Gen_Set.P20.seurat)

p10 <- plotExpression(Gen_Set.P10, 
               features = genes , #define your genes
               x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 1, 
               colour_by = "AnnotatedClusters"
               ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")

p15 <- plotExpression(Gen_Set.P15, 
                      features = genes , #define your genes
                      x = "AnnotatedClusters", #group by cluster identity
               point_alpha = 1, 
               colour_by = "AnnotatedClusters"
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "")


p20 <- plotExpression(Gen_Set.P20, 
                      features = genes , #define your genes
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