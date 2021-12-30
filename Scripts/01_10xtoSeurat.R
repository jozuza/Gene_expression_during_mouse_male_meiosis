# R script created on November, 5th, 2021
# This script will process the 10x data and save a seurat object

#load libraries
library(Matrix)
library(Seurat)
library(tidyverse)
library(R.utils)

#setwd("~/Documents/Joost_Lab_2021/ROSA26_meiosis_expression/E-MTAB-6946_10x")
data_dir <- "Raw/prep_10x/"
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz


matrix_dir = data_dir
barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")

print("reading mat")
mat <- readMM(file = matrix.path)
print("reading features")
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
print("reading barcodes")                           
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$V1
print("converting")
seurat_object = CreateSeuratObject(counts = mat)
seurat_object@meta.data$Barcode = rownames(seurat_object[[]])

# add metadata
cell_metadata <- read.csv("~/Documents/Joost_Lab_2021/ROSA26_meiosis_expression/E-MTAB-6946_10x/Raw/cell_metadata.txt", sep="")
table(cell_metadata$Sample)
ToFilter = c("P10", "P15", "P20", "P25", "P30", "P35", "P5")
cell_metadata_red <- cell_metadata %>% filter(Sample %in% ToFilter)
seurat_object <- subset(x = seurat_object, subset = Barcode %in% cell_metadata_red$Barcode)
cell_metadata_red <- cell_metadata_red[match(seurat_object@meta.data$Barcode, cell_metadata_red$Barcode),]
identical(cell_metadata_red$Barcode, seurat_object@meta.data$Barcode)

rownames(cell_metadata_red) <- cell_metadata_red$Barcode

seurat_object <- AddMetaData(object = seurat_object, metadata = cell_metadata_red)
identical(rownames(seurat_object[[]]),seurat_object@meta.data$Barcode)
saveRDS(seurat_object, file = paste0(data_dir, "SeuratObj_red"))
