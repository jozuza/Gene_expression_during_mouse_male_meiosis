[E-MTAB-6946](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6946/) scRNA-Seq: 8 timepoints (P5x2, P10, P15, P25, P30, P35)

**E-MTAB-6946 - Single-cell RNA sequencing of mouse germ cells**

|File|Size|Date|What|
|:-:|:-:|:-:|:-:
|E-MTAB-6946.processed.1.zip|2.2 MB|  9 May 2019, 11:56| cell_metadata.txt
|E-MTAB-6946.processed.2.zip|265 KB|  9 May 2019, 11:56| genes.tsv
|E-MTAB-6946.processed.3.zip|676.8 MB|9 May 2019, 11:56|raw_counts.mtx
|E-MTAB-6946.processed.4.zip|496.6 MB|9 May 2019, 11:56|raw_counts_emptyDrops.mtx
|E-MTAB-6946.processed.5.zip|4.3 MB|  9 May 2019, 11:56|cell_metadata_emptyDrops.txt
|E-MTAB-6946.processed.6.zip|257 KB|  9 May 2019, 11:56|genes_emptyDrops.tsv
---


Meiosis: P10, 15, 20 
```
#scRNA-Seq (X) or bulk RNA-Seq (B)
```


|Sample|Library|
|-|-|
|P10|do17821|
|P15|do18195|
|P20|do17824|

---
## Data 
```
|── Raw
│   ├── cell_metadata.txt
│   ├── genes.tsv
│   ├── prep_10x
│   │   ├── barcodes.tsv.gz
│   │   ├── features.tsv.gz
│   │   ├── matrix.mtx.gz
│   │   ├── SeuratObj.rds
│   │   └── SeuratObj_red
│   └── raw_counts.mtx
|── Seurat_processed
    ├── Meiosis_sce.rds
    └── Meiosis_SeuratObj.rds
```
---
## scripts
01_10xtoSeurat.R
02_Filtering_and_Clustering.R
03_ROSA_Expression.R
