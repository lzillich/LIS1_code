#Environment and relevant functions for sc analysis
# Single-Cell Analysis
library(Seurat)
library(SeuratWrappers)
library(Signac)

# Data Visualization
library(ggplot2)
library(patchwork)
library(ggpubr)
library(viridis)
library(cowplot)
library(ComplexHeatmap)
library(circlize)

# Data Manipulation
library(dplyr)
library(Matrix)

# Gene Set Enrichment and Pathway Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# Motif Analysis and Genomics
#library(chromVARmotifs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(JASPAR2020)
library(TFBSTools)
library(GenomeInfoDb)

# Connectivity Map (CMap) Analysis
library(cmapR)
library(httr)
library(jsonlite)

# Color Palettes
library(RColorBrewer)

# Stuff that I am not sorting now
library(DT)
library(scales)
library(purrr)
library(DoubletFinder)

#General parameters
#mem.maxVSize(16384*5) 


#Palettes
rna_palette <- colorRampPalette(c("#cdcdcd", "#b3cde0","#6497b1","#005b96","#03396c","#011f4b"))
diverging_palette<-colorRampPalette(c("#0057B7","#4FA8F0","#F0F0F0","#F58B3C","#D65A00"))

