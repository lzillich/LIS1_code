## Data Preparation WGCNA 
# Written by LZ 2021-07-14
# last change: LZ 2021-07-16

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)
library(DESeq2)

DF <- data.frame

setwd("/path/to/")

options(stringsAsFactors = FALSE);

  get(load("WGCNA/log_counts.Rdata"))


WGCNA_input <- read_csv("WGCNA_input.csv")

datExpr0 = as.data.frame(WGCNA_input[, -1])
rownames(datExpr0) = WGCNA_input$sample

  #check if genes have too many missings
  gsg = goodSamplesGenes(datExpr0, verbose = 3);
  gsg$allOK
  
  if (!gsg$allOK)
  {
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
  }

  sampleTree = hclust(dist(datExpr0), method = "average");
  # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
  # The user should change the dimensions if the window is too large or too small.
  sizeGrWindow(12,9)
  #pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
  par(cex = 0.6);
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
       cex.axis = 1.5, cex.main = 2)
  
  #traitData
  phe <- read_csv("WGCNA_pheno.csv")
  allTraits = phe
  
  # Form a data frame analogous to expression data that will hold the clinical traits.
  Samples = rownames(datExpr0);
  traitRows = match(Samples, allTraits$sample);
  datTraits = allTraits[traitRows, -1];
  rownames(datTraits) = allTraits$sample[traitRows] 
  
  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr0), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  datTraits$condition <- as.factor(datTraits$condition)
  traitColors = numbers2colors(as.numeric(datTraits$condition), signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  
  png(file = "WGCNA/sampleClustering.png", width = 12, height = 9);
  
  plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits$condition),
                      main = "Sample dendrogram and trait heatmap")
  dev.off()

  save(datExpr0, datTraits, file = "WGCNA/WGCNA_input.Rdata")
  
