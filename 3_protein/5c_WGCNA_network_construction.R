## Data Preparation WGCNA 
# last change: LZ 2023-11

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(varhandle)
library(matrixStats)

DF <- data.frame
setwd("/path/to/WGCNA")

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()


        lnames = load(file = "WGCNA_input.Rdata")    
        # Choose a set of soft-thresholding powers
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        
        # Call the network topology analysis function
        sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
        
        sft$powerEstimate
        
        setwd("/path/to/WGCNA_prot")
        
        net = blockwiseModules(datExpr0, power = sft$powerEstimate,
                               TOMType = "signed", minModuleSize = 30, maxBlockSize = 36000,
                               reassignThreshold = 0, mergeCutHeight = 0.25,
                               numericLabels = TRUE, pamRespectsDendro = FALSE,
                               saveTOMs = TRUE,
                               saveTOMFileBase = "prot_TOM",
                               verbose = 3)
        table(net$colors)
        
        # Convert labels to colors for plotting
        mergedColors = labels2colors(net$colors)
        # Plot the dendrogram and the module colors underneath
        png("dendogram.png", width = 1200, height = 1800)
        plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                            "Module colors",
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05)
        dev.off()
        
        moduleLabels = net$colors
        moduleColors = labels2colors(net$colors)
        MEs = net$MEs;
        geneTree = net$dendrograms[[1]];
        save(MEs, moduleLabels, moduleColors, geneTree,
             file = "networkConstruction_prot.RData")
     