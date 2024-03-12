## WGCNA - associate networks with traits
# last change: LZ 2023-11

library(WGCNA)
library(data.table)
library(dplyr)
library(readr)
library(missMethyl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

DF <- data.frame

setwd("/path/to/WGCNA")

options(stringsAsFactors = FALSE)


  # Load the expression and trait data saved in the first part
lnames = load(file = "WGCNA_input.Rdata")  
  #The variable lnames contains the names of loaded variables.
  lnames
  
  setwd("/path/to/WGCNA/WGCNA_prot")
  # Load network data saved in the second part.
  lnames = load(file = "networkConstruction_prot.RData");
  lnames
  
  #prepare traits
  datTraits$mild <- NA
  datTraits$mild[datTraits$condition=="mild"] <- 1
  datTraits$mild[datTraits$condition=="control"] <- 0
  
  datTraits$moderate <- NA
  datTraits$moderate[datTraits$condition=="moderate"] <- 1
  datTraits$moderate[datTraits$condition=="control"] <- 0
  
  datTraits$severe <- NA
  datTraits$severe[datTraits$condition=="severe"] <- 1
  datTraits$severe[datTraits$condition=="control"] <- 0
  
  # Define numbers of genes and samples
  nGenes = ncol(datExpr0);
  nSamples = nrow(datExpr0);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits[,2:4], use = "p")
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
  write.table(MEs, file = "MEs_Expr.txt", sep = ";", quote = F, row.names = T)
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  
  setwd("/path/to/WGCNA") 
  
  png("heatmap_modules_prot.png", width = 600, height = 500)
  par(mar = c(6, 11.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits[,2:4]),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 1.8,
                 cex.lab.x = 2,
                 cex.lab.y = 2,
                 zlim = c(-1,1))
  dev.off()
  

  
  # Define variable weight containing the weight column of datTrait
  
  # names (colors) of the modules
  modNames = substring(names(MEs), 3)
  geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  
  names(geneModuleMembership) = paste("MM", modNames, sep="");
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits[2:4], use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
  names(geneTraitSignificance) = paste("GS.", names(datTraits[2:4]), sep="");
  names(GSPvalue) = paste("p.GS.", names(datTraits[2:4]), sep="")
  
  mod <- DF(moduleTraitCor)
  module = "blue"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  sizeGrWindow(7, 7);
  par(mfrow = c(1,1));
  
  
  png("module_membership_gene_significance_Expr.png")
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", module, "module"),
                     ylab = "Gene significance for AUD status",
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 3, cex.lab = 3, cex.axis = 3, col = module)
  dev.off()
  
  write.table(mod, "mod_trait_cor_prot.txt", sep = ";", quote = F)
  
  #read annotation
  anno <- read_table2("/path/to/protein/UniProt_to_GeneID.tsv")
  

  dim(anno)
  names(anno)
  probes = colnames(datExpr0)
  probes2annot = match(probes, anno$From)
  # The following is the number or probes without annotation:
  sum(is.na(probes2annot))
  
  
  # Create the starting data frame
  geneInfo0 = data.frame(genx = probes,
                         geneSymbol = anno$To[probes2annot],
                         moduleColor = moduleColors,
                         geneTraitSignificance,
                         GSPvalue)
  # Order modules by their significance for weight
  modOrder = order(-abs(cor(MEs, datTraits[,4], use = "p")));
  # Add module membership information in the chosen order
  for (mod in 1:ncol(geneModuleMembership))
  { oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
  
  # Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
  geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.severe));
  geneInfo = geneInfo0[geneOrder, ]
  
  write.table(geneInfo, file = "geneInfo.txt", sep = ";", quote = F, row.names = F)
  write.table(table(geneInfo$moduleColor), file = "geneInfo_noMods.txt", sep = ";", quote = F, row.names = F)
  
 
    # GO plots using compareCluster
  GO_gene_list <- list(blue=geneInfo0$geneSymbol[geneInfo0$moduleColor == "blue"],turquoise=geneInfo0$geneSymbol[geneInfo0$moduleColor == "turquoise"])
  
  # BP 
  
  module_GO_BP <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="BP")
  module_GO_BP <-pairwise_termsim(module_GO_BP)
  
  write.table(module_GO_BP@compareClusterResult, file="/path/to/WGCNA/GO_BP.csv", sep = ",")
  
  p1 <- emapplot(module_GO_BP,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=12,pie="equal")+scale_fill_manual(values=c("blue", "turquoise"))
  ggsave("GO_modules_BP.pdf",p1,height=6,width=7)
  
  # MF
  module_GO_MF <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="MF")
  module_GO_MF <-pairwise_termsim(module_GO_MF)
  
  p2 <- emapplot(module_GO_MF,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=12,pie="equal")+scale_fill_manual(values=c("blue", "turquoise"))
  ggsave("GO_modules_MF.pdf",p2,height=6,width=7)
  
  write.table(module_GO_MF@compareClusterResult, file="/path/to/WGCNA/GO_MF.csv", sep = ",")
  
  
  # CC
  module_GO_CC <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="CC")
  module_GO_CC <-pairwise_termsim(module_GO_CC)
  
  p3 <- emapplot(module_GO_CC,legend_n=3,cex_line=0.1,cex_label_category=0.75,layout="nicely",cex_category=0.8,showCategory=12,pie="equal")+scale_fill_manual(values=c("blue", "turquoise"))
  ggsave("GO_modules_CC.pdf",p3,height=6,width=7)
  
  write.table(module_GO_CC@compareClusterResult, file="/path/to/WGCNA/GO_CC.csv", sep = ",")
  
  
  