# Analysis of protein data
# Lea Zillich 
# created: 16.11.2023

setwd("/path/to/")

library(readxl)
library(readr)
library(biomaRt)
library(org.Hs.eg.db)
library(ggplot2)
library(data.table)
library(dplyr)
library(clusterProfiler)
library(enrichplot)
library(rstatix)
library(colorRamp2)
library(ggpubr)
library(ComplexHeatmap)
library(pheatmap)

col_fun_rna = colorRamp2(c(0,0.25, 0.5,0.75, 1), c("#1B4D68","#28739C","#E4E4E4","#F7AB64","#C7361B"))
col_fun_prot = colorRamp2(c(0,0.25,0.5,0.75,1), c("#003E65","#3A98CC","#E4E4E4","#FCD8C1","#F7AB64"))

#load annotations
UniProt_to_Ensembl <- read_table2("1_Raw_data/protein/UniProt_to_Ensembl.tsv")
UniProt_to_GeneID <- read_table2("1_Raw_data/protein/UniProt_to_GeneID.tsv")

UniProt_to_GeneID <- UniProt_to_GeneID %>% 
   group_by(From) %>%
   sample_n(1)

#merge
anno <- merge(UniProt_to_GeneID, UniProt_to_Ensembl, by = "From")
colnames(anno) <- c("uniprot","GeneID","Ensembl")

#delete missings
control <- c("Abundances (Normalized): F1: Control","Abundances (Normalized): F2: Control",
                  "Abundances (Normalized): F3: Control","Abundances (Normalized): F4: Control")
mild <- c("Abundances (Normalized): F5: Sample","Abundances (Normalized): F6: Sample",
             "Abundances (Normalized): F7: Sample","Abundances (Normalized): F8: Sample")
moderate <- c("Abundances (Normalized): F9: Sample","Abundances (Normalized): F10: Sample",
             "Abundances (Normalized): F11: Sample","Abundances (Normalized): F12: Sample",
             "Abundances (Normalized): F13: Sample","Abundances (Normalized): F14: Sample")
severe <- c("Abundances (Normalized): F15: Sample","Abundances (Normalized): F16: Sample",
             "Abundances (Normalized): F17: Sample","Abundances (Normalized): F18: Sample",
             "Abundances (Normalized): F19: Sample","Abundances (Normalized): F20: Sample")


count_na <- function(x) sum(is.na(x))   

prot <- prot %>%
  mutate(control_na = apply(.[control], 1, count_na),
         mild_na = apply(.[mild], 1, count_na),
         moderate_na = apply(.[moderate], 1, count_na),
         severe_na = apply(.[severe], 1, count_na))

prot <- prot[prot$control_na<1&prot$mild_na<1&prot$moderate_na<1&prot$severe_na<1,]

prot$log_f1_ctrl <- log(prot$`Abundances (Normalized): F1: Control`)
prot$log_f2_ctrl <- log(prot$`Abundances (Normalized): F2: Control`)
prot$log_f3_ctrl <- log(prot$`Abundances (Normalized): F3: Control`)
prot$log_f4_ctrl <- log(prot$`Abundances (Normalized): F4: Control`)
prot$log_f5_mild <- log(prot$`Abundances (Normalized): F5: Sample`)
prot$log_f6_mild <- log(prot$`Abundances (Normalized): F6: Sample`)
prot$log_f7_mild <- log(prot$`Abundances (Normalized): F7: Sample`)
prot$log_f8_mild <- log(prot$`Abundances (Normalized): F8: Sample`)
prot$log_f9_mod <- log(prot$`Abundances (Normalized): F9: Sample`)
prot$log_f10_mod <- log(prot$`Abundances (Normalized): F10: Sample`)
prot$log_f11_mod <- log(prot$`Abundances (Normalized): F11: Sample`)
prot$log_f12_mod <- log(prot$`Abundances (Normalized): F12: Sample`)
prot$log_f13_mod <- log(prot$`Abundances (Normalized): F13: Sample`)
prot$log_f14_mod <- log(prot$`Abundances (Normalized): F14: Sample`)
prot$log_f15_sev <- log(prot$`Abundances (Normalized): F15: Sample`)
prot$log_f16_sev <- log(prot$`Abundances (Normalized): F16: Sample`)
prot$log_f17_sev <- log(prot$`Abundances (Normalized): F17: Sample`)
prot$log_f18_sev <- log(prot$`Abundances (Normalized): F18: Sample`)
prot$log_f19_sev <- log(prot$`Abundances (Normalized): F19: Sample`)
prot$log_f20_sev <- log(prot$`Abundances (Normalized): F20: Sample`)



count <- data.frame(t(prot[,c(41:60)]), keep.rownames=T)
colnames(count) <- prot$Accession
count <- count[,1:1648]
count$condition <- as.factor(c("control","control","control","control","mild","mild","mild","mild","moderate","moderate",
                     "moderate","moderate","moderate","moderate","severe","severe","severe","severe","severe","severe"))

save(count, file="WGCNA/log_counts.Rdata")

nprot <- length(prot$Accession)

results <- data.frame(uniprot = prot$Accession,
                      r_mild = rep(0, times = nprot),
                      p_mild = rep(0, times = nprot),
                      r_moderate = rep(0, times = nprot),
                      p_moderate = rep(0, times = nprot),
                      r_severe = rep(0, times = nprot),
                      p_severe = rep(0, times = nprot),
                      mean_control = rep(0, times = nprot),
                      mean_mild = rep(0, times = nprot),
                      mean_moderate = rep(0, times = nprot),
                      mean_severe = rep(0, times = nprot))

for(i in 1:nprot) {

   test_frame <- data.frame(count[,c(i,1649)])
   colnames(test_frame) <- c("protein","condition")
   wiltest <- wilcox_test(test_frame, protein~condition)
   wileff <- wilcox_effsize(test_frame, protein~condition)
   
   #Effektstärke und p-Wert
   results[i,2] <- wileff[1,4]
   results[i,3] <- wiltest[1,7]
   results[i,4] <- wileff[2,4]
   results[i,5] <- wiltest[2,7]   
   results[i,6] <- wileff[3,4]
   results[i,7] <- wiltest[3,7]
   
   #Mittelwerte
   results[i,8] <- mean(test_frame[test_frame$condition=="control","protein"])
   results[i,9] <- mean(test_frame[test_frame$condition=="mild","protein"])
   results[i,10] <- mean(test_frame[test_frame$condition=="moderate","protein"])
   results[i,11] <- mean(test_frame[test_frame$condition=="severe","protein"])
   
   #Richtung des Effekts
   results[i,2] <- ifelse(results[i,9]-results[i,8]<0,results[i,2]*(-1),results[i,2])
   results[i,4] <- ifelse(results[i,10]-results[i,8]<0,results[i,4]*(-1),results[i,4])
   results[i,6] <- ifelse(results[i,11]-results[i,8]<0,results[i,6]*(-1),results[i,6])
   
}

results$p_mild_adj <- p.adjust(results$p_mild, method="fdr")
results$p_moderate_adj <- p.adjust(results$p_moderate, method="fdr")
results$p_severe_adj <- p.adjust(results$p_severe, method="fdr")

#annotation of results
results_anno <- merge(results, UniProt_to_GeneID, by.x="uniprot", by.y="From")
colnames(results_anno)[colnames(results_anno)=="To"] <-"GeneID"
results_anno <- merge(results_anno, prot[,c(1,2)], by.x="uniprot",by.y="Accession", all.x=T)

write.table(results_anno,"/path/to/protein_results/results_wilcox_with_annotation.txt",sep=";",quote=F,row.names = F)


###### GO ANALYSIS ######

results_anno <- read_delim("/path/to/protein_results/results_wilcox_with_annotation.txt", delim = ";", escape_double = FALSE, trim_ws = TRUE)
DE_mild <- read_csv("/path/to/DEresults/progenitor_Mild.csv")
DE_moderate <- read_csv("/path/to/DEresults/progenitor_Moderate.csv")
DE_severe <- read_csv("/path/to/DEresults/progenitor_Severe.csv")

DE_mild <- DE_mild[DE_mild$p_val<0.05,]
DE_moderate <- DE_moderate[DE_moderate$p_val<0.05,]
DE_severe <- DE_severe[DE_severe$p_val<0.05,]

#mit EntrezgeneID annotieren
library(org.Hs.eg.db)
hs <- org.Hs.eg.db

ann_mild <- AnnotationDbi::select(hs, 
       keys = DE_mild$...1,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")

ann_mod <- AnnotationDbi::select(hs, 
                                 keys = DE_moderate$...1,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")

ann_severe <- AnnotationDbi::select(hs, 
                                    keys = DE_severe$...1,
                                    columns = c("ENTREZID", "SYMBOL"),
                                    keytype = "SYMBOL")
##GO COMPARE PROTEIN
# GO plots using compareCluster
GO_prot_list <- list(mild_protein=results_anno$GeneID[results_anno$p_mild<.05],moderate_protein=results_anno$GeneID[results_anno$p_moderate<.05],severe_protein=results_anno$GeneID[results_anno$p_severe<.05])

# BP 
prot_GO_BP <- compareCluster(geneClusters = GO_prot_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="BP")
ggsave("/path/to/protein_results/GO_prot_BP.pdf",p1,height=8,width=10)

# MF
prot_GO_MF <- compareCluster(geneClusters = GO_prot_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=1,qvalueCutoff=1,ont="MF")
write.table(prot_GO_MF@compareClusterResult, file="/path/to/protein_results/GO_MF.csv", sep = ",")

# CC
prot_GO_CC <- compareCluster(geneClusters = GO_prot_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="CC")
ggsave("/path/to/protein_results/GO_prots_CC.pdf",p3,height=8,width=10)

write.table(prot_GO_CC@compareClusterResult, file="/path/to/protein_results/GO_CC.csv", sep = ",")

#"#F184A7","#E8326D", "#F7AB64","#BE154C","#006960","#D1BCDC","#428E77","#74BA59","#BFE1D7","#21473C"
#"#C7361B","#F7AB64", "#313C48","#880000","#006960","#FCD8C1","#428E77","#778CA1","#BFE1D7","#21473C"

##GO COMPARE ALL
# GO plots using compareCluster
GO_gene_list <- list(mild_protein=results_anno$GeneID[results_anno$p_mild<.05],moderate_protein=results_anno$GeneID[results_anno$p_moderate<.05],severe_protein=results_anno$GeneID[results_anno$p_severe<.05],
                     mild_RNA=ann_mild$ENTREZID,moderate_RNA=ann_mod$ENTREZID,severe_RNA=ann_severe$ENTREZID)

# BP 
module_GO_BP <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="BP")
module_GO_BP_all <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=1,qvalueCutoff=1,ont="BP")
GO_BP <- module_GO_BP_all@compareClusterResult
module_GO_BP <-pairwise_termsim(module_GO_BP)

write.table(module_GO_BP@compareClusterResult, file="/path/to/protein_results/GO_compare/GO_BP_prot_RNA.csv", row.names = F,sep = ",")

p4 <- emapplot(module_GO_BP,cex.params=list(category_node=1,line=0.1,category_label=1),layout.params=list("nicely"),showCategory=12,pie.params=list(pie="equal",legend_n=3))+scale_fill_manual(values=c("#85CEE4","#3A98CC","#003E65","#F7AB64","#C7361B","#880000"))
ggsave("/path/to/protein_results/GO_compare/GO_BP_prot_RNA.pdf",p4,height=7.5,width=10)

p4b <- emapplot(module_GO_BP,cex.params=list(category_node=1,line=0.1,category_label=1),layout.params=list("nicely"),showCategory=5,pie.params=list(pie="equal",legend_n=3))+scale_fill_manual(values=c("#85CEE4","#3A98CC","#003E65","#F7AB64","#C7361B","#880000"))
ggsave("/path/to/protein_results/GO_compare/GO_BP_prot_RNA_reduced.pdf",p4b,height=5,width=6.5)

# MF
module_GO_MF <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=1,qvalueCutoff=1,ont="MF")
module_GO_MF_all <- compareCluster(geneClusters = GO_gene_list,fun = "enrichGO",OrgDb = org.Hs.eg.db,pvalueCutoff=0.25,ont="MF")
GO_MF <- module_GO_MF_all@compareClusterResult
module_GO_MF <-pairwise_termsim(module_GO_MF)

p5 <- emapplot(module_GO_MF,cex.params=list(category_node=1,line=0.1,category_label=1),layout.params=list("nicely"),showCategory=12,pie.params=list(pie="equal",legend_n=3))+scale_fill_manual(values=c("#85CEE4","#3A98CC","#003E65","#F7AB64","#C7361B","#880000"))
ggsave("/path/to/protein_results/GO_compare/GO_MF_prot_RNA.pdf",p5,height=8.5,width=10)

p5b <- emapplot(module_GO_MF,cex.params=list(category_node=1,line=0.1,category_label=1),layout.params=list("nicely"),showCategory=5,pie.params=list(pie="equal",legend_n=3))+scale_fill_manual(values=c("#85CEE4","#3A98CC","#003E65","#F7AB64","#C7361B","#880000"))
ggsave("/path/to/protein_results/GO_compare/GO_MF_prot_RNA_reduced.pdf",p5b,height=5,width=7)

write.table(module_GO_MF@compareClusterResult, file="/path/to/protein_results/GO_compare/GO_MF_prot_RNA.csv", row.names = F,sep = ",")



   #cadherin binding combined plot
GO_MF <- module_GO_MF@compareClusterResult
plot_cad <- GO_MF[GO_MF$ID=="GO:0045296",]
plot_cad$logP <- -log10(plot_cad$pvalue)

plot_cad <- plot_cad[,c("Cluster","logP")]
plot_cad$condition <- gsub("_.*","",plot_cad$Cluster)
plot_cad$omic <- gsub(".*_","",plot_cad$Cluster)
plot_cad$Cluster <- factor(plot_cad$Cluster,levels=c("mild_RNA","mild_protein","moderate_RNA","moderate_protein","severe_RNA","severe_protein"))


p7 <- ggplot(plot_cad,aes(x=condition,y=logP, fill=Cluster)) + ylim(0, 25) + guides(fill = guide_legend(nrow = 1)) +
   geom_bar(position="dodge",stat="identity") + geom_hline(yintercept = 1.31) + geom_hline(yintercept = 2.079181, linetype = "dashed") + 
   scale_fill_manual(values = c("#F7AB64","#85CEE4","#C7361B","#3A98CC","#880000","#003E65"),labels=c("mild RNA","mild protein","moderate RNA","moderate protein","severe RNA","severe protein")) +
   theme_classic() + theme(text=element_text(size=18), legend.text = element_text(size=12), legend.position = "bottom") + labs(title="cadherin binding",  x="condition", y = "-log10(p)") + scale_x_discrete(labels= c("mild","moderate","severe"))

 #actin binding

plot_act <- GO_MF[GO_MF$ID=="GO:0003779",]
plot_act$logP <- -log10(plot_act$pvalue)

plot_act <- plot_act[,c("Cluster","logP")]
plot_act$condition <- gsub("_.*","",plot_act$Cluster)
plot_act$omic <- gsub(".*_","",plot_act$Cluster)
plot_act$Cluster <- factor(plot_act$Cluster,levels=c("mild_RNA","mild_protein","moderate_RNA","moderate_protein","severe_RNA","severe_protein"))


p8 <- ggplot(plot_act,aes(x=condition,y=logP, fill=Cluster))+ ylim(0, 10) +
   geom_bar(position="dodge",stat="identity") + geom_hline(yintercept = 1.31) + geom_hline(yintercept = 2.079181, linetype = "dashed") + 
   scale_fill_manual(values = c("#F7AB64","#85CEE4","#C7361B","#3A98CC","#880000","#003E65"),labels=c("mild RNA","mild protein","moderate RNA","moderate protein","severe RNA","severe protein")) +
   theme_classic() + theme(text=element_text(size=14), legend.position = "none") + labs(title="actin binding",  x="condition", y = "-log10(p)") + scale_x_discrete(labels= c("mild","moderate","severe"))
  
    #actin filamentbinding
plot_act_fil <- GO_MF[GO_MF$ID=="GO:0051015",]
plot_act_fil$logP <- -log10(plot_act_fil$pvalue)

plot_act_fil <- plot_act_fil[,c("Cluster","logP")]
plot_act_fil$condition <- gsub("_.*","",plot_act_fil$Cluster)
plot_act_fil$omic <- gsub(".*_","",plot_act_fil$Cluster)
plot_act_fil$Cluster <- factor(plot_act_fil$Cluster,levels=c("mild_RNA","mild_protein","moderate_RNA","moderate_protein","severe_RNA","severe_protein"))


p9 <- ggplot(plot_act_fil,aes(x=condition,y=logP, fill=Cluster))+ ylim(0, 10) +
   geom_bar(position="dodge",stat="identity") + geom_hline(yintercept = 1.31) + geom_hline(yintercept = 2.079181, linetype = "dashed") + 
   scale_fill_manual(values = c("#F7AB64","#85CEE4","#C7361B","#3A98CC","#880000","#003E65"),labels=c("mild RNA","mild protein","moderate RNA","moderate protein","severe RNA","severe protein")) +
   theme_classic() + theme(text=element_text(size=14), legend.position = "none") + labs(title="actin filament binding",  x="condition", y = "-log10(p)") + scale_x_discrete(labels= c("mild","moderate","severe"))

   #tubulin binding
plot_tub <- GO_MF[GO_MF$ID=="GO:0015631",]
plot_tub$logP <- -log10(plot_tub$pvalue)

plot_tub <- plot_tub[,c("Cluster","logP")]
plot_tub$condition <- gsub("_.*","",plot_tub$Cluster)
plot_tub$omic <- gsub(".*_","",plot_tub$Cluster)
plot_tub$Cluster <- factor(plot_tub$Cluster,levels=c("mild_RNA","mild_protein","moderate_RNA","moderate_protein","severe_RNA","severe_protein"))


p10 <- ggplot(plot_tub,aes(x=condition,y=logP, fill=Cluster))+ ylim(0, 10) +
   geom_bar(position="dodge",stat="identity") + geom_hline(yintercept = 1.31) + geom_hline(yintercept = 2.079181, linetype = "dashed") + 
   scale_fill_manual(values = c("#F7AB64","#85CEE4","#C7361B","#3A98CC","#880000","#003E65"),labels=c("mild RNA","mild protein","moderate RNA","moderate protein","severe RNA","severe protein")) +
   theme_classic() + theme(text=element_text(size=14), legend.position = "none") + labs(title="tubulin binding",  x="condition", y = "-log10(p)") + scale_x_discrete(labels= c("mild","moderate","severe"))

   #cytoskeleton
plot_cyto <- GO_MF[GO_MF$ID=="GO:0005200",]
plot_cyto <- rbind(plot_cyto[1,],plot_cyto)
plot_cyto[1,1] <- "mild_protein"
plot_cyto[1,6] <- 1
plot_cyto$logP <- -log10(plot_cyto$pvalue)

plot_cyto <- plot_cyto[,c("Cluster","logP")]
plot_cyto$condition <- gsub("_.*","",plot_cyto$Cluster)
plot_cyto$omic <- gsub(".*_","",plot_cyto$Cluster)
plot_cyto$Cluster <- factor(plot_cyto$Cluster,levels=c("mild_RNA","mild_protein","moderate_RNA","moderate_protein","severe_RNA","severe_protein"))


p11 <- ggplot(plot_cyto,aes(x=condition,y=logP, fill=Cluster))+ ylim(0, 10) +
   geom_bar(position="dodge",stat="identity") + geom_hline(yintercept = 1.31) + geom_hline(yintercept = 2.079181, linetype = "dashed") + 
   scale_fill_manual(values = c("#F7AB64","#85CEE4","#C7361B","#3A98CC","#880000","#003E65"),labels=c("mild RNA","mild protein","moderate RNA","moderate protein","severe RNA","severe protein")) +
   theme_classic() + theme(text=element_text(size=14), legend.position = "none") + labs(title="cytoskeleton constituent",  x="condition", y = "-log10(p)") + scale_x_discrete(labels= c("mild","moderate","severe"))



ggsave("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/figures/Fig3A.pdf",
       ggarrange(p8,p9,p10,p11,nrow=1,common.legend = T, legend="bottom"),
       width = 6, height = 6)

ggsave("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/figures/GO_cadherin.pdf",
       p7,width = 3.5, height = 4)



#Heatmap
#prepare count data for heatmap
heatmap <- read_excel("protein/prots_genes_for_heatmap.xlsx")

library(tidyr)
#get p values
all_mild <- read_csv("DEresults/progenitor_all_marker_Mild.csv")
all_mod <- read_csv("DEresults/progenitor_all_marker_Moderate.csv")
all_sev <- read_csv("DEresults/progenitor_all_marker_Severe.csv")

heatmap <- heatmap[heatmap$Uniprot %in% results_anno$uniprot,]
heatmap <- heatmap[heatmap$name %in% all_mild$...1,]

results_reduced <- merge(results_anno,heatmap,by.x="uniprot",by.y="Uniprot")
results_reduced  <- results_reduced[order(results_reduced$name),]

all_mild <- all_mild[all_mild$...1 %in%heatmap$name,]
all_mild <- all_mild[order(all_mild$...1),]
all_mod <- all_mod[all_mod$...1 %in%heatmap$name,]
all_mod <- all_mod[order(all_mod$...1),]
all_sev <- all_sev[all_sev$...1 %in%heatmap$name,]
all_sev <- all_sev[order(all_sev$...1),]


names_heat <- UniProt_to_GeneID

stars_gen <- function(x){
  stars <- c("***", "**", "*","")
  vec <- c(0, 0.001, 0.01, 0.05,1.01)
  i <- findInterval(x, vec)
  stars[i]
}

z_mat <- data.frame(z_mild=round(results_reduced$r_mild,2),z_moderate=results_reduced$r_moderate,z_severe=results_reduced$r_severe)
rownames(z_mat) <- results_reduced$name
colnames(z_mat) <- c("mild", "moderate","severe")
z_mat <- as.matrix(z_mat)

p_mat <- data.frame(p_mild=results_reduced$p_mild,p_moderate=results_reduced$p_moderate,p_severe=results_reduced$p_severe)
rownames(p_mat) <- results_reduced$name
colnames(p_mat) <- c("mild", "moderate","severe")
p_mat[,1] <- stars_gen(p_mat[,1])
p_mat[,2] <- stars_gen(p_mat[,2])
p_mat[,3] <- stars_gen(p_mat[,3])
p_mat <- as.matrix(p_mat)


pdf("/path/to/figures/Heatmap_LIS1_interaction.pdf",height = 16,width = 8)
pheatmap(z_mat,color = col_fun_rna(seq(0,1,length=100)),cluster_rows = F,cluster_cols = F,fontsize = 16, angle_col = "45",cellwidth = 20,cellheight = 20,display_numbers = p_mat)
dev.off()

#RNA
z_rna <- data.frame(z_mild=round(all_mild$avg_log2FC,2),z_moderate=all_mod$avg_log2FC,z_severe=all_sev$avg_log2FC)
rownames(z_rna) <- all_mild$...1
colnames(z_rna) <- c("mild", "moderate","severe")
z_rna <- as.matrix(z_rna)

p_rna <- data.frame(p_mild=all_mild$p_val,p_moderate=all_mod$p_val,p_severe=all_sev$p_val)
rownames(p_rna) <- all_mild$...1
colnames(p_rna) <- c("mild", "moderate","severe")
p_rna[,1] <- stars_gen(p_rna[,1])
p_rna[,2] <- stars_gen(p_rna[,2])
p_rna[,3] <- stars_gen(p_rna[,3])
p_rna <- as.matrix(p_rna)

pdf("/path/to/figures/Heatmap_LIS1_interaction_RNA.pdf",height = 16,width = 8)
pheatmap(z_rna,color = col_fun_rna(seq(0,1,length=100)),cluster_rows = F,cluster_cols = F,fontsize = 16, angle_col = "45",cellwidth = 20,cellheight = 20,display_numbers = p_rna)
dev.off()


##HEATMAP protein
# Import DPE results - log2FC 0.5 cutoff and padj 0.05, only celltypes with more than 100 cells from each condition

results_heatmap <- results_anno[abs(results_anno$r_mild)>0.75|abs(results_anno$r_moderate)>0.75|abs(results_anno$r_severe)>0.75,]

ann_results <- AnnotationDbi::select(hs, 
                                     keys = as.character(results_heatmap$GeneID),
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENTREZID")

results_heatmap <- merge(results_heatmap, ann_results, by.x="GeneID", by.y="ENTREZID")

mild_top10 <- results_heatmap$SYMBOL[order(abs(results_heatmap$r_mild),decreasing = T)][c(1:10)]
moderate_top10 <- results_heatmap$SYMBOL[order(abs(results_heatmap$r_moderate),decreasing = T)][c(1:10)]
severe_top10 <- results_heatmap$SYMBOL[order(abs(results_heatmap$r_severe),decreasing = T)][c(1:10)]


# Create log2FC heatmap
hm_DEG <- as.matrix(results_heatmap[,c("r_mild","r_moderate","r_severe")])
rownames(hm_DEG) <- results_heatmap$SYMBOL
colnames(hm_DEG) <- c("mild","moderate","severe")

#reorder the rows to match with annotation labels

# subset genes for labeling in the heatmap
lab <- rownames(hm_DEG)
lab_genes <- c(mild_top10,moderate_top10,severe_top10)

#remove mitochondrial genes, pseudogenes and XIST

lab_genes <- lab_genes[-grep("RPS|RPL|XIST|MT-|AC|AL",lab_genes)]
lab[!(lab %in% lab_genes)] <- ""
lab_genes <- unique(lab_genes)

# Plot heatmap
pdf("/path/to/figures/SupFig3E_protein_hm.pdf",width=3.5,height=3)
Heatmap(hm_DEG,name="r",col=col_fun_prot(seq(0,1,length=20)),show_row_dend=F,cluster_columns = FALSE,column_names_gp  = gpar(fontsize=8),column_names_rot = 45) + rowAnnotation(link = anno_mark(at = which(rownames(hm_DEG)==lab), 
                                                                                                                                                  labels = lab[lab !=""], padding = unit(1, "mm"),link_gp = gpar(lwd = 0.1) ,labels_gp = gpar(fontsize = 6)))
dev.off()
