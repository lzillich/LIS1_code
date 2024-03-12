# Explore the results of drug signature analysis 
library(ggplot2)
library(forcats)
library(data.table)

# Mild

mi_query_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Mild/arfs/TAG/query_result.gct",skip=2)
mi_query_result <-mi_query_result[-1,]
#-log10(0.05) -> 1.30103

mi_query_result_sig <-mi_query_result[abs(as.numeric(mi_query_result$fdr_q_nlog10)) > 1.30103,]
mi_neg_corr_perturbagens <- mi_query_result_sig[order(as.numeric(mi_query_result_sig$norm_cs),decreasing = F),]

# GSEA results

mi_gsea_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Mild/gsea/TAG/arfs/NORM_CS/gsea_result.gct",skip=2)
mi_gsea_result <- mi_gsea_result[-1,]

mi_gsea_result_sig_MOA <- mi_gsea_result[mi_gsea_result$set_type == "MOA_CLASS",]
mi_gsea_result_sig_MOA <- mi_gsea_result_sig_MOA[order(mi_gsea_result_sig_MOA$fdr_q_nlog10,decreasing=T),]
mi_gsea_result_sig_PATH <- mi_gsea_result[mi_gsea_result$set_type == "PATHWAY_SET",]
mi_gsea_result_sig_PATH <- mi_gsea_result_sig_PATH[order(mi_gsea_result_sig_PATH$fdr_q_nlog10,decreasing=T),]
mi_gsea_result_sig_PCL <- mi_gsea_result[mi_gsea_result$set_type == "PCL",]
mi_gsea_result_sig_PCL <- mi_gsea_result_sig_PCL[order(mi_gsea_result_sig_PCL$fdr_q_nlog10,decreasing=T),]

# MOA

mi_neg_corr_gsea_MOA <- mi_gsea_result_sig_MOA[order(as.numeric(mi_gsea_result_sig_MOA$norm_cs),decreasing = F),]
mi_neg_corr_gsea_MOA$cell_iname[mi_neg_corr_gsea_MOA$cell_iname=="-666"]<-""
mi_neg_corr_gsea_MOA$name <- paste0(mi_neg_corr_gsea_MOA$src_set_id,"_",mi_neg_corr_gsea_MOA$cell_iname)
mi_neg_corr_gsea_MOA$name <- gsub("*_$","",mi_neg_corr_gsea_MOA$name)

# PCL

mi_neg_corr_gsea_PCL <- mi_gsea_result_sig_PCL[order(as.numeric(mi_gsea_result_sig_PCL$norm_cs),decreasing = F),]
mi_neg_corr_gsea_PCL$src_set_id <- gsub("^...","",neg_corr_gsea_PCL$src_set_id)
mi_neg_corr_gsea_PCL$cell_iname[mi_neg_corr_gsea_PCL$cell_iname=="-666"]<-""
mi_neg_corr_gsea_PCL$name <- paste0(mi_neg_corr_gsea_PCL$src_set_id,"_",mi_neg_corr_gsea_PCL$cell_iname)
mi_neg_corr_gsea_PCL$name <- gsub("*_$","",mi_neg_corr_gsea_PCL$name)

# Moderate
# Individual perturbagens 

m_query_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Moderate/arfs/TAG/query_result.gct",skip=2)
m_query_result <-m_query_result[-1,]
#-log10(0.05) -> 1.30103

m_query_result_sig <-m_query_result[abs(as.numeric(m_query_result$fdr_q_nlog10)) > 1.30103,]

m_neg_corr_perturbagens <- m_query_result_sig[order(as.numeric(m_query_result_sig$norm_cs),decreasing = F),]

# GSEA results

m_gsea_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Moderate/gsea/TAG/arfs/NORM_CS/gsea_result.gct",skip=2)
m_gsea_result <- m_gsea_result[-1,]

m_gsea_result_sig_MOA <- m_gsea_result[m_gsea_result$set_type == "MOA_CLASS",]
m_gsea_result_sig_MOA <- m_gsea_result_sig_MOA[order(m_gsea_result_sig_MOA$fdr_q_nlog10,decreasing=T),]
m_gsea_result_sig_PCL <- m_gsea_result[m_gsea_result$set_type == "PCL",]
m_gsea_result_sig_PCL <- m_gsea_result_sig_PCL[order(m_gsea_result_sig_PCL$fdr_q_nlog10,decreasing=T),]

# MOA
m_neg_corr_gsea_MOA <- m_gsea_result_sig_MOA[order(as.numeric(m_gsea_result_sig_MOA$norm_cs),decreasing = F),]
m_neg_corr_gsea_MOA$cell_iname[m_neg_corr_gsea_MOA$cell_iname=="-666"]<-""
m_neg_corr_gsea_MOA$name <- paste0(m_neg_corr_gsea_MOA$src_set_id,"_",m_neg_corr_gsea_MOA$cell_iname)
m_neg_corr_gsea_MOA$name <- gsub("*_$","",m_neg_corr_gsea_MOA$name)

# PCL
m_neg_corr_gsea_PCL <- m_gsea_result_sig_PCL[order(as.numeric(m_gsea_result_sig_PCL$norm_cs),decreasing = F),]
m_neg_corr_gsea_PCL$src_set_id <- gsub("^...","",m_neg_corr_gsea_PCL$src_set_id)
m_neg_corr_gsea_PCL$cell_iname[m_neg_corr_gsea_PCL$cell_iname=="-666"]<-""
m_neg_corr_gsea_PCL$name <- paste0(m_neg_corr_gsea_PCL$src_set_id,"_",m_neg_corr_gsea_PCL$cell_iname)
m_neg_corr_gsea_PCL$name <- gsub("*_$","",m_neg_corr_gsea_PCL$name)

# Severe 

# Individual perturbagens 

query_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Severe/arfs/TAG/query_result.gct",skip=2)
query_result <-query_result[-1,]
#-log10(0.05) -> 1.30103

query_result_sig <-query_result[abs(as.numeric(query_result$fdr_q_nlog10)) > 1.30103,]
neg_corr_perturbagens <- query_result_sig[order(as.numeric(query_result_sig$norm_cs),decreasing = F),]

gsea_result <- fread("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Severe/gsea/TAG/arfs/NORM_CS/gsea_result.gct",skip=2)
gsea_result <- gsea_result[-1,]
gsea_result_sig_MOA <- gsea_result[gsea_result$set_type == "MOA_CLASS",]
gsea_result_sig_MOA <- gsea_result_sig_MOA[order(gsea_result_sig_MOA$fdr_q_nlog10,decreasing=T),]
gsea_result_sig_PCL <- gsea_result[gsea_result$set_type == "PCL",]
gsea_result_sig_PCL <- gsea_result_sig_PCL[order(gsea_result_sig_PCL$fdr_q_nlog10,decreasing=T),]

# MOA
neg_corr_gsea_MOA <- gsea_result_sig_MOA[order(as.numeric(gsea_result_sig_MOA$norm_cs),decreasing = F),]
neg_corr_gsea_MOA$cell_iname[neg_corr_gsea_MOA$cell_iname=="-666"]<-""
neg_corr_gsea_MOA$name <- paste0(neg_corr_gsea_MOA$src_set_id,"_",neg_corr_gsea_MOA$cell_iname)
neg_corr_gsea_MOA$name <- gsub("*_$","",neg_corr_gsea_MOA$name)

# PCL
neg_corr_gsea_PCL <- gsea_result_sig_PCL[order(as.numeric(gsea_result_sig_PCL$norm_cs),decreasing = F),]
neg_corr_gsea_PCL$src_set_id <- gsub("^...","",neg_corr_gsea_PCL$src_set_id)
neg_corr_gsea_PCL$cell_iname[neg_corr_gsea_PCL$cell_iname=="-666"]<-""
neg_corr_gsea_PCL$name <- paste0(neg_corr_gsea_PCL$src_set_id,"_",neg_corr_gsea_PCL$cell_iname)
neg_corr_gsea_PCL$name <- gsub("*_$","",neg_corr_gsea_PCL$name)

# Generate plots

#perturbagens
mi_neg_corr_perturbagens <- mi_neg_corr_perturbagens[,c("pert_iname","fdr_q_nlog10","norm_cs")]
colnames(mi_neg_corr_perturbagens) <- c("compound","pval","cscore_norm")
mi_neg_corr_perturbagens$pval <- as.numeric(mi_neg_corr_perturbagens$pval)
mi_neg_corr_perturbagens$cscore_norm <- as.numeric(mi_neg_corr_perturbagens$cscore_norm)

mi_neg_corr_perturbagens <- mi_neg_corr_perturbagens[!duplicated(mi_neg_corr_perturbagens$compound),]
mi_neg_corr_perturbagens <- mi_neg_corr_perturbagens[order(mi_neg_corr_perturbagens$cscore_norm,decreasing = F),]
mi_neg_corr_perturbagens$condition <- "MILD"
colnames(mi_neg_corr_perturbagens)[c(2:4)]<-paste0(colnames(mi_neg_corr_perturbagens)[c(2:4)],"_MILD")


m_neg_corr_perturbagens <- m_neg_corr_perturbagens[,c("pert_iname","fdr_q_nlog10","norm_cs")]
colnames(m_neg_corr_perturbagens) <- c("compound","pval","cscore_norm")
m_neg_corr_perturbagens$pval <- as.numeric(m_neg_corr_perturbagens$pval)
m_neg_corr_perturbagens$cscore_norm <- as.numeric(m_neg_corr_perturbagens$cscore_norm)

m_neg_corr_perturbagens <- m_neg_corr_perturbagens[!duplicated(m_neg_corr_perturbagens$compound),]
m_neg_corr_perturbagens <- m_neg_corr_perturbagens[order(m_neg_corr_perturbagens$cscore_norm,decreasing = F),]
m_neg_corr_perturbagens$condition <- "MODERATE"
colnames(m_neg_corr_perturbagens)[c(2:4)]<-paste0(colnames(m_neg_corr_perturbagens)[c(2:4)],"_MODERATE")


neg_corr_perturbagens <- neg_corr_perturbagens[,c("pert_iname","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_perturbagens) <- c("compound","pval","cscore_norm")
neg_corr_perturbagens$pval <- as.numeric(neg_corr_perturbagens$pval)
neg_corr_perturbagens$cscore_norm <- as.numeric(neg_corr_perturbagens$cscore_norm)

neg_corr_perturbagens <- neg_corr_perturbagens[!duplicated(neg_corr_perturbagens$compound),]
neg_corr_perturbagens <- neg_corr_perturbagens[order(neg_corr_perturbagens$cscore_norm,decreasing = F),]
neg_corr_perturbagens$condition <- "SEVERE"
colnames(neg_corr_perturbagens)[c(2:4)]<-paste0(colnames(neg_corr_perturbagens)[c(2:4)],"_SEVERE")

pert <- merge(mi_neg_corr_perturbagens,m_neg_corr_perturbagens,by="compound",all=T)
pert<-merge(pert,neg_corr_perturbagens,by="compound",all=T)


pert_1 <- pert[!is.na(pert$cscore_norm_MODERATE) & !is.na(pert$cscore_norm_SEVERE) &!is.na(pert$cscore_norm_MILD) & pert$pval_MODERATE > -log10(0.05) & pert$pval_SEVERE > -log10(0.05)& pert$pval_MILD > -log10(0.05) & pert$cscore_norm_MILD<0 & pert$cscore_norm_MODERATE<0 & pert$cscore_norm_SEVERE<0,]
pert_1_df <- data.frame(compound=c(pert_1$compound,pert_1$compound,pert_1$compound),pval=c(pert_1$pval_MILD,pert_1$pval_MODERATE,pert_1$pval_SEVERE),cscore_norm=c(pert_1$cscore_norm_MILD,pert_1$cscore_norm_MODERATE,pert_1$cscore_norm_SEVERE),condition=c(pert_1$condition_MILD,pert_1$condition_MODERATE,pert_1$condition_SEVERE))


pert_1_df$condition <- factor(pert_1_df$condition,levels=c("SEVERE","MODERATE","MILD"))

stars_gen <- function(x){
  stars <- rev(c( "***", "**", "*",""))
  vec <- rev(c(1e6, 3, 2, 1.30103,0.0000000001))
  i <- findInterval(x, vec)
  stars[i]
}

pert_1_df$pval <- stars_gen(as.numeric(pert_1_df$pval)) 

pert_1_df$score <- 0
for(i in pert_1_df$compound){
  pert_1_df$score[pert_1_df$compound==i] <-  sum(pert_1_df$cscore_norm[pert_1_df$compound==i])
}
pert_1_df <- pert_1_df[order(pert_1_df$score,decreasing = F),]
pert_1_df$compound <- factor(pert_1_df$compound,levels=rev(unique(pert_1_df$compound)))

p1 <- ggplot(pert_1_df, aes(fill=condition, y=cscore_norm, x=compound,label=pval)) + coord_flip()+
  geom_bar(position="dodge", stat="identity")+ labs(x="", y="Normalized connectivity score",
                                                    title="",legend="sig")+geom_text(aes(label=pval,group=condition),position=position_dodge(width=0.9),hjust=1.5,vjust=0.7)+scale_fill_manual(values=c("#880000","#C7361B","#F7AB64")) +theme_minimal() + theme(axis.text = element_text(size=14),legend.text=element_text(size=12))
ggsave("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/Combined_perturbagen.pdf",p1,height=8,width=8)


# GSEA MOA

mi_neg_corr_gsea_MOA <- mi_neg_corr_gsea_MOA[,c("name","fdr_q_nlog10","norm_cs")]
colnames(mi_neg_corr_gsea_MOA) <- c("term","pval","cscore_norm")
mi_neg_corr_gsea_MOA$pval <- as.numeric(mi_neg_corr_gsea_MOA$pval)
mi_neg_corr_gsea_MOA$cscore_norm <- as.numeric(mi_neg_corr_gsea_MOA$cscore_norm)

mi_neg_corr_gsea_MOA <- mi_neg_corr_gsea_MOA[!duplicated(mi_neg_corr_gsea_MOA$term),]
mi_neg_corr_gsea_MOA <- mi_neg_corr_gsea_MOA[order(mi_neg_corr_gsea_MOA$cscore_norm,decreasing = F),]
mi_neg_corr_gsea_MOA$condition <- "MILD"
colnames(mi_neg_corr_gsea_MOA)[c(2:4)]<-paste0(colnames(mi_neg_corr_gsea_MOA)[c(2:4)],"_MILD")

m_neg_corr_gsea_MOA <- m_neg_corr_gsea_MOA[,c("name","fdr_q_nlog10","norm_cs")]
colnames(m_neg_corr_gsea_MOA) <- c("term","pval","cscore_norm")
m_neg_corr_gsea_MOA$pval <- as.numeric(m_neg_corr_gsea_MOA$pval)
m_neg_corr_gsea_MOA$cscore_norm <- as.numeric(m_neg_corr_gsea_MOA$cscore_norm)

m_neg_corr_gsea_MOA <- m_neg_corr_gsea_MOA[!duplicated(m_neg_corr_gsea_MOA$term),]
m_neg_corr_gsea_MOA <- m_neg_corr_gsea_MOA[order(m_neg_corr_gsea_MOA$cscore_norm,decreasing = F),]
m_neg_corr_gsea_MOA$condition <- "MODERATE"
colnames(m_neg_corr_gsea_MOA)[c(2:4)]<-paste0(colnames(m_neg_corr_gsea_MOA)[c(2:4)],"_MODERATE")

neg_corr_gsea_MOA <- neg_corr_gsea_MOA[,c("name","fdr_q_nlog10","norm_cs")]
colnames(neg_corr_gsea_MOA) <- c("term","pval","cscore_norm")
neg_corr_gsea_MOA$pval <- as.numeric(neg_corr_gsea_MOA$pval)
neg_corr_gsea_MOA$cscore_norm <- as.numeric(neg_corr_gsea_MOA$cscore_norm)

neg_corr_gsea_MOA <- neg_corr_gsea_MOA[!duplicated(neg_corr_gsea_MOA$term),]
neg_corr_gsea_MOA <- neg_corr_gsea_MOA[order(neg_corr_gsea_MOA$cscore_norm,decreasing = F),]
neg_corr_gsea_MOA$condition <- "SEVERE"
colnames(neg_corr_gsea_MOA)[c(2:4)]<-paste0(colnames(neg_corr_gsea_MOA)[c(2:4)],"_SEVERE")

pert2 <- merge(mi_neg_corr_gsea_MOA,m_neg_corr_gsea_MOA,by="term",all=T)
pert2 <-merge(pert2,neg_corr_gsea_MOA,by="term",all=T)


pert_2 <- pert2[ pert2$cscore_norm_MODERATE<0 & pert2$cscore_norm_SEVERE<0 ,]

pert_2 <- pert_2[order(pert_2$pval_MODERATE,pert_2$pval_SEVERE,decreasing = T),]

pert_2_df <- data.frame(term=c(pert_2$term,pert_2$term,pert_2$term),pval=c(pert_2$pval_MILD,pert_2$pval_MODERATE,pert_2$pval_SEVERE),cscore_norm=c(pert_2$cscore_norm_MILD,pert_2$cscore_norm_MODERATE,pert_2$cscore_norm_SEVERE),condition=c(pert_2$condition_MILD,pert_2$condition_MODERATE,pert_2$condition_SEVERE))

pert_2_df$pval[pert_2_df$pval ==0] <- 0.0001

pert_2_df$pval <- stars_gen(pert_2_df$pval)

pert_2_df$score <- 0
for(i in pert_2_df$term){
  pert_2_df$score[pert_2_df$term==i] <-  sum(pert_2_df$cscore_norm[pert_2_df$term==i])
}

pert_2_df <- pert_2_df[order(pert_2_df$cscore_norm,decreasing = F),]
pert_2_df$term <- factor(pert_2_df$term,levels=rev(unique(pert_2_df$term)))
pert_2_df$condition <- factor(pert_2_df$condition,levels=c("SEVERE","MODERATE","MILD"))
unique(pert_2_df$term[c(1:12)])
pert_2_df2 <- pert_2_df[pert_2_df$term %in% unique(pert_2_df$term[c(1:12)]),]

p2 <- ggplot(pert_2_df2, aes(fill=condition, y=cscore_norm, x=term,label=pval)) + coord_flip()+
  geom_bar(position="dodge", stat="identity")+ylim(-2.5,1.5) +labs(x="", y="Normalized connectivity score",
                                                    title="",legend="sig")+geom_text(aes(label=pval,group=condition),position=position_dodge(width=0.9),hjust=1.5,vjust=0.7)+scale_fill_manual(values=c("#880000","#C7361B","#F7AB64")) +theme_minimal() + theme(axis.text = element_text(size=14),legend.text=element_text(size=12))
ggsave("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/MOA_Severe_Moderate.pdf",p2,height=8,width=8)

# all within the same direction
pert_3 <- pert2[ pert2$cscore_norm_MODERATE<0 & pert2$cscore_norm_SEVERE<0 & pert2$cscore_norm_MILD<0,]

pert_3 <- pert_3[order(pert_3$pval_MODERATE,pert_3$pval_SEVERE,pert_3$pval_MILD,decreasing = T),]

pert_3_df <- data.frame(term=c(pert_3$term,pert_3$term,pert_3$term),pval=c(pert_3$pval_MILD,pert_3$pval_MODERATE,pert_3$pval_SEVERE),cscore_norm=c(pert_3$cscore_norm_MILD,pert_3$cscore_norm_MODERATE,pert_3$cscore_norm_SEVERE),condition=c(pert_3$condition_MILD,pert_3$condition_MODERATE,pert_3$condition_SEVERE))

pert_3_df$pval[pert_3_df$pval ==0] <- 0.0001

pert_3_df$pval <- stars_gen(pert_3_df$pval)

pert_3_df$score <- 0
for(i in pert_3_df$term){
  pert_3_df$score[pert_3_df$term==i] <-  sum(pert_3_df$cscore_norm[pert_3_df$term==i])
}

pert_3_df <- pert_3_df[order(pert_3_df$cscore_norm,decreasing = F),]
pert_3_df$term <- factor(pert_3_df$term,levels=rev(unique(pert_3_df$term)))
pert_3_df$condition <- factor(pert_3_df$condition,levels=c("SEVERE","MODERATE","MILD"))
unique(pert_3_df$term[c(1:10)])
pert_3_df2 <- pert_3_df[pert_3_df$term %in% unique(pert_3_df$term[c(1:10)]),]

p3 <- ggplot(pert_3_df2, aes(fill=condition, y=cscore_norm, x=term,label=pval)) + coord_flip()+
  geom_bar(position="dodge", stat="identity")+ylim(-2,0) +labs(x="", y="Normalized connectivity score",
                                                                   title="",legend="sig")+geom_text(aes(label=pval,group=condition),position=position_dodge(width=0.9),hjust=1.5,vjust=0.7)+scale_fill_manual(values=c("#880000","#C7361B","#F7AB64")) +theme_minimal() + theme(axis.text = element_text(size=14),legend.text=element_text(size=12))
ggsave("/zi-flstorage/group_genepi/shared-scRNAseq/1_Projekt_Andrea/drug_repurposing/CMap/MOA_Mild_AND_Severe_Moderate.pdf",p3,height=8,width=8)

