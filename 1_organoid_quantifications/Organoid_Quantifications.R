# morphological analyses 
# last update: 17.06.25

library(readr)
library(data.table)
library(readxl)
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("/path/to/")

loop <- read_excel("loop_parameters.xlsx", 
                                 sheet = "loop_sum")


##QUANTIFICATIONS

#VZ diameter
kruskal.test(loop$VZ_dia_mean~loop$condition)
pairwise.wilcox.test(loop$VZ_dia_mean,loop$condition, p.adjust="bonferroni")

com_con <- list(c("control", "mild"), c("control", "moderate"), c("control", "severe"), c("mild", "moderate"),c("mild", "severe"), c("moderate", "severe"))

p1 <- ggplot(loop,aes(x=condition, y=VZ_dia_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ylim(0,110) +
  xlab("") + labs(y="µm") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list( c("control", "mild"), c("control", "moderate"), c("control", "severe"), c("mild", "severe")))

ggsave(file="/path/to/figures/Fig1D_VZ_diameter_rev.pdf",p1, width = 5.5, height = 5)

#ventricle size
kruskal.test(loop$ven_size_mean~loop$condition)
pairwise.wilcox.test(loop$ven_size_mean,loop$condition, p.adjust="bonferroni")

com_ven <- list(c("control", "moderate"), c("control", "severe") )

p2 <- ggplot(loop,aes(x=condition, y=ven_size_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + labs(y=bquote("µm"^2)) + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list( c("control", "moderate"), c("control", "severe"),c("mild", "severe"), c("moderate", "severe")))


ggsave(file="/path/to/figures/SupFig1D_venlike.pdf",p2, width = 5.5, height = 5)

#apical membrane length
kruskal.test(loop$ap_mem_mean~loop$condition)
pairwise.wilcox.test(loop$ap_mem_mean,loop$condition, p.adjust="bonferroni")

p3 <- ggplot(loop,aes(x=condition, y=ap_mem_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + labs(y="µm") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354"))+ 
  stat_compare_means(comparisons = list(c("control", "mild"), c("control", "moderate"), c("control", "severe"),c("mild", "severe")))

ggsave(file="/path/to/figures/Figure1E_apical_membrane_length_rev.pdf",p3, width = 5.5, height = 5)

#total loop size
kruskal.test(loop$tot_loop_mean~loop$condition)
pairwise.wilcox.test(loop$tot_loop_mean,loop$condition, p.adjust="bonferroni")

p4 <- ggplot(loop,aes(x=condition, y=tot_loop_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + labs(y=bquote("µm"^2))  + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list(c("control", "severe"),c("mild", "severe"), c("moderate", "severe")))

ggsave(file="/path/to/figures/SuppFig1G_total_VZ_area.pdf",p4, width = 5.5, height = 5)

#basal membrane length
kruskal.test(loop$bas_mem_mean~loop$condition)
pairwise.wilcox.test(loop$bas_mem_mean,loop$condition, p.adjust="bonferroni")

p5 <- ggplot(loop,aes(x=condition, y=bas_mem_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + labs(y="µm") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list( c("control", "moderate"), c("control", "severe")))


ggsave(file="/path/to/figures/SuppFig1C_basal_membrane_length.pdf",p5, width = 5.5, height = 5)


#loop area
kruskal.test(loop$loop_area_mean~loop$condition)
pairwise.wilcox.test(loop$loop_area_mean,loop$condition, p.adjust="bonferroni")

p6 <- ggplot(loop,aes(x=condition, y=loop_area_mean, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("") + labs(y=bquote("µm"^2))  + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list(c("control", "severe"),c("mild", "severe"), c("moderate", "severe")))


ggsave(file="/path/to/figures/SupFig1F_VZ_loop_area.pdf",p6, width = 5.5, height = 5)


#Ac-TUB - Figure 3C - code by MG

wilcox <- org_means %>%
  group_by(Severity) %>%
  summarize(
    wilcox_p = wilcox.test(VZ_mean, MZ_mean)$p.value,
    n        = n()
  )

wilcox

org_long <- org_means %>%
  pivot_longer(
    cols = c(VZ_mean, MZ_mean),
    names_to  = "Type",
    values_to = "Density"
  )

summary_plot <- org_long %>%
  group_by(Severity, Type) %>%
  summarize(
    mean_dens = mean(Density),
    sd_dens   = sd(Density),
    .groups   = "drop"
  )

wilcox <- org_means %>%
  group_by(Severity) %>%
  summarize(
    wilcox_p = wilcox.test(VZ_mean, MZ_mean)$p.value,
    .groups  = "drop"
  ) %>%
  mutate(
    sig = case_when(
      wilcox_p < 0.001 ~ "***",
      wilcox_p < 0.01  ~ "**",
      wilcox_p < 0.05  ~ "*",
      TRUE             ~ "ns"
    )
  )

label_pos <- summary_plot %>%
  group_by(Severity) %>%
  summarize(
    y_max = max(mean_dens + sd_dens) * 1.05  # 5% headroom
  ) %>%
  left_join(wilcox, by = "Severity")

p7_new <- ggplot(summary_plot, aes(x = Severity, y = mean_dens,
                                   fill = interaction(Type,Severity))) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_dens - sd_dens,
                    ymax = mean_dens + sd_dens),
                width = .2, position = position_dodge(.9)) +
  # add p‐value stars
  geom_text(
    data = label_pos,
    aes(x = Severity, y = y_max, label = sig),
    inherit.aes = FALSE,
    size = 6
  ) +
  theme_classic() +
  ylab("AC-TUB strand density (%)") +
  scale_fill_manual(values = c(
    "#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac",
    "#7F7F7F","#6EAA5E","#414354","#469536"
  )) +
  theme(
    text = element_text(size = 20, colour = "black"),
    legend.position = "none"
  )

## Compare the VZ signal across severity grades (green bars)
pairwise.wilcox.test(
  x               = org_means$VZ_mean,
  g               = org_means$Severity,
  p.adjust.method = "BH"
)

## Compare the VZ signal across severity grades (gray bars)

pairwise.wilcox.test(
  x               = org_means$MZ_mean,
  g               = org_means$Severity,
  p.adjust.method = "BH"
)
ggsave("/path/to/figures/actin-tubulin_fig3c.pdf",p7,width=6, height=4)


#N-Cad signal Figure 3E
ncad<- read_excel("niche_disruption_quantification.xlsx", 
                                                 sheet = "NCad_sum")
ncad_bar<- read_excel("niche_disruption_quantification.xlsx", 
                  sheet = "NCad_bar")

kruskal.test(ncad$ncad~ncad$condition)
pairwise.wilcox.test(ncad$ncad,ncad$condition,p.adjust.method = "bonferroni")

p8 <- ggplot(ncad,aes(x=condition, y=ncad, fill=condition)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  #ylim(0,65) +
  xlab("") + labs(y="apical NCAD signal (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_text(size=18, colour="black"), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#BFBFBF","#7F7F7F","#414354")) + 
  stat_compare_means(comparisons = list(c("control", "mild"), c("control", "moderate"), c("control", "severe"), c("mild", "moderate"),c("mild", "severe"), c("moderate", "severe")), method = "wilcox")

ggsave("/path/to/figures/ncad_signal_fig3e.pdf",p8,width = 5,height=4.5)

#dividing cells Figure 3K
div <- read_excel("plane_of_cell_division.xlsx", 
                  sheet = "div_sum")

div_bar <- read_excel("plane_of_cell_division.xlsx", 
                  sheet = "div_bar")


library(ggpattern)
p9 <-ggplot(data=div_bar, aes(x=condition, y=mean, fill=type)) +
  geom_bar(stat="identity", color="black",linewidth=1.5)  +
  geom_col_pattern(
    aes(condition, mean, pattern_fill = type, pattern_angle = type, pattern_spacing = type), 
    fill            = '#414354',
    colour          = 'black', 
    pattern_density = 0.49) + scale_pattern_spacing_discrete(range = c(0.015, 0.01)) + theme_classic() +
    scale_pattern_fill_manual(values = c("#BFBFBF","#BFBFBF","#BFBFBF"))+ylab("dividing aRG cells (%)") +
  theme(text = element_text(size = 20,colour="black"))

ggsave("/path/to/figures/dividing_3k.pdf", p9, width = 7,height = 6)
  
###### QUANTIFICATIONS RESCUE #######

loop_treat <- read_excel("loop_parameters.xlsx", 
                   sheet = "loop_treat_sum")

#Figures 4B&4G

loop_treat$inter <- paste0(loop_treat$condition," ",loop_treat$treat)

loop_CHIR <- loop_treat[loop_treat$treat=="DMSO"|loop_treat$treat=="CHIR",]
loop_CHIR$inter <- factor(loop_CHIR$inter,levels = c("control DMSO","control CHIR","mild DMSO","mild CHIR",
                                                       "moderate DMSO","moderate CHIR","severe DMSO","severe CHIR"))

loop_epo <- loop_treat[loop_treat$treat=="DMSO"|loop_treat$treat=="EpoD",]
loop_epo$inter <- factor(loop_epo$inter,levels = c("control DMSO","control EpoD","mild DMSO","mild EpoD",
                                                     "moderate DMSO","moderate EpoD","severe DMSO","severe EpoD"))


### CHIR ###

#VZ diameter
kruskal.test(loop_CHIR$VZ_dia_mean~loop_CHIR$inter)
wilcox.test(loop_CHIR$VZ_dia_mean[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$VZ_dia_mean[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$VZ_dia_mean[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$VZ_dia_mean[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

com_chir <- list( c("control DMSO", "control CHIR"), c("mild DMSO", "mild CHIR"), c("moderate DMSO", "moderate CHIR"), c("severe DMSO","severe CHIR"))

p13a <- ggplot(loop_CHIR,aes(x=inter, y=VZ_dia_mean, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ylim(0,80) +
  xlab("") + labs(y="VZ  diameter (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

ggsave("/path/to/figures/Figure_4G.pdf",p13a, width = 6, height =5)

#apical membrane length
kruskal.test(loop_CHIR$ap_mem~loop_CHIR$inter)
wilcox.test(loop_CHIR$ap_mem[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$ap_mem[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$ap_mem[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$ap_mem[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

p13b <- ggplot(loop_CHIR,aes(x=inter, y=ap_mem, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("apical membrane length") +
  xlab("") + labs(y="apical membrane length (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

#basal membrane length
kruskal.test(loop_CHIR$bas_mem~loop_CHIR$inter)
wilcox.test(loop_CHIR$bas_mem[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$bas_mem[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$bas_mem[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$bas_mem[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

p13c <- ggplot(loop_CHIR,aes(x=inter, y=bas_mem, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("basal membrane length") +
  xlab("") + labs(y="basal membrane length (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

#venricle-like zone
kruskal.test(loop_CHIR$ven_size~loop_CHIR$inter)
wilcox.test(loop_CHIR$ven_size[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$ven_size[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$ven_size[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$ven_size[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

p13d <- ggplot(loop_CHIR,aes(x=inter, y=ven_size, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("ventricle-like zone")+
  xlab("") + labs(y=bquote("ventricle-like area µm"^2))  + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

#total loop area
kruskal.test(loop_CHIR$tot_loop~loop_CHIR$inter)
wilcox.test(loop_CHIR$tot_loop[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$tot_loop[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$tot_loop[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$tot_loop[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

p13e <- ggplot(loop_CHIR,aes(x=inter, y=tot_loop, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("total loop area") +
  xlab("") + labs(y=bquote("total loop area µm"^2)) + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

#total VZ area
kruskal.test(loop_CHIR$loop_tis~loop_CHIR$inter)
wilcox.test(loop_CHIR$loop_tis[loop_CHIR$condition=="control"]~loop_CHIR$treat[loop_CHIR$condition=="control"])
wilcox.test(loop_CHIR$loop_tis[loop_CHIR$condition=="mild"]~loop_CHIR$treat[loop_CHIR$condition=="mild"])
wilcox.test(loop_CHIR$loop_tis[loop_CHIR$condition=="moderate"]~loop_CHIR$treat[loop_CHIR$condition=="moderate"])
wilcox.test(loop_CHIR$loop_tis[loop_CHIR$condition=="severe"]~loop_CHIR$treat[loop_CHIR$condition=="severe"])

p13f <- ggplot(loop_CHIR,aes(x=inter, y=loop_tis, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  #ylim(0,80) +
  xlab("") + labs(y=bquote("total VZ area µm"^2))  + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF")) + 
  stat_compare_means(comparisons = com_chir, method = "wilcox")

p13sup <- ggarrange(p13b,p13c,p13d,p13e,nrow=2,ncol=2,common.legend = T)

ggsave("/path/to/figures/sup4C.pdf",p13sup, width=14, height=12)

##### EpoD ####

com_epo <- list( c("control DMSO", "control EpoD"), c("mild DMSO", "mild EpoD"), c("moderate DMSO", "moderate EpoD"), c("severe DMSO","severe EpoD"))

#VZ diameter
kruskal.test(loop_epo$VZ_dia_mean~loop_epo$inter)
wilcox.test(loop_epo$VZ_dia_mean[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$VZ_dia_mean[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$VZ_dia_mean[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$VZ_dia_mean[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])


p14 <- ggplot(loop_epo,aes(x=inter, y=VZ_dia_mean, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ylim(0,90) +
  xlab("") + labs(y="VZ  diameter (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

ggsave("/path/to/figures/Figure_4B.pdf",p14, width = 6, height =4)

#apical membrane length
kruskal.test(loop_epo$ap_mem~loop_epo$inter)
wilcox.test(loop_epo$ap_mem[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$ap_mem[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$ap_mem[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$ap_mem[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])

p14b <- ggplot(loop_epo,aes(x=inter, y=ap_mem, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("apical membrane length")  +
  xlab("") + labs(y="apical membrane length (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

#basal membrane length
kruskal.test(loop_epo$bas_mem~loop_epo$inter)
wilcox.test(loop_epo$bas_mem[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$bas_mem[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$bas_mem[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$bas_mem[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])

p14c <- ggplot(loop_epo,aes(x=inter, y=bas_mem, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("basal membrane length")  +
  xlab("") + labs(y="basal membrane length (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

#ventricle-like zone
kruskal.test(loop_epo$ven_size~loop_epo$inter)
wilcox.test(loop_epo$ven_size[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$ven_size[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$ven_size[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$ven_size[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])

p14d <- ggplot(loop_epo,aes(x=inter, y=ven_size, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("ventricle-like area")  +
  xlab("") + labs(y="ventricle-like area (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

#total loop area
kruskal.test(loop_epo$tot_loop~loop_epo$inter)
wilcox.test(loop_epo$tot_loop[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$tot_loop[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$tot_loop[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$tot_loop[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])

p14e <- ggplot(loop_epo,aes(x=inter, y=tot_loop, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("total loop area")  +
  xlab("") + labs(y="total loop area (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

#total VZ area
loop_epo$loop_tis <- as.numeric(loop_epo$loop_tis)
kruskal.test(loop_epo$loop_tis~loop_epo$inter)
wilcox.test(loop_epo$loop_tis[loop_epo$condition=="control"]~loop_epo$treat[loop_epo$condition=="control"])
wilcox.test(loop_epo$loop_tis[loop_epo$condition=="mild"]~loop_epo$treat[loop_epo$condition=="mild"])
wilcox.test(loop_epo$loop_tis[loop_epo$condition=="moderate"]~loop_epo$treat[loop_epo$condition=="moderate"])
wilcox.test(loop_epo$loop_tis[loop_epo$condition=="severe"]~loop_epo$treat[loop_epo$condition=="severe"])

p14f <- ggplot(loop_epo,aes(x=inter, y=loop_tis, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ggtitle("total VZ area") +
  xlab("") + labs(y="total VZ area (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank()) +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox") 

p14sup <- ggarrange(p14b,p14c,p14d,p14e,nrow=2,ncol=2,common.legend = T)

ggsave("figures/Sup4B.pdf",p14sup,width=14, height=12)

#Figure 4A - AC TUB density
epod_plot <- read_excel("ac_tub.xlsx", 
                        sheet = "EpoD_sum_MZ")
epod <- read_excel("ac_tub.xlsx", 
                   sheet = "EpoD_all")

#basal
epod$cat<-factor(paste0(epod$condition,".",epod$treat))
kruskal.test(epod$MZ,epod$cat)
wilcox.test(epod$MZ[epod$condition=="control"]~epod$treat[epod$condition=="control"])
wilcox.test(epod$MZ[epod$condition=="mild"]~epod$treat[epod$condition=="mild"])
wilcox.test(epod$MZ[epod$condition=="moderate"]~epod$treat[epod$condition=="moderate"])
wilcox.test(epod$MZ[epod$condition=="severe"]~epod$treat[epod$condition=="severe"])


p11 <- ggplot(data=epod_plot, aes(x=condition, y=mean, fill=interaction(treat,condition))) +
  geom_bar(stat="identity", color="black",position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9))  + theme_classic() + ylab("AC-TUB basal density (%)") +
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) +
  theme(text = element_text(size = 20,colour="black"),legend.position="none")

ggsave("figures/Fig4A.pdf",p11,width=6,height=5)


#CHIR rescue - dividing cells
div_chir <- read_excel("plane_of_cell_division.xlsx", 
                    sheet = "div_CHIR")
chir_bar <- read_excel("plane_of_cell_division.xlsx", 
                       sheet = "div_chir_bar")

chisq.test(div_chir[div_chir$condition=="control",3:5])
chisq.test(div_chir[div_chir$condition=="mild",3:5])
chisq.test(div_chir[div_chir$condition=="moderate",3:5])
chisq.test(div_chir[div_chir$condition=="severe",3:5])

chir_bar$inter <-paste0(chir_bar$condition," ",chir_bar$treat)
chir_bar$inter <- factor(chir_bar$inter,levels = c("control DMSO","control CHIR","mild DMSO","mild CHIR",
                                                   "moderate DMSO","moderate CHIR","severe DMSO","severe CHIR"))

p12 <- ggplot(data=chir_bar, aes(x=inter, y=perc_div, fill=type)) +
  geom_bar(stat="identity", color="black",linewidth=1.5) +
  geom_col_pattern(
    aes(inter, perc_div, pattern_fill = inter, pattern_angle = type, pattern_spacing = type), 
    fill="grey50",
    colour          = 'black', 
    pattern_density = 0.49) + scale_pattern_spacing_discrete(range = c(0.015, 0.01)) + theme_classic() + ylab("dividing aRG cells (%)") +
  scale_pattern_fill_manual(values = c("#FFFFFF","#E5E5FF","#BFBFBF","#B2B2FF","#7F7F7F","#7F7FFF","#414354","#6666FF"))+
  theme(text = element_text(size = 20,colour="black"),axis.text.x = element_blank(),legend.position="none") + xlab("")

ggsave("figures/Figure4H.pdf",p12, width = 7,height = 5)


#NCAD signal - rescue

ncad_res <- read_excel("niche_disruption_quantification.xlsx", 
                       sheet = "resc_all")

ncad_res_bar <- read_excel("niche_disruption_quantification.xlsx", 
                       sheet = "NCad_bar_resc")

ncad_res$inter <- paste0(ncad_res$condition," ",ncad_res$treat)
ncad_res$inter <- factor(ncad_res$inter,levels = c("control DMSO","control EpoD","mild DMSO","mild EpoD",
                                                   "moderate DMSO","moderate EpoD","severe DMSO","severe EpoD"))

kruskal.test(ncad_res$ncad~ncad_res$inter)

#wilcoxon test within condition
wilcox.test(ncad_res$ncad[ncad_res$condition=="control"]~ncad_res$treat[ncad_res$condition=="control"])
wilcox.test(ncad_res$ncad[ncad_res$condition=="mild"]~ncad_res$treat[ncad_res$condition=="mild"])
wilcox.test(ncad_res$ncad[ncad_res$condition=="moderate"]~ncad_res$treat[ncad_res$condition=="moderate"])
wilcox.test(ncad_res$ncad[ncad_res$condition=="severe"]~ncad_res$treat[ncad_res$condition=="severe"])


p15a <- ggplot(ncad_res,aes(x=inter, y=ncad, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ylim(0,25) +
  xlab("") + labs(y="apical NCAD signal (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#DBEAD5","#BFBFBF","#B7D5Ac","#7F7F7F","#6EAA5E","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")

p15b <- ggplot(ncad_res,aes(x=inter, y=ncad, fill=inter)) +
  geom_boxplot() +
  geom_jitter(color="black", position=position_jitter(0.1)) +
  theme_classic() +
  ylim(0,25) +
  xlab("") + labs(y="apical NCAD signal (µm)") + 
  theme(text = element_text(size=18),axis.title.x=element_blank(), axis.text.y=element_text(size=18, colour="black"),axis.text.x=element_blank(), axis.ticks.x = element_blank(),legend.position = "none") +   
  scale_fill_manual(values=c("#FFFFFF","#469536","#BFBFBF","#469536","#7F7F7F","#469536","#414354","#469536")) + 
  stat_compare_means(comparisons = com_epo, method = "wilcox")


ggsave("figures/Figure_4D.pdf",p15a, width = 6, height =5)


