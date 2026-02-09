library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(ggplot2)
library(CytoTRACE2)
####################
setwd('I:\\肿瘤脑转移\\Experiment\\data')
pbmc <- readRDS('./pbmc_cnv_rather.rds')

############################
Idents(pbmc) <- pbmc$Cancer_annotation                                       
pbmc <- subset(pbmc, idents = c('Cancer'))                                       

#Plasma <- subset(pbmc,idents = 'PC')
Plasma <- pbmc
meta <- Plasma@meta.data
                                       
Plasma.CytoTRACE <- cytotrace2(Plasma,species = "human",is_seurat = TRUE)
meta$CytoTRACE2_Score <- Plasma.CytoTRACE$CytoTRACE2_Score
meta$CytoTRACE2_Potency <- Plasma.CytoTRACE$CytoTRACE2_Potency
meta$CytoTRACE2_Relative <- Plasma.CytoTRACE$CytoTRACE2_Relative
meta$preKNN_CytoTRACE2_Score <- Plasma.CytoTRACE$preKNN_CytoTRACE2_Score
meta$preKNN_CytoTRACE2_Potency <- Plasma.CytoTRACE$preKNN_CytoTRACE2_Potency
annotation <- data.frame(phenotype=Plasma.CytoTRACE@meta.data$Cancer_annotation) %>% set_rownames(., colnames(Plasma.CytoTRACE))
plots <- plotData(cytotrace2_result = Plasma.CytoTRACE, 
                  annotation = annotation,
                  expression_data = Plasma.CytoTRACE,is_seurat = TRUE
                  )

############################################## Barplot
Plot_matrix <- meta
Plot_matrix$group <- factor(Plot_matrix$group,levels = c("T","L","M"))
my_comparisons <- list(c("T", "L"),c("T","M"),c("M","L"))
p1 <- ggplot(Plot_matrix,aes(group, CytoTRACE2_Score)) +   ggsignif::geom_signif(
    color = "black",
    comparisons = list(c("T", "L"),c("T", "M"),c("L", "M")),
    test = wilcox.test,
    step_increase = -0.2,
    textsize = 0.5*6
  )+
     geom_boxplot(aes(fill = group),outlier.shape = NA,alpha = 1)+
     theme_classic()+ylab("CytoTRACE score") + ylim(0,0.5) + xlab('')
setwd('G:\\YouRui\\abc')
pdf(file="./Tissue_CytoTRACE2.pdf",width=2.55,height=2.20)
p1
dev.off()                                       
                                       
###################################################
Plot_matrix <- meta
#Plot_matrix <- Plot_matrix[which(Plot_matrix$Isotypes != "None"),]
#Plot_matrix$Isotypes <- factor(Plot_matrix$Isotypes,levels=c("IgM","IgD","IgG4","IgA2","IgG2","IgA1","IgG1","IgG3"))
p2 <- ggplot(Plot_matrix,aes(B_SubCluster, CytoTRACE2_Score)) + stat_compare_means()+
     geom_boxplot(aes(fill = B_SubCluster),outlier.shape = NA,alpha = 1)+
     theme_classic()+ylab("CytoTRACE score") + ylim(0,0.3) + xlab('') + theme(axis.text.x = element_text(angle = 90, size = 10,vjust = 0.5,hjust = 1)) + NoLegend()
      


