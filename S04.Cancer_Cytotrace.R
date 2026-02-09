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


      


