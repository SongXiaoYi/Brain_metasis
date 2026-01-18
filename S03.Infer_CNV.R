library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(ggplot2)
#######################
library(fastCNV)
#########################
setwd('I:\\肿瘤脑转移\\Experiment\\data')
sce.integrated <- readRDS('./pbmc.rds')

scColon1 <- fastCNV(seuratObj = sce.integrated, sampleName = "sce.integrated", referenceVar = "seurat_clusters", referenceLabel = c("6", "7", "8","9", "10","11"), printPlot = T)
scColon1 <- fastCNV(seuratObj = sce.integrated, sampleName = "sce.integrated", referenceVar = "seurat_clusters", referenceLabel = c("0"), printPlot = T)

common_theme <- theme(
  plot.title = element_text(size = 10),
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  axis.title = element_text(size = 8),
  axis.text = element_text(size = 6)
)
 
FeaturePlot(scColon1, features = "cnv_fraction", reduction = "umap" ) & common_theme |
  DimPlot(scColon1, reduction = "umap", group.by =  "seurat_clusters",label = TRUE) & common_theme


ggplot(FetchData(scColon1, vars = c("seurat_clusters", "cnv_fraction")), 
       aes(seurat_clusters, cnv_fraction, fill = seurat_clusters)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"))


#############################
library(scales)
 
FeaturePlot(scColon1, features = "20.p_CNV")  +
  scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1), 
                       rescaler = function(x, to = c(0, 1), from = NULL) {
                         rescale_mid(x, to = to, mid = 0)
                       }) +
  common_theme |
FeaturePlot(scColon1, features = "X.q_CNV") +
  scale_color_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1), 
                       rescaler = function(x, to = c(0, 1), from = NULL) {
                         rescale_mid(x, to = to, mid = 0)
                       }) +
  common_theme


############################################
scColon1 <- CNVClassification(scColon1)
 
DimPlot(scColon1, group.by = "20.p_CNV_classification") &
  scale_color_manual(values = c(gain = "red", no_alteration = "grey", loss = "blue")) &
  common_theme |
DimPlot(scColon1, group.by = "X.q_CNV_classification") &
  scale_color_manual(values = c(gain = "red", no_alteration = "grey", loss = "blue")) &
  common_theme

##########################################
DimPlot(scColon1, group.by = "cnv_clusters") + common_theme

###############################################
abcs[abcs$Clusters %in% c('1','2','3','7','8','9','11','12','18','19','21','23'),'Cancer_type'] <- 'Cancer'
scColon1$Cancer_annotation <- abcs$Cancer_type

DimPlot(scColon1, group.by = "Cancer_annotation")
saveRDS(scColon1,'./pbmc_cnv_rather.rds')
