#加载包
#####################################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(ggplot2)
#####################################
setwd('I:\\肿瘤脑转移\\Experiment\\data')
##################
GSM5645891_Breast1 <- readRDS("./GSM5645891.rds")
GSM5645892_Breast2 <- readRDS("./GSM5645892.rds")
GSM5645893_Breast3 <- readRDS("./GSM5645893.rds")
GSM5645894_Lung1 <- readRDS("./GSM5645894.rds")
GSM5645895_Lung2 <- readRDS("./GSM5645895.rds")
GSM7475325_BRBMET_3 <- readRDS("./GSM7475325.rds")
GSM7475326_BRBMET_87 <- readRDS("./GSM7475326.rds")
GSM7475327_LUBMET_7 <- readRDS("./GSM7475327.rds")
LMBT_20250714 <- readRDS("./LMBT_20250714.rds")
LMBT_20250722 <- readRDS("./LMBT_20250722.rds")
LMBT_20250724 <- readRDS("./LMBT_20250724.rds")
LMBT_20250820 <- readRDS("./LMBT_20250820.rds")



GSM5645891_Breast1 <- subset(GSM5645891_Breast1,subset = Doublet == "Singlet")
GSM5645892_Breast2 <- subset(GSM5645892_Breast2,subset = Doublet == "Singlet")
GSM5645893_Breast3 <- subset(GSM5645893_Breast3,subset = Doublet == "Singlet")
GSM5645894_Lung1 <- subset(GSM5645894_Lung1,subset = Doublet == "Singlet")
GSM5645895_Lung2 <- subset(GSM5645895_Lung2,subset = Doublet == "Singlet")
GSM7475325_BRBMET_3 <- subset(GSM7475325_BRBMET_3,subset = Doublet == "Singlet")
GSM7475326_BRBMET_87 <- subset(GSM7475326_BRBMET_87,subset = Doublet == "Singlet")
GSM7475327_LUBMET_7 <- subset(GSM7475327_LUBMET_7,subset = Doublet == "Singlet")
LMBT_20250714 <- subset(LMBT_20250714,subset = Doublet == "Singlet")
LMBT_20250722 <- subset(LMBT_20250722,subset = Doublet == "Singlet")
LMBT_20250724 <- subset(LMBT_20250724,subset = Doublet == "Singlet")
LMBT_20250820 <- subset(LMBT_20250820,subset = Doublet == "Singlet")


cellid  <- c('GSM5645891_Breast1','GSM5645892_Breast2','GSM5645893_Breast3','GSM5645894_Lung1','GSM5645895_Lung2','GSM7475325_BRBMET_3','GSM7475326_BRBMET_87','GSM7475327_LUBMET_7'
            ,'LMBT_20250714','LMBT_20250722','LMBT_20250724','LMBT_20250820')
sceList <- list(GSM5645892_Breast2,GSM5645893_Breast3,GSM5645894_Lung1,GSM5645895_Lung2,GSM7475325_BRBMET_3,GSM7475326_BRBMET_87,GSM7475327_LUBMET_7
            ,LMBT_20250714,LMBT_20250722,LMBT_20250724,LMBT_20250820)

pbmc.big <- merge(GSM5645891_Breast1, y = sceList, add.cell.ids = cellid, project = "10X_GC")
saveRDS(pbmc.big, file = "sce.anchors.rds")

rm(list = ls())
gc()

#################################
sce.integrated <- readRDS('./sce.anchors.rds')
sce.integrated <- NormalizeData(object = sce.integrated, normalization.method = "LogNormalize", scale.factor = 10000)
sce.integrated <- FindVariableFeatures(object = sce.integrated, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(sce.integrated)
sce.integrated <- ScaleData(sce.integrated, features = all.genes)
sce.integrated <- RunPCA(object= sce.integrated,npcs = 20,pc.genes=VariableFeatures(object = sce.integrated))   #这里PCA设置为前100个

#开始去批次
sce.integrated$orig.ident <- sce.integrated$sample
sce.integrated <- RunHarmony(sce.integrated,"orig.ident",dims.use = 1:20,max.iter.harmony = 15,max.iter.cluster = 20)   #去批次这一部看结果
sce.integrated <- RunUMAP(sce.integrated, dims = 1:20, min.dist = 0.2,reduction = "harmony")   #尝试调整min.dist会得到更好的结果
sce.integrated <- FindNeighbors(sce.integrated, reduction = "harmony",dims = 1:20)

sce.integrated <- FindClusters(sce.integrated, resolution = 0.4)   

#sce.integrated <- RunTSNE(sce.integrated, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)    #尝试调整max_iter会得到更好的结果
sce.integrated <- RunUMAP(sce.integrated, dims = 1:20, min.dist = 0.75)   #尝试调整min.dist会得到更好的结果

#画图
DimPlot(sce.integrated,reduction = "umap", label = T) + NoLegend()
saveRDS(sce.integrated,'./pbmc.rds')
