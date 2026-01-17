#加载包
#####################################
library(DoubletFinder)
library(Seurat)
library(dplyr)
library(stringr)
library(cowplot)
library(harmony)
library(patchwork)
library(data.table)
#####################################
setwd('I:\\肿瘤脑转移\\Experiment\\data')
#开始处理
###################################################  GSM5645891
GSM5645891.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE186344/GSM5645891/")
dense.size <- object.size(x = as.matrix(x = GSM5645891.data))
sparse.size <- object.size(x = GSM5645891.data)
GSM5645891 <- CreateAssayObject(counts = GSM5645891.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645891 <- CreateSeuratObject(GSM5645891, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645891@meta.data$sample<- "GSM5645891"
GSM5645891@meta.data$tissue<- "Breast1"
GSM5645891[["percent.mt"]] <- PercentageFeatureSet(object = GSM5645891, pattern = "^MT-")     #注意MT
GSM5645891 <- subset(x = GSM5645891, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM5645891 <- NormalizeData(object = GSM5645891, normalization.method = "LogNormalize", scale.factor = 10000)
GSM5645891 <- FindVariableFeatures(object = GSM5645891, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM5645891)
GSM5645891 <- ScaleData(GSM5645891, features = all.genes)
GSM5645891 <- RunPCA(object= GSM5645891,npcs = 20,pc.genes=VariableFeatures(object = GSM5645891))  #这里取前20个PCs
pcSelect=20
GSM5645891 <- FindNeighbors(object = GSM5645891, dims = 1:pcSelect) 
GSM5645891 <- FindClusters(object = GSM5645891, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM5645891 <- RunUMAP(GSM5645891, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM5645891, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM5645891@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM5645891@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM5645891 <- doubletFinder_v3(GSM5645891, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM5645891 <- doubletFinder_v3(GSM5645891, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM5645891@meta.data$Doublet <- GSM5645891@meta.data$DF.classifications_0.25_0.08_198
table(GSM5645891@meta.data$DF.classifications_0.25_0.08_198)
#保存
saveRDS(GSM5645891, file = "./GSM5645891.rds")
#删除占内存的
rm(list = ls())
gc()

###############################################################
###################################################  GSM5645892
GSM5645892.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE186344/GSM5645892/")
dense.size <- object.size(x = as.matrix(x = GSM5645892.data))
sparse.size <- object.size(x = GSM5645892.data)
GSM5645892 <- CreateAssayObject(counts = GSM5645892.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645892 <- CreateSeuratObject(GSM5645892, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645892@meta.data$sample<- "GSM5645892"
GSM5645892@meta.data$tissue<- "Breast2"
GSM5645892[["percent.mt"]] <- PercentageFeatureSet(object = GSM5645892, pattern = "^MT-")     #注意MT
GSM5645892 <- subset(x = GSM5645892, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM5645892 <- NormalizeData(object = GSM5645892, normalization.method = "LogNormalize", scale.factor = 10000)
GSM5645892 <- FindVariableFeatures(object = GSM5645892, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM5645892)
GSM5645892 <- ScaleData(GSM5645892, features = all.genes)
GSM5645892 <- RunPCA(object= GSM5645892,npcs = 20,pc.genes=VariableFeatures(object = GSM5645892))  #这里取前20个PCs
pcSelect=20
GSM5645892 <- FindNeighbors(object = GSM5645892, dims = 1:pcSelect) 
GSM5645892 <- FindClusters(object = GSM5645892, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM5645892 <- RunUMAP(GSM5645892, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM5645892, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM5645892@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM5645892@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM5645892 <- doubletFinder_v3(GSM5645892, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM5645892 <- doubletFinder_v3(GSM5645892, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM5645892@meta.data$Doublet <- GSM5645892@meta.data$DF.classifications_0.25_0.08_198
table(GSM5645892@meta.data$DF.classifications_0.25_0.3_276)
#保存
saveRDS(GSM5645892, file = "./GSM5645892.rds")
#删除占内存的
rm(list = ls())
gc()
####################
###################################################  GSM5645893
GSM5645893.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE186344/GSM5645893/")
dense.size <- object.size(x = as.matrix(x = GSM5645893.data))
sparse.size <- object.size(x = GSM5645893.data)
GSM5645893 <- CreateAssayObject(counts = GSM5645893.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645893 <- CreateSeuratObject(GSM5645893, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645893@meta.data$sample<- "GSM5645893"
GSM5645893@meta.data$tissue<- "Breast2"
GSM5645893[["percent.mt"]] <- PercentageFeatureSet(object = GSM5645893, pattern = "^MT-")     #注意MT
GSM5645893 <- subset(x = GSM5645893, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM5645893 <- NormalizeData(object = GSM5645893, normalization.method = "LogNormalize", scale.factor = 10000)
GSM5645893 <- FindVariableFeatures(object = GSM5645893, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM5645893)
GSM5645893 <- ScaleData(GSM5645893, features = all.genes)
GSM5645893 <- RunPCA(object= GSM5645893,npcs = 20,pc.genes=VariableFeatures(object = GSM5645893))  #这里取前20个PCs
pcSelect=20
GSM5645893 <- FindNeighbors(object = GSM5645893, dims = 1:pcSelect) 
GSM5645893 <- FindClusters(object = GSM5645893, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM5645893 <- RunUMAP(GSM5645893, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM5645893, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM5645893@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM5645893@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM5645893 <- doubletFinder_v3(GSM5645893, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM5645893 <- doubletFinder_v3(GSM5645893, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM5645893@meta.data$Doublet <- GSM5645893@meta.data$DF.classifications_0.25_0.07_298
table(GSM5645893@meta.data$DF.classifications_0.25_0.07_298)
#保存
saveRDS(GSM5645893, file = "./GSM5645893.rds")
#删除占内存的
rm(list = ls())
gc()

###################################################  GSM5645894
GSM5645894.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE186344/GSM5645894/")
dense.size <- object.size(x = as.matrix(x = GSM5645894.data))
sparse.size <- object.size(x = GSM5645894.data)
GSM5645894 <- CreateAssayObject(counts = GSM5645894.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645894 <- CreateSeuratObject(GSM5645894, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645894@meta.data$sample<- "GSM5645894"
GSM5645894@meta.data$tissue<- "Breast2"
GSM5645894[["percent.mt"]] <- PercentageFeatureSet(object = GSM5645894, pattern = "^MT-")     #注意MT
GSM5645894 <- subset(x = GSM5645894, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM5645894 <- NormalizeData(object = GSM5645894, normalization.method = "LogNormalize", scale.factor = 10000)
GSM5645894 <- FindVariableFeatures(object = GSM5645894, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM5645894)
GSM5645894 <- ScaleData(GSM5645894, features = all.genes)
GSM5645894 <- RunPCA(object= GSM5645894,npcs = 20,pc.genes=VariableFeatures(object = GSM5645894))  #这里取前20个PCs
pcSelect=20
GSM5645894 <- FindNeighbors(object = GSM5645894, dims = 1:pcSelect) 
GSM5645894 <- FindClusters(object = GSM5645894, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM5645894 <- RunUMAP(GSM5645894, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM5645894, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM5645894@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM5645894@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM5645894 <- doubletFinder_v3(GSM5645894, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM5645894 <- doubletFinder_v3(GSM5645894, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM5645894@meta.data$Doublet <- GSM5645894@meta.data$DF.classifications_0.25_0.14_304
table(GSM5645894@meta.data$DF.classifications_0.25_0.14_304)
#保存
saveRDS(GSM5645894, file = "./GSM5645894.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  GSM5645895
GSM5645895.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE186344/GSM5645895/")
dense.size <- object.size(x = as.matrix(x = GSM5645895.data))
sparse.size <- object.size(x = GSM5645895.data)
GSM5645895 <- CreateAssayObject(counts = GSM5645895.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645895 <- CreateSeuratObject(GSM5645895, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM5645895@meta.data$sample<- "GSM5645895"
GSM5645895@meta.data$tissue<- "Breast2"
GSM5645895[["percent.mt"]] <- PercentageFeatureSet(object = GSM5645895, pattern = "^MT-")     #注意MT
GSM5645895 <- subset(x = GSM5645895, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM5645895 <- NormalizeData(object = GSM5645895, normalization.method = "LogNormalize", scale.factor = 10000)
GSM5645895 <- FindVariableFeatures(object = GSM5645895, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM5645895)
GSM5645895 <- ScaleData(GSM5645895, features = all.genes)
GSM5645895 <- RunPCA(object= GSM5645895,npcs = 20,pc.genes=VariableFeatures(object = GSM5645895))  #这里取前20个PCs
pcSelect=20
GSM5645895 <- FindNeighbors(object = GSM5645895, dims = 1:pcSelect) 
GSM5645895 <- FindClusters(object = GSM5645895, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM5645895 <- RunUMAP(GSM5645895, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM5645895, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM5645895@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM5645895@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM5645895 <- doubletFinder_v3(GSM5645895, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM5645895 <- doubletFinder_v3(GSM5645895, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM5645895@meta.data$Doublet <- GSM5645895@meta.data$DF.classifications_0.25_0.3_93
table(GSM5645895@meta.data$DF.classifications_0.25_0.3_93)
#保存
saveRDS(GSM5645895, file = "./GSM5645895.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  GSM7475325
GSM7475325.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE234832/GSM7475325/")
dense.size <- object.size(x = as.matrix(x = GSM7475325.data))
sparse.size <- object.size(x = GSM7475325.data)
GSM7475325 <- CreateAssayObject(counts = GSM7475325.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475325 <- CreateSeuratObject(GSM7475325, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475325@meta.data$sample<- "GSM7475325"
GSM7475325@meta.data$tissue<- "Breast2"
GSM7475325[["percent.mt"]] <- PercentageFeatureSet(object = GSM7475325, pattern = "^MT-")     #注意MT
GSM7475325 <- subset(x = GSM7475325, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM7475325 <- NormalizeData(object = GSM7475325, normalization.method = "LogNormalize", scale.factor = 10000)
GSM7475325 <- FindVariableFeatures(object = GSM7475325, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM7475325)
GSM7475325 <- ScaleData(GSM7475325, features = all.genes)
GSM7475325 <- RunPCA(object= GSM7475325,npcs = 20,pc.genes=VariableFeatures(object = GSM7475325))  #这里取前20个PCs
pcSelect=20
GSM7475325 <- FindNeighbors(object = GSM7475325, dims = 1:pcSelect) 
GSM7475325 <- FindClusters(object = GSM7475325, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM7475325 <- RunUMAP(GSM7475325, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM7475325, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM7475325@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM7475325@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM7475325 <- doubletFinder_v3(GSM7475325, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM7475325 <- doubletFinder_v3(GSM7475325, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM7475325@meta.data$Doublet <- GSM7475325@meta.data$DF.classifications_0.25_0.3_16
table(GSM7475325@meta.data$DF.classifications_0.25_0.3_16)
#保存
saveRDS(GSM7475325, file = "./GSM7475325.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  GSM7475326
GSM7475326.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE234832/GSM7475326/")
dense.size <- object.size(x = as.matrix(x = GSM7475326.data))
sparse.size <- object.size(x = GSM7475326.data)
GSM7475326 <- CreateAssayObject(counts = GSM7475326.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475326 <- CreateSeuratObject(GSM7475326, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475326@meta.data$sample<- "GSM7475326"
GSM7475326@meta.data$tissue<- "Breast2"
GSM7475326[["percent.mt"]] <- PercentageFeatureSet(object = GSM7475326, pattern = "^MT-")     #注意MT
GSM7475326 <- subset(x = GSM7475326, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM7475326 <- NormalizeData(object = GSM7475326, normalization.method = "LogNormalize", scale.factor = 10000)
GSM7475326 <- FindVariableFeatures(object = GSM7475326, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM7475326)
GSM7475326 <- ScaleData(GSM7475326, features = all.genes)
GSM7475326 <- RunPCA(object= GSM7475326,npcs = 20,pc.genes=VariableFeatures(object = GSM7475326))  #这里取前20个PCs
pcSelect=20
GSM7475326 <- FindNeighbors(object = GSM7475326, dims = 1:pcSelect) 
GSM7475326 <- FindClusters(object = GSM7475326, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM7475326 <- RunUMAP(GSM7475326, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM7475326, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM7475326@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM7475326@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM7475326 <- doubletFinder_v3(GSM7475326, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM7475326 <- doubletFinder_v3(GSM7475326, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM7475326@meta.data$Doublet <- GSM7475326@meta.data$DF.classifications_0.25_0.25_11
table(GSM7475326@meta.data$DF.classifications_0.25_0.25_11)
#保存
saveRDS(GSM7475326, file = "./GSM7475326.rds")
#删除占内存的
rm(list = ls())
gc()

####################################
###################################################  GSM7475327
GSM7475327.data <- Read10X(data.dir = "I:/肿瘤脑转移/GSE234832/GSM7475327/")
dense.size <- object.size(x = as.matrix(x = GSM7475327.data))
sparse.size <- object.size(x = GSM7475327.data)
GSM7475327 <- CreateAssayObject(counts = GSM7475327.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475327 <- CreateSeuratObject(GSM7475327, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
GSM7475327@meta.data$sample<- "GSM7475327"
GSM7475327@meta.data$tissue<- "Breast2"
GSM7475327[["percent.mt"]] <- PercentageFeatureSet(object = GSM7475327, pattern = "^MT-")     #注意MT
GSM7475327 <- subset(x = GSM7475327, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
GSM7475327 <- NormalizeData(object = GSM7475327, normalization.method = "LogNormalize", scale.factor = 10000)
GSM7475327 <- FindVariableFeatures(object = GSM7475327, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(GSM7475327)
GSM7475327 <- ScaleData(GSM7475327, features = all.genes)
GSM7475327 <- RunPCA(object= GSM7475327,npcs = 20,pc.genes=VariableFeatures(object = GSM7475327))  #这里取前20个PCs
pcSelect=20
GSM7475327 <- FindNeighbors(object = GSM7475327, dims = 1:pcSelect) 
GSM7475327 <- FindClusters(object = GSM7475327, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
GSM7475327 <- RunUMAP(GSM7475327, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(GSM7475327, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- GSM7475327@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(GSM7475327@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
GSM7475327 <- doubletFinder_v3(GSM7475327, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
GSM7475327 <- doubletFinder_v3(GSM7475327, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
GSM7475327@meta.data$Doublet <- GSM7475327@meta.data$DF.classifications_0.25_0.08_198
table(GSM7475327@meta.data$DF.classifications_0.25_0.3_36)
#保存
saveRDS(GSM7475327, file = "./GSM7475327.rds")
#删除占内存的
rm(list = ls())
gc()

##########################################
###################################################  LMBT_20250714
LMBT_20250714.data <- Read10X(data.dir = "I:/肿瘤脑转移/自测LCBM/LMBT_20250714/")
dense.size <- object.size(x = as.matrix(x = LMBT_20250714.data))
sparse.size <- object.size(x = LMBT_20250714.data)
LMBT_20250714 <- CreateAssayObject(counts = LMBT_20250714.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250714 <- CreateSeuratObject(LMBT_20250714, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250714@meta.data$sample<- "LMBT_20250714"
LMBT_20250714@meta.data$tissue<- "Breast2"
LMBT_20250714[["percent.mt"]] <- PercentageFeatureSet(object = LMBT_20250714, pattern = "^MT-")     #注意MT
LMBT_20250714 <- subset(x = LMBT_20250714, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
LMBT_20250714 <- NormalizeData(object = LMBT_20250714, normalization.method = "LogNormalize", scale.factor = 10000)
LMBT_20250714 <- FindVariableFeatures(object = LMBT_20250714, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(LMBT_20250714)
LMBT_20250714 <- ScaleData(LMBT_20250714, features = all.genes)
LMBT_20250714 <- RunPCA(object= LMBT_20250714,npcs = 20,pc.genes=VariableFeatures(object = LMBT_20250714))  #这里取前20个PCs
pcSelect=20
LMBT_20250714 <- FindNeighbors(object = LMBT_20250714, dims = 1:pcSelect) 
LMBT_20250714 <- FindClusters(object = LMBT_20250714, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
LMBT_20250714 <- RunUMAP(LMBT_20250714, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(LMBT_20250714, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- LMBT_20250714@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(LMBT_20250714@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
LMBT_20250714 <- doubletFinder_v3(LMBT_20250714, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
LMBT_20250714 <- doubletFinder_v3(LMBT_20250714, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
LMBT_20250714@meta.data$Doublet <- LMBT_20250714@meta.data$DF.classifications_0.25_0.3_337
table(LMBT_20250714@meta.data$DF.classifications_0.25_0.3_337)
#保存
saveRDS(LMBT_20250714, file = "./LMBT_20250714.rds")
#删除占内存的
rm(list = ls())
gc()

###################################################  LMBT_20250722
LMBT_20250722.data <- Read10X(data.dir = "I:/肿瘤脑转移/自测LCBM/LMBT_20250722/")
dense.size <- object.size(x = as.matrix(x = LMBT_20250722.data))
sparse.size <- object.size(x = LMBT_20250722.data)
LMBT_20250722 <- CreateAssayObject(counts = LMBT_20250722.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250722 <- CreateSeuratObject(LMBT_20250722, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250722@meta.data$sample<- "LMBT_20250722"
LMBT_20250722@meta.data$tissue<- "Breast2"
LMBT_20250722[["percent.mt"]] <- PercentageFeatureSet(object = LMBT_20250722, pattern = "^MT-")     #注意MT
LMBT_20250722 <- subset(x = LMBT_20250722, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
LMBT_20250722 <- NormalizeData(object = LMBT_20250722, normalization.method = "LogNormalize", scale.factor = 10000)
LMBT_20250722 <- FindVariableFeatures(object = LMBT_20250722, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(LMBT_20250722)
LMBT_20250722 <- ScaleData(LMBT_20250722, features = all.genes)
LMBT_20250722 <- RunPCA(object= LMBT_20250722,npcs = 20,pc.genes=VariableFeatures(object = LMBT_20250722))  #这里取前20个PCs
pcSelect=20
LMBT_20250722 <- FindNeighbors(object = LMBT_20250722, dims = 1:pcSelect) 
LMBT_20250722 <- FindClusters(object = LMBT_20250722, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
LMBT_20250722 <- RunUMAP(LMBT_20250722, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(LMBT_20250722, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- LMBT_20250722@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(LMBT_20250722@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
LMBT_20250722 <- doubletFinder_v3(LMBT_20250722, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
LMBT_20250722 <- doubletFinder_v3(LMBT_20250722, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
LMBT_20250722@meta.data$Doublet <- LMBT_20250722@meta.data$DF.classifications_0.25_0.27_318
table(LMBT_20250722@meta.data$DF.classifications_0.25_0.27_318)
#保存
saveRDS(LMBT_20250722, file = "./LMBT_20250722.rds")
#删除占内存的
rm(list = ls())
gc()

###################################################  LMBT_20250724
LMBT_20250724.data <- Read10X(data.dir = "I:/肿瘤脑转移/自测LCBM/LMBT_20250724/")
dense.size <- object.size(x = as.matrix(x = LMBT_20250724.data))
sparse.size <- object.size(x = LMBT_20250724.data)
LMBT_20250724 <- CreateAssayObject(counts = LMBT_20250724.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250724 <- CreateSeuratObject(LMBT_20250724, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250724@meta.data$sample<- "LMBT_20250724"
LMBT_20250724@meta.data$tissue<- "Breast2"
LMBT_20250724[["percent.mt"]] <- PercentageFeatureSet(object = LMBT_20250724, pattern = "^MT-")     #注意MT
LMBT_20250724 <- subset(x = LMBT_20250724, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
LMBT_20250724 <- NormalizeData(object = LMBT_20250724, normalization.method = "LogNormalize", scale.factor = 10000)
LMBT_20250724 <- FindVariableFeatures(object = LMBT_20250724, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(LMBT_20250724)
LMBT_20250724 <- ScaleData(LMBT_20250724, features = all.genes)
LMBT_20250724 <- RunPCA(object= LMBT_20250724,npcs = 20,pc.genes=VariableFeatures(object = LMBT_20250724))  #这里取前20个PCs
pcSelect=20
LMBT_20250724 <- FindNeighbors(object = LMBT_20250724, dims = 1:pcSelect) 
LMBT_20250724 <- FindClusters(object = LMBT_20250724, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
LMBT_20250724 <- RunUMAP(LMBT_20250724, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(LMBT_20250724, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- LMBT_20250724@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(LMBT_20250724@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
LMBT_20250724 <- doubletFinder_v3(LMBT_20250724, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
LMBT_20250724 <- doubletFinder_v3(LMBT_20250724, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
LMBT_20250724@meta.data$Doublet <- LMBT_20250724@meta.data$DF.classifications_0.25_0.08_198
table(LMBT_20250724@meta.data$DF.classifications_0.25_0.005_197)
#保存
saveRDS(LMBT_20250724, file = "./LMBT_20250724.rds")
#删除占内存的
rm(list = ls())
gc()


###################################################  LMBT_20250820
LMBT_20250820.data <- Read10X(data.dir = "I:/肿瘤脑转移/自测LCBM/LMBT_20250820/")
dense.size <- object.size(x = as.matrix(x = LMBT_20250820.data))
sparse.size <- object.size(x = LMBT_20250820.data)
LMBT_20250820 <- CreateAssayObject(counts = LMBT_20250820.data, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250820 <- CreateSeuratObject(LMBT_20250820, min.cells = 3, min.features = 200, names.delim = "_", project = "10X_GC")
LMBT_20250820@meta.data$sample<- "LMBT_20250820"
LMBT_20250820@meta.data$tissue<- "Breast2"
LMBT_20250820[["percent.mt"]] <- PercentageFeatureSet(object = LMBT_20250820, pattern = "^MT-")     #注意MT
LMBT_20250820 <- subset(x = LMBT_20250820, subset = nFeature_RNA > 200 & percent.mt < 10)     #这里取了线粒体10%
LMBT_20250820 <- NormalizeData(object = LMBT_20250820, normalization.method = "LogNormalize", scale.factor = 10000)
LMBT_20250820 <- FindVariableFeatures(object = LMBT_20250820, selection.method = "vst", nfeatures = 2500)   #这里设置为前2500
all.genes <- rownames(LMBT_20250820)
LMBT_20250820 <- ScaleData(LMBT_20250820, features = all.genes)
LMBT_20250820 <- RunPCA(object= LMBT_20250820,npcs = 20,pc.genes=VariableFeatures(object = LMBT_20250820))  #这里取前20个PCs
pcSelect=20
LMBT_20250820 <- FindNeighbors(object = LMBT_20250820, dims = 1:pcSelect) 
LMBT_20250820 <- FindClusters(object = LMBT_20250820, resolution = 0.4)    #注意每个样本的resolution的参数都不一样，需要根据细胞数量设定
LMBT_20250820 <- RunUMAP(LMBT_20250820, dims = 1:pcSelect) #UMAP
#去多胞
sweep.res.list_SM <- paramSweep_v3(LMBT_20250820, PCs = 1:20)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM, GT = FALSE)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- LMBT_20250820@meta.data$res.0.6
homotypic.prop <- modelHomotypic(annotations)  
nExp_poi <- round(0.031*length(LMBT_20250820@meta.data$orig.ident)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
LMBT_20250820 <- doubletFinder_v3(LMBT_20250820, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE)
LMBT_20250820 <- doubletFinder_v3(LMBT_20250820, PCs = 1:20, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value)	 
LMBT_20250820@meta.data$Doublet <- LMBT_20250820@meta.data$DF.classifications_0.25_0.24_237
table(LMBT_20250820@meta.data$DF.classifications_0.25_0.24_237)
#保存
saveRDS(LMBT_20250820, file = "./LMBT_20250820.rds")
#删除占内存的
rm(list = ls())
gc()

#########################################################


