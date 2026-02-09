######################### Monocel
library(monocle)
################################
#library(ClusterGVis)
#####################
library(Seurat)
library(dplyr)
library(stringr)
#library(harmony)
library(patchwork)
library(ggplot2)
library(ggbump)
library(knitr)
#######################
suppressMessages(library(scMEGA))
suppressMessages(library(ArchR))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(harmony))
suppressMessages(library(dplyr))
suppressMessages(library(cowplot))
suppressMessages(library(harmony))
suppressMessages(library(Nebulosa))
suppressMessages(library(ggpubr))
#suppressMessages(library(Ipaper))
suppressMessages(library(phateR))
suppressMessages(library(viridis))
#suppressMessages(library(SeuratDisk))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(destiny))
suppressMessages(library(plotly))
######################
setwd('I:\\肿瘤脑转移\\Experiment\\data')
pbmc <- readRDS('Cytotrace_pbmc.rds')

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))   #这里PCA设置为前100个

#开始去批次
pbmc$orig.ident <- pbmc$sample
pbmc <- RunHarmony(pbmc,"orig.ident",dims.use = 1:20,max.iter.harmony = 15,max.iter.cluster = 20)   #去批次这一部看结果
pbmc <- RunUMAP(pbmc, dims = 1:20, min.dist = 0.2,reduction = "harmony")   #尝试调整min.dist会得到更好的结果
pbmc <- FindNeighbors(pbmc, reduction = "harmony",dims = 1:20)

pbmc <- FindClusters(pbmc, resolution = 0.4)   


annotation <- data.frame(phenotype=pbmc@meta.data$seurat_clusters) %>% set_rownames(., colnames(pbmc))
plots <- plotData(cytotrace2_result = pbmc, 
                  annotation = annotation,
                  expression_data = pbmc,is_seurat = TRUE
                  )

#pbmc <- RunTSNE(pbmc, dims = 1:20, tsne.method = "FIt-SNE", nthreads = 4, max_iter = 2000)    #尝试调整max_iter会得到更好的结果
pbmc <- RunUMAP(pbmc, dims = 1:20, min.dist = 0.75)
Idents(pbmc) <- pbmc$seurat_clusters
##################################################### Subcluster                                   
pbmc_subset <- subset(pbmc, idents = c(5,3,4,2,0))                                       

Idents(pbmc_subset) <- pbmc_subset$seurat_clusters
pbmc <- pbmc_subset

new.cluster.ids <- c("Differentiated", "Unipotent", "Multipotent", "Oligopotent", "Pluripotent")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$Diff <- Idents(pbmc)
########################
coembed <- pbmc
##############################
df <- coembed@meta.data %>%
    as.data.frame() %>%
    subset(., Diff %in% c("Differentiated", "Unipotent", "Multipotent", "Oligopotent", "Pluripotent"))

coembed <- coembed[, rownames(df)]
cols <- ArchR::paletteDiscrete(unique(coembed@meta.data[, "Diff"]))
###################################
traj1 <- c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent")
#######################
coembed <- RunUMAP(coembed, 
               dims = 1:20, 
               reduction = 'harmony',
               reduction.name = "harmony_v2",
               reduction.ke = 'harmony_v2',
              verbose = FALSE,
                    min.dist = 0.4)
####################################################
coembed <- AddTrajectory(object = coembed, 
                          trajectory = c("Differentiated", "Unipotent", "Oligopotent", "Multipotent", "Pluripotent"),
                          group.by = "Diff", 
                          reduction = "harmony_v2",
                          dims = 1:2, 
                          use.all = TRUE)                                   

coembed <- coembed[, !is.na(coembed$Trajectory)]

p1 <- DimPlot(coembed, reduction = "harmony_v2", 
              group.by = "Diff", cols = cols) +
              xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("Cluster")


p2 <- TrajectoryPlot(object = coembed, 
                    reduction = "harmony_v2",
                    continuousSet = "blueYellow",
                    size = 1,
                   addArrow = TRUE) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    ggtitle("Trajectory")
#########################
trajRNA <- GetTrajectory(coembed, assay = 'RNA', trajectory.name = "Trajectory", 
        groupEvery = 1, slot = "data", smoothWindow = 7, 
        log2Norm = TRUE)

groupMatRNA <- suppressMessages(TrajectoryHeatmap(trajRNA, 
        varCutOff = 0.9, pal = paletteContinuous(set = "horizonExtra"), 
        limits = c(-2, 2), returnMatrix = TRUE))

Gene_set_change <- rownames(groupMatRNA)
saveRDS(groupMatRNA, file = 'Change_genes.rds')
########################################
library(Mfuzz)
mfuzz_class <- new('ExpressionSet', exprs = groupMatRNA)

c <- 3
m <- mestimate(mfuzz_class)
cl <- mfuzz(mfuzz_class, c = c, m = m)

nrow <- ceiling(c/3)

mfuzz.plot(mfuzz_class,cl,mfrow=c(nrow,3),new.window= FALSE)
#mfuzz.plot(gene.s,cl,mfrow=c(2,3),new.window= FALSE,time.labels=as.vector(colnames(gene_tpm)))
mfuzz.plot2(mfuzz_class,cl,mfrow=c(nrow,3),time.labels=colnames(groupMatRNA),x11 = FALSE)

select_gene <- names(cl$cluster)[cl$cluster == 2]

setwd('I:\\肿瘤脑转移')
sati <- read.csv('./神经活性基因.csv')
library(homologene)
abc <- sati$gene
abc <- mouse2human(abc)


commen <- intersect(select_gene,abc$humanGene)


#########################################

