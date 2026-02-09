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
        groupEvery = 1, slot = "data", smoothWindow = 11, 
        log2Norm = TRUE)

groupMatRNA <- suppressMessages(TrajectoryHeatmap(trajRNA, 
        varCutOff = 0.9, pal = paletteContinuous(set = "horizonExtra"), 
        limits = c(-2, 2), returnMatrix = TRUE))

Gene_set_change <- rownames(groupMatRNA)
saveRDS(groupMatRNA, file = 'Change_genes.rds')
########################################
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

trajRNA <- trajRNA[select_gene[1000:1464],]
ht <- TrajectoryHeatmap(trajRNA,
                        varCutOff = 0.9,
                        pal = paletteContinuous(set = "horizonExtra"),
                        limits = c(-2, 2),labelMarkers = commen, labelTop = 25)
#########################################
library(Seurat)      # 用于单细胞数据分析
library(ggplot2)     # 用于数据可视化
library(ggpubr)      # 用于添加统计信息到图形
library(ggExtra)     # 用于添加边缘分布图
library(dplyr) 
####################################
####################################
plot_data1 <- FetchData(pbmc, 
                      vars = c("PGAP1", "SLIT2"),
                      slot = "data")
plot_data1_filtered <- plot_data1[plot_data1$PGAP1 > 0 & plot_data1$SLIT2 > 0, ]

options(repr.plot.height=8, repr.plot.width=8)  # 设置图形大小
p1_filtered <- ggplot(plot_data1_filtered, aes(x = PGAP1, y = SLIT2)) +
    geom_point(size = 1, alpha = 0.5, color = "grey30") +  # 添加散点
    geom_smooth(method = "lm", color = "blue", se = TRUE, fill = "grey90") +  # 添加拟合线和置信区间
    labs(x = "PGAP1 expression", y = "SLIT2 expression") +  # 设置轴标签
    # 设置主题和样式
    theme_bw() +
    theme(
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),  # 主网格线
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.3),  # 次网格线
        panel.background = element_rect(fill = "white"),  # 背景颜色
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # 边框
        axis.ticks = element_line(color = "black", linewidth = 0.3),  # 刻度线
        axis.text = element_text(color = "black", size = 10),  # 轴文本
        axis.title = element_text(color = "black", size = 16),  # 轴标题
        axis.ticks.length = unit(0.2, "cm"),  # 刻度线长度
        plot.background = element_rect(fill = "white")  # 图形背景
    ) +
    # 设置x轴刻度
    scale_x_continuous(
        breaks = seq(0, max(plot_data1_filtered$PGAP1), by = 1),
        minor_breaks = seq(0, max(plot_data1_filtered$PGAP1), by = 0.5)
    ) +
    # 设置y轴刻度
    scale_y_continuous(
        limits = c(0, max(plot_data1_filtered$SLIT2)),
        breaks = seq(0, max(plot_data1_filtered$SLIT2), by = 1),
        minor_breaks = seq(0, max(plot_data1_filtered$SLIT2), by = 0.5)
    ) +
    # 添加相关性统计信息
    stat_cor(method = "pearson",
             label.x.npc = "left",
             label.y.npc = "top",
             size = 4)

# 为过滤后的散点图添加边缘分布图
p1_filtered_with_marginal <- ggMarginal(p1_filtered,
                             type = "density",  # 设置边缘图类型为密度图
                             margins = "both",  # 显示两个轴的边缘分布
                             size = 5,
                             xparams = list(fill = "orange", alpha = 1),  # x轴边缘图样式
                             yparams = list(fill = "blue", alpha = 1))    # y轴边缘图样式
p1_filtered_with_marginal  # 显示图形


