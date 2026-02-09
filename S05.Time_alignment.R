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

#########################################
ht <- Heatmap(as.matrix(groupMatRNA),
              rect_gp = gpar(col = "black", lwd = 0.5),
             #col = ArchR::paletteContinuous("solarExtra", n = 100),
             cluster_rows = TRUE,
             cluster_columns = FALSE,
             name = "RNA")
#############################################
save(groupMatRNA, abc, file = 'Select_gene.Rdata')
#############################################
library(TCGAbiolinks)
abc <- cl$cluster
abc <- cbind(names(abc),abc)
###############
clustering <- as.data.frame(abc)
colnames(clustering)[1] <- 'gene'
Cluster1 <- clustering[which(clustering$abc == '1'),'gene']
Cluster2 <- clustering[which(clustering$abc == '2'),'gene']
Cluster3 <- clustering[which(clustering$abc == '3'),'gene']
########################
ansEA <- TCGAanalyze_EAcomplete(
    TFname = "Phase3",
    RegulonList = Cluster3
)  

####### Plot the pathway
TCGAvisualize_EAbarplot(
    tf = rownames(ansEA$ResBP),
    GOBPTab = ansEA$ResBP,
    GOCCTab = ansEA$ResCC,
    GOMFTab = ansEA$ResMF,
    PathTab = ansEA$ResPat,
    nRGTab = Cluster3,
    nBar = 10,
    text.size = 2,
    fig.width = 25,
    fig.height = 15
)
##############################################
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db) ##加载小鼠
library(org.Hs.eg.db) ##加载人类
######################################
##########################
GO_A = bitr(Cluster1, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

YO_bp = enrichGO(gene=GO_A$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 

YO_UP_GO_A_result <- as.data.frame(YO_bp@result)

color_1 <- c("blue","red")
ggplot(YO_UP_GO_A_result[1:10,])+
  geom_point(aes(RichFactor,Description,
                 color = pvalue,
                 size = Count))+
  labs(x="GeneRatio",y="GO description") + 
  labs(title="")+
  scale_color_gradient(low = color_1[1],high=color_1[2],name="Pvalue")+
  theme_bw() + ylab('') + ggtitle('Cluster1')
##########################################
GO_A = bitr(Cluster2, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

YO_bp = enrichGO(gene=GO_A$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 

YO_UP_GO_A_result <- as.data.frame(YO_bp@result)

color_1 <- c("blue","red")
ggplot(YO_UP_GO_A_result[c(1:2,4:11),])+
  geom_point(aes(RichFactor,Description,
                 color = pvalue,
                 size = Count))+
  labs(x="GeneRatio",y="GO description") + 
  labs(title="")+
  scale_color_gradient(low = color_1[1],high=color_1[2],name="Pvalue")+
  theme_bw() + ylab('') + ggtitle('Cluster2')
##########################################
GO_A = bitr(Cluster3, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

YO_bp = enrichGO(gene=GO_A$ENTREZID, OrgDb = org.Hs.eg.db,ont= "BP",pAdjustMethod = "BH",pvalueCutoff= 0.05,qvalueCutoff= 0.2,readable= TRUE) 

YO_UP_GO_A_result <- as.data.frame(YO_bp@result)

color_1 <- c("blue","red")
ggplot(YO_UP_GO_A_result[1:10,])+
  geom_point(aes(RichFactor,Description,
                 color = pvalue,
                 size = Count))+
  labs(x="GeneRatio",y="GO description") + 
  labs(title="")+
  scale_color_gradient(low = color_1[1],high=color_1[2],name="Pvalue")+
  theme_bw() + ylab('') + ggtitle('Cluster3')
#########################################  Trajectory
gene_list <- list(
    "Cluster1" = Cluster1[1:10],
    "Cluster2" = Cluster2[35:45],
    "Cluster3" = Cluster3[1:10]
)

abc <- groupMatRNA
plot_df <- reshape2::melt(abc)
plot_df$pseudotime <- rep(seq(1:100),2277)
dir_for_result <- 'G:/YouRui/CD8Tcell/Try'

pbmc <- coembed
all.genes <- rownames(pbmc)                                       
pbmc <- ScaleData(pbmc, features = all.genes)

for (i_name in names(gene_list)) {
    genes <- gene_list[[i_name]]
    plot_df <- pbmc@assays$RNA@scale.data[genes, ]
    colnames(plot_df) <- colnames(pbmc)
    plot_df <- reshape2::melt(plot_df)
    plot_df$pseudotime <- pbmc$Trajectory[match(plot_df$Var2, colnames(pbmc))]
    plot_df <- plot_df[!is.na(plot_df$value), ]

    ggplot(data = plot_df, mapping = aes(x = pseudotime, y = value, color = Var1, fill = Var1)) +
        # geom_point() +
        ggsci::scale_color_npg(name = "") +
        ggsci::scale_fill_npg(name = "") +
        geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), alpha = 0.2, size = 0.5) +
        cowplot::theme_cowplot() +
        theme(
            axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 7),
            text = element_text(size = 7, family = "ArialMT"),
            plot.margin = unit(c(1, 1, 1, 1), "char"),
            plot.title = element_text(hjust = 0.5, size = 8, family = "ArialMT", face = "plain"),
            axis.line = element_line(linetype = 1, color = "black", size = 0.3),
            axis.ticks = element_line(linetype = 1, color = "black", size = 0.3),
            legend.position = "none"
        ) +
        ylab("Scaled expression") +
        xlab("Pseudo-time") +
        ggtitle(i_name)
    ggsave(file.path(dir_for_result, paste0("S6E.", i_name, "_gene_trends.pdf")),
        width = 1.8,
        height = 2
    )
}      
