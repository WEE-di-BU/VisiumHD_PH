# 加载必要的R包（保留原所有包，新增所需包）
library(Seurat)
library(dplyr)
library(AUCell)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(data.table)
library(harmony)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(ggrepel)
library(rjags)
library(infercnv)
library(CellChat)
library(ComplexHeatmap)
library(scRNAtoolVis)
library(RColorBrewer)
library(MetBrewer)
library(tidyverse)
library(SingleR)       # 新增
library(celldex)       # 新增
library(sctransform)   # 新增
library(clustree)      # 新增

# 定义颜色方案（完全保留原代码）
mycolor <- c(
  "#84B8D7", "#4391A9", "#8DC594", "#72BF5A", "#B49D99", "#F6955A", "#D9CE99", 
  "#FDAE53", "#408DBF", "#D8AC60", "#969D62", "#92CF72", "#659E47", "#F47A79", 
  "#E93A3B", "#B294C7", "#815AA8", "#EBD57C", "#C48244", "#B15928", "#E94330", 
  "#E4986B", "#F06C45", "#FE9E37", "#984ea3", "#377eb8", "#4daf4a", "#e41a1c", 
  "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#3d7bb0", 
  "#a2c8db", "#5aa554", "#afd195", "#d83f36", "#e99997", "#f2bc7c", "#e68740", 
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3", "#476D87", 
  "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A", "#8C549C", "#585658", 
  "#9FA3A8", "#E0D4CA", "#a6cee3", "#fdbf6f", "#b2df8a", "#ffff99", "#33a02c", 
  "#1f78b4", "#ff7f00", "#cab2d6", "#6a3d9a", "#fb9a99", "#e31a1c", "#7db0e9", 
  "#3d7bb0", "#57a65b", "#7fcd87", "#b9364e", "#c56585", "#7e4878", "#513c65"
)
colors22 <- c(
  "C0-CD4_T cell" = "#84B8D7", 
  "C1-Follicular B cell" = "#4391A9", 
  "C2-NK" = "#8DC594", 
  "C3-Stem cell" = "#72BF5A", 
  "C4-Plasma cell" = "#B49D99", 
  "C5-CD14_Myeloid cell" = "#F6955A", 
  "C6-NK" = "#D9CE99", 
  "C7-Ciliated cell" = "#FDAE53", 
  "C8-Transit Amplifying Cell" = "#408DBF", 
  "C9-TAM" = "#D8AC60", 
  "C10-Goblet cell" = "#969D62", 
  "C11-Germinal center B cell" = "#92CF72", 
  "C12-Progenitor cell" = "#659E47", 
  "C13-Fibroblast" = "#F47A79", 
  "C14-Transit amplifying cell" = "#E93A3B", 
  "C15-Brush cell" = "#B294C7", 
  "C16-Paneth cell" = "#815AA8", 
  "C17-Progenitor cell" = "#EBD57C", 
  "C18-Enteroendocrine cell" = "#C48244", 
  "C19-Mast cell" = "#B15928", 
  "C20-Endothelial cell" = "#E94330", 
  "C21-Stem cell" = "#E4986B"
)

# 数据加载和预处理
setwd('D:/R_Projects/CRC_paper-main/')
data_dirs <- c(
  "data/Patient0_carcinoma",
  "data/Patient1_adenoma",
  "data/Patient1_carcinoma"
)

CRC <- CreateSeuratObject(
  counts = Read10X(data.dir = data_dirs),
  min.cells = 5,
  min.features = 300
)

sample_map <- c(
  "1" = "P0_carcinoma",
  "2" = "P1_adenoma",
  "3" = "P1_carcinoma"
)
CRC$sample <- factor(
  CRC$orig.ident,
  levels = c(1, 2, 3),
  labels = c("P0_carcinoma", "P1_adenoma", "P1_carcinoma")
)

# 数据质控和过滤
CRC[["percent.mt"]] <- PercentageFeatureSet(CRC, pattern = "^MT-")
CRC = subset(CRC, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
CRC$sample_type <- ifelse(CRC$sample %in% c('P1_normal', 'P2_normal', 'P3_normal'), 'normal',
                          ifelse(CRC$sample %in% c('P1_adenoma', 'P2_adenoma', 'P3_adenoma1', 'P3_adenoma2'), 'adenoma',
                                 ifelse(CRC$sample %in% c('P0_carcinoma', 'P1_carcinoma', 'P2_carcinoma', 'P3_carcinoma'), 'carcinoma',
                                        ifelse(CRC$sample %in% c('P2_paracancer'), 'paracancer', NA))))

######################################################
### 技术改进1：使用SCTransform替代常规归一化 ###
######################################################

# 原代码：
# CRC<- NormalizeData(CRC)
# CRC = FindVariableFeatures(CRC,selection.method = "vst", nfeatures = 2000, verbose = F)
# CRC <- ScaleData(CRC)

# 改进代码：
CRC <- SCTransform(CRC, 
                   vars.to.regress = c("percent.mt", "nCount_RNA"),
                   verbose = FALSE)

# 创建基础图，报告中图1
plot <- VariableFeaturePlot(CRC)
print(plot)

# 获取top10可变基因
top10_var_genes <- head(VariableFeatures(CRC), 10)

# 添加标签
plot + LabelPoints(plot = plot, 
                   points = top10_var_genes, 
                   repel = TRUE) +
  labs(title = "Top 10 variable genes after SCTransform")

######################################################
### 技术改进2：使用Leiden算法进行聚类 ###
######################################################

# 保留原Harmony批次校正
CRC <- RunPCA(CRC, verbose = FALSE)
CRC <- RunHarmony(CRC, "sample", plot_convergence = TRUE)


# 改进代码：
CRC <- RunUMAP(CRC, reduction = "harmony", dims = 1:20)
CRC <- FindNeighbors(CRC, reduction = "pca", dims = 1:20, 
                     graph.name = "sct_leiden")

# 测试不同分辨率
resolutions <- c(0.4, 0.6, 0.8, 1.0)
for (res in resolutions) {
  CRC <- FindClusters(CRC, 
                      resolution = res, 
                      algorithm = 4,  # 4代表Leiden算法
                      graph.name = "sct_leiden")
}

# 验证聚类列已添加
grep("sct_leiden_res", colnames(CRC@meta.data), value = TRUE)

# 查看聚类结果分布
table(CRC$sct_leiden_res.1)

# 可视化特定分辨率的聚类
DimPlot(CRC, group.by = "sct_leiden_res.1", label = TRUE) + NoLegend()

# 可视化不同分辨率下的聚类，论文中第二张图
clustree_plot <- clustree(CRC, prefix = "sct_leiden_res.") +
  ggtitle("Clustering at different resolutions (Leiden algorithm)") + theme(legend.position = "bottom")
print(clustree_plot)

# 选择最佳分辨率（示例使用0.6）
selected_resolution <- 0.6
CRC$seurat_clusters <- CRC[[paste0("sct_leiden_res.", selected_resolution)]]
Idents(CRC) <- CRC$seurat_clusters

# 可视化聚类结果，论文图三
DimPlot(CRC, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
  ggtitle(paste("Leiden clustering at resolution", selected_resolution))

######################################################
### 技术改进3：整合SingleR自动化注释 ###
######################################################

# 原代码使用COSG和手动注释：
# library(COSG)
# COSG_markers <- cosg(CRC,groups='all',assay='RNA',slot='data',mu=1,n_genes_user=20)
# CRC_markergene = as.data.frame(COSG_markers$names)
# 然后手动定义细胞类型...

# 改进代码：
# 下载参考数据集（首次运行需要下载）
hpca.ref <- HumanPrimaryCellAtlasData()
blueprint.ref <- BlueprintEncodeData()

# 准备SingleR输入
sce <- as.SingleCellExperiment(CRC)

# 使用HPCA参考进行注释
pred.hpca <- SingleR(test = sce, 
                     ref = hpca.ref, 
                     labels = hpca.ref$label.main,
                     assay.type.test = "logcounts")

# 使用Blueprint参考进行注释
pred.blueprint <- SingleR(test = sce, 
                          ref = blueprint.ref, 
                          labels = blueprint.ref$label.main,
                          assay.type.test = "logcounts")

# 将注释结果添加到Seurat对象
CRC$hpca_labels <- pred.hpca$labels
CRC$blueprint_labels <- pred.blueprint$labels

# 比较两种参考的注释结果
annotation_compare <- table(HPCA = CRC$hpca_labels, 
                            Blueprint = CRC$blueprint_labels)
print(annotation_compare)

# 创建一致性注释
CRC$consensus_labels <- ifelse(CRC$hpca_labels == CRC$blueprint_labels,
                               CRC$hpca_labels,
                               paste("Discrepancy:", CRC$hpca_labels, "vs", CRC$blueprint_labels))

# 可视化自动化注释结果，报告图5
p1 <- DimPlot(CRC, group.by = "hpca_labels", label = TRUE, repel = TRUE) + 
  ggtitle("HPCA annotation") + NoLegend()
p2 <- DimPlot(CRC, group.by = "blueprint_labels", label = TRUE, repel = TRUE) + 
  ggtitle("Blueprint annotation") + NoLegend()
p1 + p2

######################################################
### 结合自动注释和手动验证 ###
######################################################

# 保留原COSG标记基因识别（用于验证）
library(COSG)
COSG_markers <- cosg(CRC, groups='all', assay='SCT', slot='data', mu=1, n_genes_user=20)
CRC_markergene <- as.data.frame(COSG_markers$names)

# 可视化标记基因（使用SCT assay），报告图6
DotPlot(CRC, 
        features = unique(unlist(head(CRC_markergene, 5))), 
        assay = "SCT",
        group.by = "seurat_clusters") + 
  RotatedAxis() +
  scale_color_gradientn(colors = rev(brewer.pal(11, "RdBu")))

# 基于自动注释和标记基因进行最终注释
# 这里保留原手动注释框架，但结合自动注释结果
CRC_seurat_cluster2 <- c(
  '0' = 'CD4_T_cell',
  '1' = 'Follicular_B_cell',
  '2' = 'NK',
  '3' = 'Stem_cell',
  '4' = 'Plasma_cell',
  '5' = 'CD14_Myeloid_cell',
  '6' = 'NK_activated',
  '7' = 'Ciliated_cell',
  '8' = 'Transit_Amplifying_Cell',
  '9' = 'TAM',
  '10' = 'Goblet_cell',
  '11' = 'Germinal_center_B_cell',
  '12' = 'Progenitor_cell',
  '13' = 'Fibroblast',
  '14' = 'Transit_amplifying_cell',
  '15' = 'Brush_cell',
  '16' = 'Paneth_cell',
  '17' = 'Progenitor_cell_2',
  '18' = 'Enteroendocrine_cell',
  '19' = 'Mast_cell',
  '20' = 'Endothelial_cell',
  '21' = 'Stem_cell_2'
)

# 添加最终注释
CRC[['seurat_clusters2']] <- unname(CRC_seurat_cluster2[as.character(CRC$seurat_clusters)])
CRC$new_labels <- paste0('C', CRC$seurat_clusters, '-', CRC$seurat_clusters2)

# 排序因子水平
a <- names(table(CRC$new_labels))
b <- gsub("-.*$", "", a)
c <- a[order(as.numeric(gsub("C", "", b)))]
CRC$celltype <- factor(CRC$new_labels, levels = c)

# 最终可视化，报告中图四
final_plot <- DimPlot(CRC, 
                      group.by = "celltype", 
                      label = TRUE, 
                      pt.size = 0.5,
                      label.size = 3,
                      repel = TRUE,
                      cols = mycolor) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  ggtitle("Final cell type annotation with Leiden clustering")

print(final_plot)


###############################################
### 结果保存与可视化 ###
###############################################
qsave(CRC, file = "CRC_annotated_improved.qs")

# 生成注释报告
annotation_report <- data.frame(
  Cluster = CRC$seurat_clusters,
  HPCA = CRC$hpca_labels,
  Blueprint = CRC$blueprint_labels,
  Final_Annotation = CRC$celltype
)

write.csv(annotation_report, file = "cell_annotation_report.csv", row.names = FALSE)

# 保存重要可视化结果
ggsave("D:/R_Projects/CRC_paper-main/figures/final_umap_annotated.png", plot = final_plot, width = 10, height = 8, dpi = 300)
ggsave("D:/R_Projects/CRC_paper-main/figures/clustree_resolution.png", plot = clustree_plot, width = 12, height = 8, dpi = 300)
