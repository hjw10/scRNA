#rm(list = ls())
#V5_path = "/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/seurat5/"
#.libPaths(V5_path)
#.libPaths()
library(stringr)
library(Seurat)
library(Augur)
library(dplyr)
library(patchwork)
library(viridis)
library(BiocParallel)
library(ggplot2)

load("/public/workspace/stu21230110/scRNA/06_scRNA_harmony_singler202400810.Rdata")
outdir = '/public/workspace/stu21230110/design/08Augur/results'

# check
scRNA = scRNA_harmony_singler
table(scRNA$orig.ident)
table(scRNA$Annotation)


# 此函数的主要目的是评估并计算每种细胞类型在特定条件（如实验组别或时间点）下的表现，通常通过计算区域下曲线（AUC）来量化。它是单一条件下的细胞类型表现的度量，可以帮助识别在特定实验设置下哪些细胞类型表现最为显著或重要。
sc_augur <- calculate_auc(scRNA,cell_type_col="Annotation",label_col="group",n_threads=4)

saveRDS(sc_augur,"augur.rds")
#sc_augur <- readRDS("augur.rds")
head(sc_augur$AUC, 5)

# A tibble: 11 × 2
   cell_type                  auc
   <chr>                    <dbl>
 1 Epithelial_cells         0.766
 2 B_cells                  0.759
 3 Tissue_stem_cells        0.665
 4 T_cells                  0.652
 5 Smooth_muscle_cells      0.634
 6 Monocyte                 0.615
 7 NK_cells                 0.604
 8 Stromal_cells            0.587
 9 Endothelial_cells        0.582
10 Myeloid_progenitor_cells 0.581
11 Macrophage               0.564

# 明确哪个细胞群体的RNA水平变化最大
# plot_lollipop函数是基于ggplot框架的,所以很容易就可以美化图片
# 还可以提取原始函数进行修改
setwd("results")
pdf("sc_augur.pdf")
plot_lollipop(sc_augur) +
    geom_segment(aes(xend = cell_type,yend = 0.5), size = 1) +
    geom_point(size = 3,aes(color = cell_type)) + # 增加点的大小和颜色映射
    scale_color_manual(values = c(
        "B_cells"="#F39B7FFF",
        "Endothelial_cells"="#E377C2FF",
        "Epithelial_cells"="#3C5488FF",
        "Macrophage"="#7E6148FF",
        "Monocyte"="#00A087FF",
        "Myeloid_progenitor_cells"="#8494FF",
        "NK_cells"="#4DBBD5FF",
        "Smooth_muscle_cells"="#9467BDFF",
        "Stromal_cells"="#E64B35FF",
        "T_cells"="#91D1C2FF",
        "Tissue_stem_cells"="#B09C85FF"
))
dev.off()

# 结合umap
# rank模式
pdf('rank_umap.pdf')
plot_umap(sc_augur,
            scRNA,
            mode="rank",
            reduction="umap",
            palette="cividis",
#            augur_mode="default",
            cell_type_col="Annotation")
dev.off()

# default模式
pdf('default_umap.pdf')
plot_umap(sc_augur,
            scRNA,
            mode="default",
            reduction="umap",
            palette="viridis", #"viridis", "plasma", "magma", "inferno" "YlGnBu"
#            augur_mode="default",
            cell_type_col="Annotation")
dev.off()

# 这个函数使用来自 calculate_auc 的输出进行更深入的分析，比较两个不同条件下的AUC分数，以识别统计上显著的差异。通过置换测试（permute），该函数评估两个条件下细胞类型优先级的变化是否显著，超过了随机变化的可能性

# 拆分数据-笔者假设batch是不同的药物
scRNA.list <- SplitObject(scRNA, split.by = "batch")

# 初始化结果列表
augur_results <- list()
permuted_results <- list()

# 计算每个批次的AUC
for (i in seq_along(scRNA.list)) {
  # 从列表中获取当前的Seurat对象
  current_scRNA <- scRNA.list[[i]]

  # 正常模式计算AUC
  normal_augur <- calculate_auc(current_scRNA,
                                cell_type_col = "celltype",
                                label_col = "group",
                                augur_mode = "default",  # 使用默认模式
                                n_threads = 8)
  augur_results[[i]] <- normal_augur

  # 置换模式计算AUC
  permuted_augur <- calculate_auc(current_scRNA,
                                  cell_type_col = "celltype",
                                  label_col = "group",
                                  augur_mode = "permute",  # 使用置换模式
                                  n_threads = 8)
  permuted_results[[i]] <- permuted_augur
}

# 保存结果
qsave(augur_results, "augur_normal_list.qs")
qsave(permuted_results, "augur_permuted_list.qs")

# 如果需要比较第一个和第二个批次的差异
res <- calculate_differential_prioritization(augur_results[[1]],
                                             augur_results[[2]],
                                             permuted_results[[1]],
                                             permuted_results[[2]],
                                             n_subsamples = 50,
                                             n_permutations = 1000)
res
# # A tibble: 5 × 9
#   cell_type    auc.x auc.y delta_auc     b     m     z    pval    padj
#   <fct>        <dbl> <dbl>     <dbl> <int> <int> <dbl>   <dbl>   <dbl>
# 1 Adipocytes   0.753 0.671   -0.0824   998  1000 -3.13 0.00599 0.00599
# 2 CD8+ T-cells 0.734 0.639   -0.0954  1000  1000 -4.45 0.00200 0.00500
# 3 Fibroblasts  0.711 0.790    0.0786     0  1000  4.06 0.00200 0.00500
# 4 Monocytes    0.619 0.679    0.0603     1  1000  2.93 0.00400 0.00500
# 5 CD4+ T-cells 0.588 0.645    0.0573     1  1000  3.20 0.00400 0.00500

plot_scatterplot(augur_results[[1]],augur_results[[2]],top_n = 6)+
    geom_point(size = 3, aes(color = cell_type))   # 增加点的大小和颜色映射

