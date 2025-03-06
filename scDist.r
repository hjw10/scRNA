library(stringr)
library(Seurat)
library(scDist)
library(dplyr)
library(patchwork)
library(viridis)
library(ggplot2)
library(BiocParallel)

load("/public/workspace/stu21230110/scRNA/06_scRNA_harmony_singler202400810.Rdata")
options(future.globals.maxSize = 100000000000)

# check
scRNA = scRNA_harmony_singler
table(scRNA$orig.ident)
table(scRNA$singleRnew)

scRNA <- SCTransform(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)
scRNA <- RunPCA(scRNA, verbose = FALSE)
scRNA <- RunUMAP(scRNA, dims = 1:30, verbose = FALSE)
scRNA <- FindNeighbors(scRNA, dims = 1:30, verbose = FALSE)
scRNA <- FindClusters(scRNA, verbose = FALSE, resolution = 0.05)
pdf("scRNA.SCT.DimPlot.pdf")
DimPlot(scRNA,label = TRUE)
dev.off()

dat <- list(Y=scRNA@assays$SCT@scale.data %>% as.data.frame(),
            meta.data=scRNA@meta.data %>% as.data.frame())
dim(dat$Y)
head(dat$meta.data)
out <- scDist(dat$Y,dat$meta.data,
              d = 20,
              fixed.effects="group",
              random.effects="orig.ident",
              clusters="singleRnew")

saveRDS(out,"scDist.rds")
saveRDS(scRNA,"scRNA.SCT.rds")

out$results

#原始矩阵
                            Dist. 95% CI (low) 95% CI (upper)       p.val
B_cells                   8.19566      0.00000       34.59698 0.706572934
Endothelial_cells         0.00000      0.00000        0.00000 0.869981300
Epithelial_cells         19.43759     14.96266       31.89521 0.009239908
Macrophage                0.00000      0.00000        0.00000 0.987170128
Monocyte                  0.00000      0.00000        0.00000 0.924330757
Myeloid_progenitor_cells 14.54827      0.00000       26.61902 0.628533715
NK_cells                  0.00000      0.00000       10.22942 0.455495445
Smooth_muscle_cells       0.00000      0.00000        0.00000 0.753642464
Stromal_cells             0.00000      0.00000        0.00000 0.766452335
T_cells                   0.00000      0.00000       35.05571 0.504964950
Tissue_stem_cells        50.65466      0.00000       84.73033 0.348216518

#SCT
                             Dist. 95% CI (low) 95% CI (upper)     p.val
B_cells                   6.370027     3.146784      10.163078 0.3819162
Endothelial_cells        12.400859     0.000000      21.772109 0.5417246
Epithelial_cells          0.000000     0.000000       0.000000 0.7902921
Macrophage                0.000000     0.000000       0.000000 0.6769232
Monocyte                  1.894509     0.000000       7.011164 0.5477145
Myeloid_progenitor_cells  4.542387     2.526224       6.357726 0.3296867
NK_cells                  0.000000     0.000000       0.000000 0.6404836
Smooth_muscle_cells       0.000000     0.000000       0.000000 0.8912511
Stromal_cells             0.000000     0.000000       0.000000 0.4737253
T_cells                   3.940336     0.000000       6.670926 0.5814642
Tissue_stem_cells         0.000000     0.000000       0.000000 0.6277737



setwd("results")
pdf("SCT.displot.pdf")
DistPlot(out)
dev.off()

pdf("SCT.epi_disGenes.pdf")
distGenes(out,cluster="Epithelial_cells")
dev.off()

pdf("SCT.epi_plotBetas.pdf")
plotBetas(out,cluster="Epithelial_cells")
dev.off()

pdf("SCT.FDRDisPlot.pdf")
FDRDistPlot(out)
dev.off()
