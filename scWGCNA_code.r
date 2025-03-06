##########################################################15hdWGCNA-20240622
setwd("/public/workspace/stu21230110/design/")
# 推荐本地安装
#devtools::install_local('hdWGNCA.zip')
#BiocManager::install('harmony',update=F,ask=F)
library(hdWGCNA)
#加载单细胞分析包
library(Seurat)
#加载作图包
library(tidyverse)
library(cowplot)
library(patchwork)
#加载共表达网络分析包
library(WGCNA)
load("01RDada/06_scRNA_harmony_singler202400810.Rdata")


choose_cell <- c("Macrophage","Stromal_cells","Tissue_stem_cells") ###巨噬细胞 基质细胞 组织干细胞 
# 读取单细胞数据集，读取自己的
scRNA_harmony<-scRNA_harmony_singler

if (!dir.exists("06hdWGCNA")){
dir.create("06hdWGCNA")
}

setwd("./06hdWGCNA")#切换工作目录

if (!dir.exists("Macrophage")){
dir.create("Macrophage")
}
if (!dir.exists("Stromal_cells")){
dir.create("Stromal_cells")
}
if (!dir.exists("Tissue_stem_cells")){
dir.create("Tissue_stem_cells")



# 提细胞亚群,重新降维聚类

#巨噬细胞 基质细胞 干细胞
#setwd("/public/workspace/stu21230110/design/06hdWGCNA/Macrophage")
#setwd("/public/workspace/stu21230110/design/06hdWGCNA/Stromal_cells")
setwd("/public/workspace/stu21230110/design/06hdWGCNA/Tissue_stem_cells")
#scRNA<-subset(scRNA_harmony,singleRnew=="Macrophage")
#scRNA<-subset(scRNA_harmony,singleRnew=="Stromal_cells")
scRNA<-subset(scRNA_harmony,singleRnew=="Tissue_stem_cells")

####需要对scRNA的metadata
scRNA@meta.data$orig.ident2<-scRNA@meta.data$orig.ident
scRNA@meta.data$orig.ident<-as.character(scRNA@meta.data$group)
scRNA@meta.data$celltype<-scRNA@meta.data$singleRnew

scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)
scRNA<- RunPCA(scRNA, features = VariableFeatures(scRNA))
DimPlot(scRNA, reduction = "pca", group.by = "orig.ident")
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:20)
scRNA <- FindClusters(scRNA,resolution=0.3)
scRNA <- RunUMAP(scRNA, dims = 1:20)

#


# 瞄一眼，小亚群
pdf("01dimplot.pdf", height=10, width=18)
DimPlot(scRNA, label=TRUE) 
dev.off()
#过滤出至少在5%的细胞中表达的基因
scRNA <- SetupForWGCNA(
  scRNA,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Bio_com" # the name of the hdWGCNA experiment
)


#构建metacells!!这一步非常重要，WGCNA对数据的稀疏性非常敏感，与普通转录组的WGCNA分析相比
# 单细胞的稀疏矩阵的解决方法是WGCNA流程的问题核心

# construct metacells  in each group
scRNA<- MetacellsByGroups(
  seurat_obj = scRNA,k=20,
  max_shared = 10,
  # group.by一般关注的是组织类型和细胞类型!这边组织类型是orig.ident，CT正常，PR疾病
  group.by = c("celltype",'orig.ident'), # 也可以选择别的groupby
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
scRNA <- NormalizeMetacells(scRNA)
metacell_obj <- GetMetacellObject(scRNA)


#转置表达矩阵
# 安全起见，另起一个对象，以巨噬细胞细胞为例
seurat_obj  <- SetDatExpr(
  scRNA,
#  group_name = "Macrophage", 
#  group_name = "Stromal_cells",
  group_name = "Tissue_stem_cells",# 选择感兴趣恶的细胞类群！
  group.by='celltype' # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
)


#选择softpower
seurat_obj <- TestSoftPowers(
  seurat_obj,
  setDatExpr = FALSE, # 这边不用转置了，前面已转置
)


# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
pdf("02wrap_plots.pdf", height=10, width=18)
wrap_plots(plot_list, ncol=2)
dev.off()
# assemble with patchwork




#查看powerTable
power_table <- GetPowerTable(seurat_obj)
head(power_table)


#构建共表达网络
softpower=7  # 根据自己的改，选5也没问题
# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=softpower,
  group.by='celltype', 
  #group_name='Macrophage',
  #group_name='Stromal_cells',
  group_name='Tissue_stem_cells',
  setDatExpr = F)


#可视化WGCNA网络

pdf("03Dendrogram_Plot.pdf", height=10, width=10)
PlotDendrogram(seurat_obj, main='hdWGCNA Dendrogram')
dev.off()

#(可选)获取TOM矩阵，可用来进行其他高级分析
TOM <- GetTOM(seurat_obj)


#计算模块协调特征
#记得scale一下 or else harmony throws an error:
seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = GetWGCNAGenes(seurat_obj),
  
)
# 计算ME，根据组织类型分组
# harmony必须biocManager安装，不可以用github安装！！！
library(harmony)

seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident" #harmony对象
)


seurat_obj <- ModuleConnectivity(seurat_obj)
# plot genes ranked by kME for each module
#可视化每个模块中，按照kME打分的基因

pdf("04KME_Plot.pdf", height=10, width=30)
PlotKMEs(seurat_obj, ncol=3)
dev.off()




# 获取hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 25)
head(hub_df)
pdf("05KME_Plot.pdf", height=10, width=30)
PlotKMEs(seurat_obj)
dev.off()
#记得保存上面hdWGNCA关键分析过程！！！
#saveRDS(seurat_obj, file='hdWGCNA_object_Macrophage.rds')
#saveRDS(seurat_obj, file='hdWGCNA_object_Stromal_cells.rds')
saveRDS(seurat_obj, file='hdWGCNA_object_Tissue_stem_cells.rds')

intersect(hub_df$gene_name,macrophage_marker[which(macrophage_marker$p_val_adj < 1),1])
 [1] "NPC2"     "LYZ"      "B2M"      "TYMP"     "C1QC"     "TYROBP"
 [7] "ACTB"     "AIF1"     "HLA-DRB1" "TMSB4X"   "IFI30"    "CD74"
[13] "HLA-DQA1" "HLA-DRA"  "FTL"      "RPS15"    "RPL15"    "RPS19"
[19] "RPS10"    "FCGBP"    "LAPTM4A"  "WFDC2"    "DPP6"     "IGLC2"
[25] "PGR"      "PALLD"    "IGFBP7"   "FBLN1"    "SPARC"    "ECM1"
[31] "RAMP1"    "SDK1"     "MDK"      "IGFBP5"   "KCNIP4"   "MGP"
[37] "COL1A2"   "COL3A1"   "COL1A1"   "SFRP4"    "MMP11"



####------------一些可视化-----------------------
## 模块间的相关性
library(igraph)
library(qgraph)
# 载入保存的
#seurat_obj=readRDS('e:/writing/benke/hdWGCNA/hdWGCNA_object.rds')
# 画模块间相关性图
pdf("06ModuleCorrelogram_Plot.pdf", height=10, width=10)
ModuleCorrelogram(seurat_obj, sig.level = 0.001, pch.cex=2)
dev.off()
# 由于识别到了每个模块的hub基因，可以去计算hub基因打分
# compute gene scoring for the top 25 hub genes by kME for each module
# (方法一)with Seurat method------不用此方法
#seurat_obj <- ModuleExprScore(
#  seurat_obj,
#  n_genes = 25,
#  method='Seurat'
#)
# compute gene scoring for the top 25 hub genes by kME for each module
# (方法二)with UCell method #推荐这种方法
# 由于Ucell刚刚更新，所以4.1.x的同学请用本地安装,依赖包自行安装
#devtools::install_local("d:/R/UCell-1.3.zip")
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# featureplot
# 瞄一眼
pdf("07featureplot.pdf", height=10, width=10)
DimPlot(scRNA, label=TRUE,split.by = 'orig.ident') 
dev.off()


plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE ,# order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf("08plot_list.pdf", height=10, width=10)
wrap_plots(plot_list)
dev.off()

### dotplot
# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods)

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

# plot output

pdf("09output_plot.pdf", height=10, width=10)
p
dev.off()
######公众号缺少了一部分
seurat_obj$orig.ident<-as.factor(seurat_obj$orig.ident)
cur_traits <- c('orig.ident','nCount_RNA')

seurat_obj <- ModuleTraitCorrelation(
  seurat_obj,
  traits = cur_traits,
  features = "hMEs",
  cor_method = "pearson",
  group.by='orig.ident'
)
mt_cor <- GetModuleTraitCorrelation(seurat_obj)
names(mt_cor)
pdf("10tModuleTraitCorrelation_plot.pdf", height=10, width=10)
PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
   label_symbol = 'stars',
    text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = '#fc9272',
   mid_color = '#ffffbf',
   low_color = '#9ecae1',
   plot_max = 0.2,
  combine=T 
  )
dev.off()


######################################第三节
#单细胞也能实现WGCNA分析！快来学习hdWGCNA！（三）hdWGCNA配套富集/网络可视化大放送！（三）

#加载seurat数据和包
# single-cell analysis package
##setwd('e:/writing/benke/hdWGCNA/')
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
#install.packages('enrichR')
library(enrichR)

#BiocManager::install('GeneOverlap',update = F,ask = F)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
#set.seed(12345)

# load the Zhou et al snRNA-seq dataset
#seurat_obj <- readRDS('hdWGCNA_object.rds')

#GO富集分析
# enrichr databases to test
dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021')

# 富集分析，会逐个模块分析
seurat_obj <- RunEnrichr(
  seurat_obj,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(seurat_obj)

# make GO term plots作图，在文件夹下生成！
pdf("11EnrichrBarPlot.pdf", height=10, width=10)


EnrichrBarPlot(
  seurat_obj,
  outdir = "enrichr_plots", # name of output directory
  n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
  plot_size = c(5,7), # width, height of the output .pdfs
  logscale=TRUE # do you want to show the enrichment as a log scale?
)
dev.off()
#气泡图
# GO_Biological_Process_2021
pdf("12EnrichrBarPlot.pdf", height=10, width=10)
EnrichrDotPlot(
  seurat_obj,
#  mods = c("turquoise","black"), # use all modules (this is the default behavior)
  database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)
dev.off()

#气泡图
# GO_Cellular_Component_2021
pdf("13EnrichrBarPlot.pdf", height=10, width=10)
EnrichrDotPlot(
  seurat_obj,
#  mods = c("turquoise","black",'blue'), # use all modules (this is the default behavior)
  database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)
dev.off()
#气泡图
# GO_Biological_Process_2021

pdf("14EnrichrBarPlot.pdf", height=10, width=10)
EnrichrDotPlot(
  seurat_obj,
#  mods = c("turquoise","black",'red'), # use all modules (this is the default behavior)
  database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
  n_terms=6
  # number of terms for each module
)

dev.off()
#差异基因重叠分析
## 这个分析帮助我们看到，哪些模块可能是相似的
# compute cell-type marker genes with Seurat:
# 常规方法计算差异基因/特征基因
Idents(seurat_obj) <- seurat_obj$seurat_clusters
markers <- Seurat::FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  logfc.threshold=1
)

# compute marker gene overlaps
overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

#条形图
# overlap barplot, produces a plot for each cell type
plot_list <- OverlapBarPlot(overlap_df)

# stitch plots with patchwork


pdf("15OverlapBarPlot.pdf", height=10, width=40)
wrap_plots(plot_list, ncol=4)
dev.off()
#气泡图
# plot odds ratio of the overlap as a dot plot
pdf("16OverlapBarPlot.pdf", height=10, width=10)
OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')
dev.off()
#----------------------------
#网络可视化
# network analysis & visualization package:
# network analysis & visualization package:
library(igraph)

#可视化每个模块的网络图
pdf("17ModuleNetworkPlot.pdf", height=10, width=10)
ModuleNetworkPlot(seurat_obj)
dev.off()
#组合网络图，在文件夹下生成
# hubgene network
pdf("18ModuleNetworkPlot.pdf", height=10, width=10)
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = "all"
)

dev.off()
#UMAP可视化
## 利用hub基因，重新UMAP，如此可以获得分群明显的图
g <- HubGeneNetworkPlot(seurat_obj,  return_graph=TRUE)

seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)


# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf("19ggplot.pdf", height=10, width=10)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()
dev.off()
pdf("20ggplot.pdf", height=10, width=10)
ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE

)
dev.off()
#监督UMAP
g <- ModuleUMAPPlot(seurat_obj,  return_graph=TRUE)
# run supervised UMAP:
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 10,
  n_neighbors=15,
  min_dist=0.1,
  supervised=TRUE,
  target_weight=0.5
)

# get the hub gene UMAP table from the seurat object
umap_df <- GetModuleUMAP(seurat_obj)

# plot with ggplot
pdf("21ggplot.pdf", height=10, width=10)
ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
    color=umap_df$color, # color each point by WGCNA module
    size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

dev.off()
rm scRNA_harmony_singler2
rm scRNA_harmony
#save.image(file="draw_Macrophage.RData")
#save.image(file="draw_Stromal_cells.RData")
save.image(file="draw_Tissue_stem_cells.RData")