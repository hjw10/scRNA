singularity shell /public/workspace/stu21230110/design/design/scRNA_20231227.sif
R

setwd("/public/workspace/stu21230110/design/")
###加载所需要的包
library(future)
plan("multiprocess",workers=20)
###future.globals.maxSize= X,x的单位是字节，下面这句代码是8个G
options(future.globals.maxSize= 80000*1024^2) 
library(data.table)
library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
#install.packages("devtools")
#library(devtools)
#install_github("immunogenomics/harmony")
library(harmony)
#BiocManager::install("SingleCellExperiment")
library(SingleCellExperiment)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(DoubletFinder)

set.seed(123456)

if (!dir.exists("02QC")){
dir.create("02QC")
}

#双细胞函数
seurat_standard_normalize_and_scale <- function(colon, cluster, cluster_resolution){
  # colon is seurat object, 
  colon <- NormalizeData(colon, normalization.method = "LogNormalize", scale.factor = 10000)
  colon <- FindVariableFeatures(colon, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(colon)
  colon <- ScaleData(colon, features = all.genes)
  colon <- RunPCA(colon, features = VariableFeatures(object = colon))
  if (cluster){
    colon <- FindNeighbors(colon, dims = 1:20)
    colon <- FindClusters(colon, resolution = cluster_resolution)
  }
  colon <- RunUMAP(colon, dims = 1:20)
  return(colon)
}

make_seurat_object_and_doublet_removal <- function(data_directory, project_name){
  # function for basic seurat based qc and doubletfinder based doublet removal
  colon.data <- Read10X(data.dir = data_directory)
  currentSample <- CreateSeuratObject(counts = colon.data, project = project_name, min.cells = 3, min.features = 200)
  currentSample[["percent.mt"]] <- PercentageFeatureSet(currentSample, pattern = "^MT-")
  #计算核糖体基因比例
  currentSample[["percent.rb"]] <- PercentageFeatureSet(currentSample, pattern = "^Rp[sl]")
  #计算红细胞基因比例
  HB.genes <- rownames(currentSample)[grep("^Hb[^(p)]", rownames(currentSample),ignore.case = T)]

  currentSample[["percent.HB"]]<-PercentageFeatureSet(currentSample, features=HB.genes) 
# qc plot-pre filtering
  pdf(paste0("02QC/qc_plots_", project_name, "_prefiltered.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05))
  dev.off()
  pdf(paste0("02QC/qc_plots_", project_name, "_prefiltered_no_points.pdf"))
  print(VlnPlot(currentSample, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0))
  dev.off()
  
  # filter everything to 400 unique genes/cell
  minGene=200
  maxGene=6000
  maxUMI=50000
  pctMT=20
  pctHB=1

  ### 数据质控并绘制小提琴图
  currentSample <- subset(currentSample, subset = nCount_RNA < maxUMI & nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene & percent.mt < pctMT & percent.HB < pctHB)
  #currentSample <- subset(currentSample, subset = nFeature_RNA > 300   & percent.mt<20 & nCount_RNA>1000)
  
  # Normalize and make UMAP
  currentSample <- seurat_standard_normalize_and_scale(currentSample, FALSE)
  
  # Run doublet finder
  nExp_poi <- round(0.04*length(currentSample@meta.data$orig.ident)*length(currentSample@meta.data$orig.ident)/10000)  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  seu_colon <- doubletFinder_v3(currentSample, PCs = 1:20, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  print(head(seu_colon@meta.data))
  
  # rename columns
  seu_colon$doublet.class <- seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]]
  seu_colon[[paste0("DF.classifications_0.25_0.09_",nExp_poi)]] <- NULL
  pann <- grep(pattern="^pANN", x=names(seu_colon@meta.data), value=TRUE)
  seu_colon$pANN <- seu_colon[[pann]]
  seu_colon[[pann]] <- NULL
  
  # plot pre and post doublet finder results
  pdf(paste0("02QC/UMAP_pre_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", group.by = "doublet.class", cols = c("#D51F26", "#272E6A")))
  dev.off()
  seu_colon <- subset(seu_colon, subset = doublet.class != "Doublet")
  pdf(paste0("02QC/UMAP_post_double_removal", project_name, ".pdf"))
  print(DimPlot(seu_colon, reduction = "umap", cols = c("#D51F26")))
  dev.off()
  
  # Remove extra stuff and return filtered Seurat object
  seu_colon <- DietSeurat(seu_colon, counts=TRUE, data=TRUE, scale.data=FALSE, assays="RNA")
  return(seu_colon)
}

#去除双细胞
data_directory=c("./data/scRNA/C1/filtered_feature_bc_matrix/",
                 "./data/scRNA/C2/filtered_feature_bc_matrix/",
                 "./data/scRNA/C3/filtered_feature_bc_matrix/",
                 "./data/scRNA/T1/filtered_feature_bc_matrix/",
                 "./data/scRNA/T2/filtered_feature_bc_matrix/",
                 "./data/scRNA/T3/filtered_feature_bc_matrix/")

project_name<-c("C1","C2","C3","T1","T2","T3")
samples <- project_name

sample1 <- make_seurat_object_and_doublet_removal(data_directory[1], samples[1])

seu_list <- sample1
for (i in 2:length(samples)){
   
  sc.i = make_seurat_object_and_doublet_removal(data_directory[i], samples[i])
  seu_list=merge(seu_list,sc.i)
  
}

scRNA_harmony=seu_list
scRNA_harmony  <- NormalizeData(scRNA_harmony ) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
###问题 harmony之后的数据在哪里？

###一定要指定harmony
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:40) %>% FindClusters(resolution =1)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:40)
pdf("02QC/UMAP_scRNA_harmony.pdf")
DimPlot(scRNA_harmony , reduction = "umap",label = T) 
dev.off()
pdf("02QC/UMAP_scRNA_harmony_orig.ident.pdf")
DimPlot(scRNA_harmony, reduction = "umap", split.by ='orig.ident')
dev.off()
#运行时发现没有group，下面3行命令不要运行
pdf("02QC/UMAP_scRNA_harmony_group.pdf")
DimPlot(scRNA_harmony, reduction = "umap", split.by ='group',label = T)
dev.off()
pdf("02QC/UMAP_scRNA_harmony_orig.ident2.pdf")
DimPlot(scRNA_harmony, reduction = "umap", group.by='orig.ident')
dev.off()
table(scRNA_harmony$orig.ident)

scRNA_harmony_raw<-scRNA_harmony
if (!dir.exists("01RDada")){
dir.create("01RDada")
}

if (!file.exists("01RDada/01_scRNA_harmony_raw.RData")){
save(scRNA_harmony_raw,file="01RDada/01_scRNA_harmony_raw.RData")
}

library(SeuratDisk)
library(Seurat)      
#seurat2h5seurat中间过渡	
SaveH5Seurat(seurat.obj,filename="test.h5seurat", overwrite = TRUE)
#数据转为最终h5ad格式
Convert("test.h5seurat", dest = "h5ad", overwrite = TRUE)

#######################################################################################

markers <- FindAllMarkers(object = scRNA_harmony, 
                           only.pos = TRUE,
                           logfc.threshold = 0.25)   
 
DefaultAssay(scRNA_harmony) <- "RNA" 
##注：虽然这个函数的默认设置是从 "RNA "插槽中获取数据，但我们鼓励你运行上面的这行代码，
##以防止在分析的上游某处改变了默认分析。原始计数和归一化计数都存储在这个位置中，用于查找标记的函数将自动提取原始计数。 
 
#install.packages('BiocManager')
#BiocManager::install('multtest')
#install.packages('metap')
 
#colnames(scRNA_harmony@meta.data)
#marker2= FindConservedMarkers(scRNA_harmony,
#                      ident.1 = 1,
#                      grouping.var = "orig.ident",
#                      only.pos = TRUE,
#                      min.diff.pct = 0.25,
#                      min.pct = 0.25,
#                      logfc.threshold = 0.25)
 
#############################################################################################03差异基因

####单细胞转录组基础分析四：细胞类型鉴定 ####
#####如果你还想关了rstudio休息 那你就要按照这节课的开头继续设置路径、加载r包、加载数据
####一定记住 每次打开rstudio时候要先设置路径，然后加载r包，最后加数据
####细胞类型的注释一般有三种方法.1、利用marker基因查找网站进行注释  2、使用singler进行注释 3、根据已有的生物学知识或者文献，按照dotplot来注释。
##现在使用方法一寻找marker基因使用网站注释 找marker基因有以下方法三选一，建议第一种
#默认wilcox方法

##3 %>% 这是通道函数 起传递左右  可以自己百度深入理解一下，我这里只告诉你他是起传递作用的

all.markers = markers %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_7<-top10 %>% select(cluster=7)
###将marker基因保存一下
if (!dir.exists("03DEG")){
dir.create("03DEG")
}

write.csv(all.markers, "03DEG/diff_genes_wilcox.csv", row.names = F)
###把top10marker基因保存一下
write.csv(top10, "03DEG/top10_diff_genes_wilcox.csv", row.names = F)

##top10基因绘制热图
top10_genes = CaseMatch(search = as.vector(top10$gene), match = rownames(scRNA_harmony)) 
##把marker基因用热图展示出来
plot1 = DoHeatmap(scRNA_harmony, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)

ggsave("03DEG/top8_markers.pdf", plot=plot1, width=8, height=6) 

#############################################################################################04singleR
##CellMarker：http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
###PanglaoDB：https://panglaodb.se/index.html
####把每一类的marker基因输入到这两个网站中，然后看网站给出的细胞注释信息
##选择基因作图展示-----文章中某些图需要
#挑选部分基因
#select_genes <- c('COL2A1',"PTPRC","EPCAM")
#vlnplot展示
#p1 <- VlnPlot(scRNA_harmony, features = select_genes, pt.size=0, group.by="celltype", ncol=2)
#p1
#ggsave("selectgenes_VlnPlot.png", p1, width=6 ,height=8)
#featureplot展示
#p2 <- FeaturePlot(scRNA_harmony, features = select_genes, reduction = "tsne", label=T, ncol=2)
#p2
#ggsave("selectgenes_FeaturePlot.png", p2, width=8 ,height=12)
#p3=p1|p2
#p3
#ggsave("selectgenes.png", p3, width=10 ,height=8)

#第二种方法用SingleR鉴定细胞类型
#BiocManager::install("SingleR")
library(SingleR)
###下载好数据库后，把ref_Human_all.Rdata加载到环境中，这样算是对数据库的加载，就可以按照singler的算法来对细胞亚群进行定义了。
load("ref_Human_all.RData")
###我们可以看到在环境中多了一个叫ref_Human_all的文件 大小为113mb  这个就是数据库
####然后我们把环境中的ref_Human_all赋值与refdata
refdata <- ref_Human_all
###把rna的转录表达数据提取

testdata <- GetAssayData(scRNA_harmony, slot="data")
###把scRNA数据中的seurat_clusters提取出来，注意这里是因子类型的
clusters <- scRNA_harmony@meta.data$seurat_clusters
###开始用singler分析
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
###制作细胞类型的注释文件
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = FALSE)
###保存一下
if (!dir.exists("04singleR")){
dir.create("04singleR")
}
write.csv(celltype,"04singleR/celltype_singleR.csv",row.names = FALSE)
##把singler的注释写到metadata中 有两种方法
###方法一----没有用
#scRNA1@meta.data$celltype = "NA"
#for(i in 1:nrow(celltype)){
#  scRNA1@meta.data[which(scRNA1@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
###因为我把singler的注释加载到metadata中时候，命名的名字叫celltype，所以画图时候，group.by="celltype"
#DimPlot(scRNA1, group.by="celltype", label=T, label.size=5, reduction='tsne')
###方法二：-----用的此方法
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F) 
scRNA_harmony@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
###因为我把singler的注释加载到metadata中时候，命名的名字叫singleR，所以画图时候，group.by="singleR"
#DimPlot(scRNA_harmony, group.by="singleR", label=T, label.size=5, reduction='tsne')
###我们可以看到  两种方法得到的结果都是一样的，但是我比较喜欢第二种方法
p<-DimPlot(scRNA_harmony, group.by="singleR", label=T, label.size=3.5, reduction='umap',split.by="orig.ident")
ggsave("04singleR/umap_celltype_singler.pdf", p, width=30 ,height=6)

p<-DimPlot(scRNA_harmony, group.by="singleR", label=T, label.size=3.5, reduction='umap')
ggsave("04singleR/umap_celltype_singler_all.pdf", p, width=10 ,height=6)


##########################################################################05HCL包注释
library(scHCL)
if (!dir.exists("05scHCL")){
dir.create("05scHCL")
}
hcl_result <- scHCL(scdata = scRNA_harmony@assays$RNA@data, numbers_plot = 3)
save(hcl_result,file="05scHCL/01scRNA_harmony_scHCL.RData")
scRNA_harmony@meta.data$scHCL<-as.factor(hcl_result$scHCL)

p<-DimPlot(scRNA_harmony, group.by="scHCL", label=T, label.size=3.5, reduction='umap')
####由于注释的细胞类型很多，图片大小需要调整，并且加入limitsize = FALSE参数，要不然会报错
ggsave("05scHCL/01umap_celltype_scHCL2.pdf", p, width=120 ,height=10,limitsize = FALSE)
#结合scMCA结果，修改singleR注释结果，修改celltype
save.image(file="01RDada/celltype_scHCL.RData")
p<-FeaturePlot(scRNA_harmony, features = c("Il33"), split.by = "orig.ident")
ggsave("IL33-umap_celltype_group.pdf", p, width=10,height=6)

p<-DimPlot(scRNA_harmony, group.by="scHCL", label=T, label.size=3.5, reduction='umap')
####由于注释的细胞类型很多，图片大小需要调整，并且加入limitsize = FALSE参数，要不然会报错
ggsave("05scMCA/01umap_celltype_sHCL3.pdf", p, width=120 ,height=30,limitsize = FALSE)
##########
#结合scHCL结果，修改singleR注释结果，修改celltype，

###颜色
sample_color <- c("#BCBD22FF","#4DBBD5FF","#00A087FF","#91D1C2FF","#F39B7FFF","#F7B6D2FF","#7E6148FF","#E64B35FF","#3C5488FF","#E377C2FF","#8491B4FF","#B09C85FF")
sample_color <- c("#F39B7FFF","#E377C2FF","#3C5488FF","#7E6148FF","#00A087FF","#8494FF","#4DBBD5FF", "#91D1C2FF","#9467BDFF","#E64B35FF","#B09C85FF")
load("01RDada/06_scRNA_harmony_singler202400810.Rdata")
pdf("04singleR/umap_celltype_singler_20240810.pdf")
DimPlot(scRNA_harmony_singler, group.by="Annotation", label=T, label.size=3.5, reduction='umap',cols=sample_color)
dev.off()
