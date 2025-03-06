######################################################chat分组重新跑
########对照之前须先分别跑正常和重症组
#################################################################################10 首先跑对照组
#cellchat   20240102重新跑screen -S hc11
#单数据集
load("01RDada/06_scRNA_harmony_singler202400810.Rdata")
####
scRNA_harmony<-scRNA_harmony_singler
library(Seurat)#载入R包
library(ggplot2)
library(CellChat)
library(patchwork)
library(NMF)
if (!dir.exists("07CellChat")){
dir.create("07CellChat")
}

setwd("/public/workspace/stu21230110/design/07CellChat/C")
# 创建cellchat对象
seurat=subset(scRNA_harmony,group=="C")#读取seurat对象
data.input = as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))#得到标准化表达矩阵
meta = seurat@meta.data # 获得实验信息
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleRnew")#创建cellchat对象

levels(cellchat@idents) # 显示细胞类型
groupSize <- as.numeric(table(cellchat@idents)) # 每种细胞类型的细胞数目统计
groupSize

CellChatDB <- CellChatDB.human # 若是小鼠，请用 CellChatDB.mouse 数据库
showDatabaseCategory(CellChatDB)#显示数据库中的配受体分类

cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))#丢弃未使用的因子水平，在我们这个例子数据中需要该步骤


cellchat@DB <- CellChatDB #数据库赋值
cellchat <- subsetData(cellchat) # 必须，即时是整个数据集
cellchat <- identifyOverExpressedGenes(cellchat)#鉴定过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#鉴定过表达互作
cellchat <- projectData(cellchat, PPI.human)#将基因表达数据映射PPI互作网络

#################################Inference of cell-cell communication network#####################################
cellchat <- computeCommunProb(cellchat)#计算互作概率
cellchat <- filterCommunication(cellchat, min.cells = 1)#过滤掉低于min.cells的细胞类型
cellchat <- computeCommunProbPathway(cellchat)#计算pathway水平的互作概率，通过总结所有的配受体对
cellchat <- aggregateNet(cellchat)#得到整体细胞互作网络
save(cellchat,file="cellchat_hc.RData")
cellchat@netP$pathways#显示显著的pathways
######一定要看结果

[1] "COLLAGEN"   "MK"         "LAMININ"    "MIF"        "APP"
 [6] "MHC-I"      "FN1"        "MHC-II"     "CD99"       "PARs"
[11] "CLEC"       "PTN"        "CADM"       "PTPRM"      "GALECTIN"
[16] "TENASCIN"   "VISFATIN"   "THY1"       "CCL"        "NRXN"
[21] "GAS"        "CD45"       "SEMA3"      "SPP1"       "CD46"
[26] "ITGB2"      "ICAM"       "IGF"        "TGFb"       "JAM"
[31] "NOTCH"      "PECAM1"     "THBS"       "VEGF"       "TNF"
[36] "ncWNT"      "EDN"        "CDH"        "BMP"        "EPHA"
[41] "ADGRE5"     "HH"         "SEMA4"      "ANNEXIN"    "PDGF"
[46] "COMPLEMENT" "ALCAM"      "CD6"        "ANGPT"      "CDH5"
[51] "CDH1"       "CXCL"       "ESAM"       "NECTIN"     "HSPG"
[56] "ANGPTL"     "RESISTIN"   "MPZ"        "CD40"       "CHEMERIN"
[61] "IFN-II"     "EGF"        "BAFF"       "WNT"        "APRIL"
[66] "EPHB"       "PROS"       "NRG"        "IL16"       "SEMA6"
[71] "CNTN"       "OCLN"







pdf("01-1count_weight.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")#互作数目
dev.off()

pdf("01-2count_weight.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")#互作强度
dev.off()







mat <- cellchat@net$weight#互作强度值
pdf("02allcell_count_weight.pdf")
par(mfrow = c(3,4), xpd=TRUE) # 设定绘图区，3行4列
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))#构建1个0值矩阵
  mat2[i, ] <- mat[i, ]#替换掉矩阵对应行
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])#网络图
}
dev.off()
# 层次图
levels(cellchat@idents) 

 [1] "B_cells"                  "Endothelial_cells"
 [3] "Epithelial_cells"         "Macrophage"
 [5] "Monocyte"                 "Myeloid_progenitor_cells"
 [7] "NK_cells"                 "Smooth_muscle_cells"
 [9] "Stromal_cells"            "T_cells"
[11] "Tissue_stem_cells"

pathways.show <- "CCL"#指定信号通路
vertex.receiver = 6 # 指定目标细胞类型
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "circle")#绘图
# 所有显著pathways
pathways.show.all <- cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  # 信号pathway以及对应的每个配受体对可视化，网络图，弦图，热图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "circle",out.format=c("pdf"))#网络图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "chord",out.format=c("pdf"),height=12)#弦图
  pdf(file=paste0("03-1",pathways.show.all[i], "_heatmap.pdf"))#热图
  print(netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds"))
  dev.off()
  # 计算每个配受体对对信号通路的贡献大小
  print("LR contribution!")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])#绘图
  ggsave(filename=paste0("03-2",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
}

#对感兴趣的细胞类型，气泡图显示所有显著LR对
pdf("03-0netVisual_bubble_Macrophage.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 4, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Stromal_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 9, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Tissue_stem_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 11, remove.isolate = FALSE)
dev.off()
#绘制感兴趣信号通路相关配受体对的基因表达图

pdf("03-1-plotGeneExpression.pdf",width=8,height=10)
plotGeneExpression(cellchat, signaling = "CCL")
dev.off()


# 网络分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 计算信号通路的网络中心得分
# 热图可视化细胞所扮演的主要角色，亦可用散点图来表示
for(i in 1:length(pathways.show.all)){
  pdf(file=paste0("04-1",pathways.show.all[i], "_major_signaling_roles_heatmap.pdf"),width=10,height=8)#热图
  print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10))
  dev.off()
  pdf(file=paste0("04-2",pathways.show.all[i], "_major_signaling_roles_scatter.pdf"),width=8,height=8)#散点图
  print(netAnalysis_signalingRole_scatter(cellchat,signaling = pathways.show.all[i]))
  dev.off()
}

# 所有信号通路汇总分析，看哪个信号通路贡献最大
pdf("05netAnalysis_signalingRole_heatmap.pdf",height=32,width=18)
par(mfrow = c(1,2))#
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width=10,height=30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width=10,height=30)


ht1 + ht2
dev.off()
#保存结果
df.net <- subsetCommunication(cellchat)#得到细胞互作数据框
write.csv(df.net,"06-net_lr.csv")#输出文件
df.netp <- subsetCommunication(cellchat,slot.name = "netP")#得到信号通路水平的细胞互作结果
write.csv(df.netp,"07-net_pathway.csv")#输出文件
cellchat_HC<-cellchat
saveRDS(cellchat_HC, file = "HC_cellchat_single_dataset.rds")#保存cellchat对象为rds文件




##########################################################################################11 跑重症组
#cellchat  20221026重新跑screen -S pd11
#单数据集
load("01RDada/06_scRNA_harmony_singler202400810.Rdata")
####
scRNA_harmony<-scRNA_harmony_singler
library(Seurat)#载入R包
library(ggplot2)
library(CellChat)
library(patchwork)
library(NMF)
if (!dir.exists("11chat_pd")){
dir.create("11chat_pd")
}

setwd("/public/workspace/stu21230110/design/07CellChat/P")
# 创建cellchat对象
seurat=subset(scRNA_harmony,group=="P")#读取seurat对象
data.input = as.matrix(GetAssayData(seurat, assay = "RNA", slot = "data"))#得到标准化表达矩阵
meta = seurat@meta.data # 获得实验信息
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleRnew")#创建cellchat对象

levels(cellchat@idents) # 显示细胞类型
groupSize <- as.numeric(table(cellchat@idents)) # 每种细胞类型的细胞数目统计
groupSize

CellChatDB <- CellChatDB.human # 若是小鼠，请用 CellChatDB.mouse 数据库
#showDatabaseCategory(CellChatDB)#显示数据库中的配受体分类

cellchat@idents = droplevels(cellchat@idents, exclude = setdiff(levels(cellchat@idents),unique(cellchat@idents)))#丢弃未使用的因子水平，在我们这个例子数据中需要该步骤


cellchat@DB <- CellChatDB #数据库赋值
cellchat <- subsetData(cellchat) # 必须，即时是整个数据集
cellchat <- identifyOverExpressedGenes(cellchat)#鉴定过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat)#鉴定过表达互作
cellchat <- projectData(cellchat, PPI.human)#将基因表达数据映射PPI互作网络

#################################Inference of cell-cell communication network#####################################
cellchat <- computeCommunProb(cellchat)#计算互作概率
cellchat <- filterCommunication(cellchat, min.cells = 1)#过滤掉低于min.cells的细胞类型
cellchat <- computeCommunProbPathway(cellchat)#计算pathway水平的互作概率，通过总结所有的配受体对
cellchat <- aggregateNet(cellchat)#得到整体细胞互作网络
save(cellchat,file="cellchat_P.RData")
cellchat@netP$pathways#显示显著的pathways
######一定要看结果



 [1] "COLLAGEN"   "MIF"        "LAMININ"    "MK"         "APP"
 [6] "MHC-I"      "FN1"        "CD99"       "MHC-II"     "PARs"
[11] "VISFATIN"   "CLEC"       "PTPRM"      "ADGRE5"     "CD45"
[16] "PTN"        "TENASCIN"   "GALECTIN"   "ITGB2"      "THY1"
[21] "CCL"        "ICAM"       "PECAM1"     "NRXN"       "IGF"
[26] "VEGF"       "ANNEXIN"    "CXCL"       "THBS"       "TGFb"
[31] "SPP1"       "CD22"       "JAM"        "EPHA"       "CD46"
[36] "CDH"        "NOTCH"      "ncWNT"      "CADM"       "GAS"
[41] "EGF"        "BMP"        "PDGF"       "IFN-II"     "NECTIN"
[46] "SEMA3"      "ANGPTL"     "CDH1"       "EDN"        "ANGPT"
[51] "ALCAM"      "CD6"        "CDH5"       "SEMA4"      "CD226"
[56] "ESAM"       "SELL"       "COMPLEMENT" "TIGIT"      "GRN"
[61] "NCAM"       "PROS"       "MPZ"        "EPHB"       "WNT"
[66] "FGF"        "HGF"        "BAFF"       "CHEMERIN"   "NRG"
[71] "SEMA6"      "SELPLG"     "SN"         "OCLN"







pdf("01-1count_weight.pdf")
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")#互作数目
dev.off()

pdf("01-2count_weight.pdf")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")#互作强度
dev.off()

mat <- cellchat@net$weight#互作强度值
pdf("02allcell_count_weight.pdf")
par(mfrow = c(3,4), xpd=TRUE) # 设定绘图区，3行4列
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))#构建1个0值矩阵
  mat2[i, ] <- mat[i, ]#替换掉矩阵对应行
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])#网络图
}
dev.off()
# 层次图
levels(cellchat@idents) 

 [1] "B_cells"                  "Endothelial_cells"
 [3] "Epithelial_cells"         "Macrophage"
 [5] "Monocyte"                 "Myeloid_progenitor_cells"
 [7] "NK_cells"                 "Smooth_muscle_cells"
 [9] "Stromal_cells"            "T_cells"
[11] "Tissue_stem_cells"

pathways.show <- "CCL"#指定信号通路
vertex.receiver = 6 # 指定目标细胞类型
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "circle")#绘图
# 所有显著pathways
pathways.show.all <- cellchat@netP$pathways

for (i in 1:length(pathways.show.all)) {
  # 信号pathway以及对应的每个配受体对可视化，网络图，弦图，热图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "circle",out.format=c("pdf"))#网络图
  netVisual(cellchat, signaling = pathways.show.all[i], layout = "chord",out.format=c("pdf"),height=12)#弦图
  pdf(file=paste0("03-1",pathways.show.all[i], "_heatmap.pdf"))#热图
  print(netVisual_heatmap(cellchat, signaling = pathways.show.all[i], color.heatmap = "Reds"))
  dev.off()
  # 计算每个配受体对对信号通路的贡献大小
  print("LR contribution!")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])#绘图
  ggsave(filename=paste0("03-2",pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
}

#对感兴趣的细胞类型，气泡图显示所有显著LR对
pdf("03-0netVisual_bubble_Macrophage.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 4, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Stromal_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 9, remove.isolate = FALSE)
dev.off()
pdf("03-0netVisual_bubble_Tissue_stem_cells.pdf",width=10,height=20)
netVisual_bubble(cellchat, sources.use = c(1:5,7:11), targets.use = 11, remove.isolate = FALSE)
dev.off()
#绘制感兴趣信号通路相关配受体对的基因表达图

pdf("03-1-plotGeneExpression.pdf",width=8,height=10)
plotGeneExpression(cellchat, signaling = "CCL")
dev.off()


# 网络分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # 计算信号通路的网络中心得分
# 热图可视化细胞所扮演的主要角色，亦可用散点图来表示
for(i in 1:length(pathways.show.all)){
  pdf(file=paste0("04-1",pathways.show.all[i], "_major_signaling_roles_heatmap.pdf"),width=10,height=8)#热图
  print(netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10))
  dev.off()
  pdf(file=paste0("04-2",pathways.show.all[i], "_major_signaling_roles_scatter.pdf"),width=8,height=8)#散点图
  print(netAnalysis_signalingRole_scatter(cellchat,signaling = pathways.show.all[i]))
  dev.off()
}

# 所有信号通路汇总分析，看哪个信号通路贡献最大
pdf("05netAnalysis_signalingRole_heatmap.pdf",height=32,width=18)
par(mfrow = c(1,2))#
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width=10,height=30)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width=10,height=30)


ht1 + ht2
dev.off()
#保存结果
df.net <- subsetCommunication(cellchat)#得到细胞互作数据框
write.csv(df.net,"06-net_lr.csv")#输出文件
df.netp <- subsetCommunication(cellchat,slot.name = "netP")#得到信号通路水平的细胞互作结果
write.csv(df.netp,"07-net_pathway.csv")#输出文件
cellchat_P<-cellchat
saveRDS(cellchat_P, file = "P_cellchat_single_dataset.rds")#保存cellchat对象为rds文件

#########################################################12cellchat多数据集
#拥有相同细胞构成的多个数据集比较分析
library(CellChat)#载入R包
library(patchwork)
setwd("/public/workspace/stu21230110/design/07CellChat")
if (!dir.exists("chat_multi")){
dir.create("chat_multi")
}
setwd("./chat_multi")#切换工作目录


cellchat.C <- readRDS("../C/HC_cellchat_single_dataset.rds")#读取ctrl的cellchat对象
cellchat.P <- readRDS("../P/P_cellchat_single_dataset.rds")#读取treat的cellchat对象
object.list <- list(C= cellchat.C, P= cellchat.P)#创建列表
cellchat <- mergeCellChat(object.list, add.names = names(object.list))#合并cellchat对象

pdf(file="01num_weight.pdf",width=16,height=8)
#互作数目和强度比较
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))#互作数目
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")#互作强度
gg1 + gg2
dev.off()
#circle plot，网络图
pdf(file="02circle_plot.pdf",width=16,height=8)
par(mfrow = c(1,2), xpd=TRUE)#设置绘图区域,1行2列
netVisual_diffInteraction(cellchat, weight.scale = T)#互作数目
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")#互作强度
dev.off()
#heatmap，热图
pdf(file="03heatmap.pdf",width=16,height=8)
gg1 <- netVisual_heatmap(cellchat)#互作数目
gg2 <- netVisual_heatmap(cellchat, measure = "weight")#互作强度
gg1 + gg2
dev.off()
#上述差异分析只适合两组比较，如果超过两组，可以采用下面的方式来比较分析。
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))#获得多个数据集间最大的边权重
#par(mfrow = c(1,2), xpd=TRUE)#设置绘图区域，1行2列
for (i in 1:length(object.list)) {#依次绘制每个数据对象的互作网络图
  p<-netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
  #ggsave(filename=paste0("04-",i, "_network.pdf"), plot=p, width = 8, height = 8, units = 'in', dpi = 300)#保存图片
  pdf(file=paste0("04-",i, "_network.pdf"), width = 8, height = 8)
  p
  dev.off()
}

#发送（outgoing）和接收（incoming）信号角色比较
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})#计算互作数目
weight.MinMax <- c(min(num.link), max(num.link)) # 获取互作数目的最小值和最大值，用来标准化图形的点大小
#gg <- list()#新建一个空列表
for (i in 1:length(object.list)) {#绘制气泡图
  p <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  ggsave(filename=paste0("05-",i, "_dotplot.pdf"), plot=p, width = 8, height = 8, units = 'in', dpi = 300)#保存图片

}
#patchwork::wrap_plots(plots = gg)
#对感兴趣细胞类型的信号通路变化进行分析
#gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "CD4 T")
#gg1
#两组信号通路差异分析
pdf(file="06deg.pdf",width=16,height=8)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()
#发送（outgoing）和接收（incoming）信号通路比较
library(ComplexHeatmap)
i = 1
#合并所有数据集鉴定的信号通路 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 30, height = 25)#第一个对象热图
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 30, height = 25)#第二个对象热图
pdf(file="07two.pdf",width=40,height=26)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))#合并两个绘图对象
dev.off()
#配受体对差异分析
#netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:10),  comparison = c(1, 2), angle.x = 45)
#分别显示互作增加或减少
target <- c(4,9,11)
pdf(file="08interaction_deg_Macrophage.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 4,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 4,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08interaction_deg_Stromal_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 9,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 9,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08interaction_deg_Tissue_stem_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 11,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:4,6:10), targets.use = 11,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()

pdf(file="08-2interaction_deg_Macrophage.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08-2interaction_deg_Stromal_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()
pdf(file="08-2interaction_deg_Tissue_stem_cells.pdf",width=16,height=12)
gg1 <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in patient", angle.x = 45, remove.isolate = T)#互作增加
gg2 <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:4,6:10),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in patient", angle.x = 45, remove.isolate = T)#互作减少
gg1 + gg2
dev.off()

gg1$data#得到配受体对数据
gg2$data
#画基因表达
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("C", "P")) # 设置因子水平
pdf(file="09gene.pdf",width=16,height=8)
plotGeneExpression(cellchat, signaling = "CCL", split.by = "datasets", colors.ggplot = T)#画基因表达图
dev.off()
#保存rds文件
saveRDS(cellchat, file = "cellchat_comparisonAnalysis_P_vs_C.rds")