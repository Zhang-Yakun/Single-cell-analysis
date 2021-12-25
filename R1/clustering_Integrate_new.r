suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(BiocParallel))
suppressMessages(require(BiocNeighbors))
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)

# Seurat (anchors and CCA) ------------------------------------------------

setwd("e:/Rwork/GSE144945/")

# 批量读取数据 --------------------------------------------------------------------

rm(list=ls())

assays <- dir("e:/Rwork/GSE144945/single cell data/")
dir <- paste0("e:/Rwork/GSE144945/single cell data/", assays)
#The vector that stores the path

samples_name <- paste0("P",1:10)
# Sample ID


# 批量创建seurat对象 ------------------------------------------------------------

scRNAlist <- list()

for(i in 1:length(dir)){
  
  counts <- Read10X(data.dir = dir[i])
  
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = samples_name[i],
                                       min.cells = 3)
  #给细胞barcode加前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = samples_name[i])
  
  ## QC
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  #计算线粒体基因比例
  
  # filter
  scRNAlist[[i]]<-subset(scRNAlist[[i]],subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)
  
  #NormalizeData
  scRNAlist[[i]]<-NormalizeData(scRNAlist[[i]],normalization.method = "LogNormalize",scale.factor = 10000)
  
  #FindVariableFeatures
  scRNAlist[[i]]<-FindVariableFeatures(scRNAlist[[i]],selection.method="vst",nfeatures=2000)
  
}
names(scRNAlist) <- samples_name


# 使用FindIntegrationAnchors函数来识别锚点（anchors），该函数的输入数据是Seurat对象的列表。 ----------

test.list <- scRNAlist[c("P1", "P2", "P3", "P4","P5","P6","P7","P8","P9","P10")]
test.anchors <- FindIntegrationAnchors(object.list = test.list, dims = 1:30)


# 然后将这些锚（anchors）传递给IntegrateData函数，该函数返回Seurat对象。 ------------------------

test.integrated <- IntegrateData(anchorset = test.anchors, dims = 1:30)

table(test.integrated$orig.ident)
#统计每个patient的细胞数

saveRDS(test.integrated,file = "e:/Rwork/GSE144945/Integrate_Anchor/integrated.rds")
#save(test.integrated, file = "scRNA_anchors.Rdata")
#保存整合后的scRNA-Seurat

test.integrated<-readRDS("e:/Rwork/GSE144945/Integrate_Anchor/integrated.rds")

test.integrated@assays
##运行IntegrateData之后，Seurat对象将包含一个具有整合（或“批次校正”）表达矩阵的新Assay。请注意，原始值（未校正的值）仍存储在“RNA”分析的对象中，因此可以来回切换。

DefaultAssay(test.integrated) <- "integrated"


# 运行标准流程并进行可视化 ------------------------------------------------------------
setwd("e:/Rwork/GSE144945/Integrate_Anchor")

########  VariableFeature

test.integrated<-FindVariableFeatures(test.integrated,selection.method="vst",nfeatures=2000)
top10<-head(VariableFeatures(test.integrated),10)
plot1<-VariableFeaturePlot(test.integrated)
#散点图，X轴：平均表达值；Y轴：标准化后的方差，显示前2000个高变基因

plot2<-LabelPoints(plot = plot1,points = top10,repel = TRUE)
#在散点图plot1中加入基因名

plot1+plot2


###### Scale & PCA

test.integrated <- ScaleData(test.integrated, verbose = FALSE)
test.integrated <- RunPCA(test.integrated, npcs = 30, verbose = FALSE)

VizDimplot<-VizDimLoadings(test.integrated,dims = 1:2,reduction="pca")

dir.create("PCA")
ggsave("PCA/PC_1&2_genes.pdf",plot = VizDimplot)

PCAplot<-DimPlot(test.integrated,reduction = "pca",group.by = "orig.ident")
#可视化PCA聚类
ggsave("PCA/PCA_batch.pdf",plot = PCAplot)

PCA_heatmap<-DimHeatmap(test.integrated,dims = 1:15,cells = 500,balanced = TRUE)
#PC1:15的基因在500个细胞中表达值的展示
ggsave("PCA/PC15_heatmap.pdf",plot = PCA_heatmap)


###### clustering

test.integrated<-FindNeighbors(test.integrated,dims=1:20)
test.integrated<-FindClusters(test.integrated,resolution = 0.5)

head(Idents(test.integrated),5)
#查看每个细胞对应的聚类数ID

head(test.integrated@meta.data)
table(test.integrated@meta.data$seurat_clusters)
#查看每个cluster中有多少细胞，cluster—0中的细胞数最多

test.integrated<-RunTSNE(test.integrated,dims=1:20)
test.integrated<-RunUMAP(test.integrated,dims = 1:20)

tsne_plot<-DimPlot(test.integrated,reduction = "tsne",group.by = "orig.ident")
umap_plot<-DimPlot(test.integrated,reduction = "umap",group.by = "orig.ident")

dir.create("Clustering")
ggsave("Clustering/tsne.pdf",plot = tsne_plot)
ggsave("Clustering/umap.pdf",plot = umap_plot)

tsne_plot2<-DimPlot(test.integrated,reduction = "tsne",label = T,label.size=5)
tsne_plot3<-DimPlot(test.integrated,reduction = "tsne",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/tsne2.pdf",plot = tsne_plot2)
ggsave("Clustering/tsne3.pdf",plot = tsne_plot3)

umap_plot2<-DimPlot(test.integrated,reduction = "umap",label = T,label.size=5)
umap_plot3<-DimPlot(test.integrated,reduction = "umap",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/umap2.pdf",plot = umap_plot2)
ggsave("Clustering/umap3.pdf",plot = umap_plot3)


#######UMAP
test.integrated <- RunUMAP(test.integrated, reduction = "pca", dims = 1:20)
p1 <- DimPlot(test.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(test.integrated, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p1, p2)

########STNE
test.integrated <- RunTSNE(test.integrated, reduction = "pca", dims = 1:20)
p3 <- DimPlot(test.integrated, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(test.integrated, reduction = "tsne", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p3, p4)

saveRDS(test.integrated,file = "e:/Rwork/GSE144945/Integrate_Anchor/integrated_cluster.rds")

# DEG analysis ---------------------------------------------------------

##test.integrated<-readRDS("e:/Rwork/GSE144945/Integrate_Anchor/integrated_cluster.rds")

#识别差异基因
scRNA.pos.markers<-FindAllMarkers(test.integrated,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)


#only.pos=TRUE表示只保留上调的差异基因

write.csv(scRNA.pos.markers,file = "e:/Rwork/GSE144945/Integrate_anchor/DEG analysis/FindAllMarkers/Allmarkers-pos-C17.csv")
#该markers文件用于后边的细胞类型注释
head(scRNA.pos.markers)

Top_pos_genes<-scRNA.pos.markers%>%group_by(cluster)%>%top_n(n=1,wt=avg_log2FC)
#找出各个cluster中的TOP差异基因，按logFC值排序，n=2展示各个cluster的TOP2

# 展示候选marker基因的表达量 --------------------------------------------------------
#从上一部分“输出marker基因”输出的cell_type_markers.csv文件里选出每类p_val_adj排名第一的基因，例如第II类(class为1)的ENSG00000205364。

library("RColorBrewer")
display.brewer.all(type = "qual")

top<-Top_pos_genes$gene%>%unique()

FeaturePlot(test.integrated, reduction = "tsne",
            features = top[1:9], 
            ncol  = 3, #画在2列里
            cols = c("#66C2A5","grey","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            #do.return =T, #返回ggplot2对象，便于组图
            #no.legend = F, #显示图例
            pt.size = 1
) 
ggsave("pos_marker_C0-C8.pdf")

# 细胞类型手动注释Cell types annotation ------------------------------------------------------------------

#读入参考marker文件
ref_marker<-fread("E:/Rwork/GSE144945/Integrate_anchor/cell_type_annotation/RGEP_ref_marker.txt")%>%as.data.frame()
colnames(ref_marker)<-c("CellType","Markers")

#cluster-specific up-markers
pbmc.markers<-fread("e:/Rwork/GSE144945/Integrate_anchor/DEG analysis/FindAllMarkers/Allmarkers-pos-C17.csv")

######超几何检验#######


p<-data.frame()

for(i in 1:length(unique(pbmc.markers$cluster))){
  
  cluster_gene<-pbmc.markers$gene[pbmc.markers$cluster==unique(pbmc.markers$cluster)[i]]
  #cluster_i中特异上调的基因向量
  
  for(j in 1:dim(ref_marker)[1]){
    
    mark<-strsplit(ref_marker$Markers[j],",")[[1]]
    cell_type<-rep(ref_marker$CellType[j],length(mark))
    ref_df<-data.frame(CellType=cell_type,MarkerGene=mark)
    #cell_type_j中的marker基因
    
    N<-union(strsplit(ref_marker$Markers,",")%>%unlist(),pbmc.markers$gene)%>%length()
    #背景基因数
    q<-intersect(cluster_gene,ref_df$MarkerGene)%>%length()
    #celltype1和cluster1的交集基因数
    m<-length(ref_df$MarkerGene)
    #celltype中的基因数
    n<-N-m
    k<-length(cluster_gene)
    #cluster中的基因数
    
    p_value<-1-phyper(q,m,n,k)
    
    df_p<-data.frame(cluster=unique(pbmc.markers$cluster)[i],
                     celltype=ref_marker$CellType[j],
                     phyper_pvalue=p_value)
    p<-rbind(p,df_p)
    
  }
}

annotation_p<-subset(p,p$phyper_pvalue<0.5)

test<-p %>% group_by(cluster) %>% dplyr::mutate(min_rank(phyper_pvalue))
annotation_top<-subset(test,test$`min_rank(phyper_pvalue)`==1)
annotation_top$cluster<-paste(rep("C",length(annotation_top$cluster)),annotation_top$cluster,sep = "")

annotation_top$celltype<-paste(annotation_top$celltype,annotation_top$cluster,sep = "_")

#####重新定义细胞类型######

new.cluster.ids<-annotation_top$celltype

names(new.cluster.ids)<-levels(test.integrated)

Idents(test.integrated)
#原来细胞聚类的名称

test.integrated<-RenameIdents(test.integrated,new.cluster.ids)
Idents(test.integrated)
head(test.integrated@meta.data)
#更改后的聚类名

DimPlot(test.integrated,reduction = "tsne",pt.size = 0.5)
DimPlot(test.integrated,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
#可视化
saveRDS(test.integrated,file = "E:/Rwork/GSE144945/Integrate_Anchor/cell_type_annotation/scRNA_Annotation.rds")


# 统计各类免疫细胞的数量 -------------------------------------------------------------

scRNA.integrated<-readRDS("E:/Rwork/GSE144945/Integrate_Anchor/cell_type_annotation/scRNA_Annotation.rds")

#每个细胞的免疫类型
cel.type<-as.character(Idents(scRNA.integrated))
cell.type<-str_split(cel.type,"_",simplify = TRUE)
cells<-table(cell.type[,1])%>%as.data.frame()
colnames(cells)<-c("cluster","number")
cells$id<-1:length(cells$cluster)

#cells是circos的输入


# 绘制免疫细胞的数量circos_barplot -------------------------------------------------

cir_bar_data<-cells
# 绘制基础环状条形图

label_data <- cir_bar_data

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
#查看标签数据
head(label_data)

# 绘制带标签的环状条形图
p <- ggplot(label_data, aes(x=as.factor(id), y=number)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("skyblue", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-2000,10000) +
  
  # Custom the theme: no axis title and no cartesian grid
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
  ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=number+100, label=cluster, hjust=hjust), 
            color="black", fontface="bold",alpha=0.6, size=2.5, 
            angle= label_data$angle, inherit.aes = FALSE ) 

p2 <- ggplot(label_data, aes(x=as.factor(id), y=number, fill=cluster)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(stat="identity", alpha=0.5) +
  ylim(-2000,10000) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=number+100, label=cluster, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data$angle, inherit.aes = FALSE ) 

