library(Seurat)
library(SeuratData)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggridges)
library(dplyr)
library(data.table)
library(dplyr)
library(ClusterProfiler)

setwd("e:/Rwork/GSE144945/")

# 批量读取数据 --------------------------------------------------------------------

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
  
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = samples_name[i])
  #给细胞barcode加前缀，防止合并后barcode重名
 
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  #计算线粒体基因比例
 
}

###给列表命名，并保存数据
dir.create("Integrate")
setwd("./Integrate")

names(scRNAlist) <- samples_name
system.time(saveRDS(scRNAlist,file = "scRNAlist.rds"))


# 用merge()将scRNAlist合成一个Seurat对象 ------------------------------------------

scRNA <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])

table(scRNA$orig.ident)
#统计每个patient的细胞数

save(scRNA, file = "scRNA_merge.Rdata")
#保存整合后的scRNA-Seurat

load("E:/Rwork/GSE144945/Integrate/scRNA_merge.Rdata")
# QC ----------------------------------------------------------------------

##绘制质控指标的小提琴图

#设置可能用到的主题
theme.set2=theme(axis.title.x = element_blank())
#设置绘图元素
plot.features=c("nFeature_RNA","nCount_RNA","percent.mt")
group="orig.ident"

#QC前的小提琴图
plots=list()

#seq_along()创建一个从1开始，长度为其输入值的序列
for(i in seq_along(plot.features)){
  
  plots[[i]]=VlnPlot(scRNA,group.by = group,pt.size = 0,
                     features = plot.features[i])+theme.set2+NoLegend()
}

violin<- wrap_plots(plots=plots,ncol=3)
violin
dir.create("QC")
ggsave("QC/vlnplot_before_qc.pdf",plot = violin, width = 9,height = 8)

##质量控制，绘制QC后的小提琴图
scRNA<- subset(scRNA, subset=nFeature_RNA>200 & nFeature_RNA<2500 & percent.mt<5)

plots=list()

for(i in seq_along(plot.features)){
  
  plots[[i]]=VlnPlot(scRNA,group.by = group,pt.size = 0,
                     features = plot.features[i])+theme.set2+NoLegend()
}

violin<-wrap_plots(plots=plots,ncol = 3)
violin
ggsave("QC/vlnplot_after_qc.pdf",plot = violin,width = 10,height = 8)


# 查看批次效应（对merge后的Seurat对象进行标准化和降维） ------------------------------------------------------------------

#标准化，高变特征筛选，归一化
scRNA<-NormalizeData(scRNA) %>% FindVariableFeatures(selection.method="vst",nfeatures=2000) %>% ScaleData()

#PCA
scRNA<-RunPCA(scRNA,features=VariableFeatures(scRNA))

#确定聚类时用的PC数
VizDimplot<-VizDimLoadings(scRNA,dims = 1:2,reduction="pca")

dir.create("PCA")
ggsave("PCA/PC_1&2_genes.pdf",plot = VizDimplot)

PCAplot<-DimPlot(scRNA,reduction = "pca",group.by = "orig.ident")
#可视化PCA聚类
ggsave("PCA/PCA_batch.pdf",plot = PCAplot)

PCA_heatmap<-DimHeatmap(scRNA,dims = 1:15,cells = 500,balanced = TRUE)
#PC1:15的基因在500个细胞中表达值的展示
ggsave("PCA/PC1:15_heatmap.pdf",plot = PCA_heatmap)


###确定聚类使用的PC个数

scRNA<-JackStraw(scRNA,num.replicate = 100)

scRNA<-ScoreJackStraw(scRNA,dims=1:20)
JackStrawPlot(scRNA,dims=1:15)
ElbowPlot(scRNA)
#确定使用PC个数


# clustering --------------------------------------------------------------

scRNA<-FindNeighbors(scRNA,dims=1:20)
scRNA<-FindClusters(scRNA,resolution = 0.65)

head(Idents(scRNA),5)
#查看每个细胞对应的聚类数ID

head(scRNA@meta.data)
table(scRNA@meta.data$seurat_clusters)
#查看每个cluster中有多少细胞，cluster—0中的细胞数最多

scRNA<-RunTSNE(scRNA,dims=1:20)
scRNA<-RunUMAP(scRNA,dims = 1:20)

tsne_plot<-DimPlot(scRNA,reduction = "tsne",group.by = "orig.ident")
umap_plot<-DimPlot(scRNA,reduction = "umap",group.by = "orig.ident")

dir.create("Clustering")
ggsave("Clustering/tsne.pdf",plot = tsne_plot)
ggsave("Clustering/umap.pdf",plot = umap_plot)

tsne_plot2<-DimPlot(scRNA,reduction = "tsne",label = T,label.size=5)
tsne_plot3<-DimPlot(scRNA,reduction = "tsne",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/tsne2.pdf",plot = tsne_plot2)
ggsave("Clustering/tsne3.pdf",plot = tsne_plot3)

umap_plot2<-DimPlot(scRNA,reduction = "umap",label = T,label.size=5)
umap_plot3<-DimPlot(scRNA,reduction = "umap",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/umap2.pdf",plot = umap_plot2)
ggsave("Clustering/umap3.pdf",plot = umap_plot3)


# DEG analysis ---------------------------------------------------------

#识别差异基因
scRNA.pos.markers<-FindAllMarkers(scRNA,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#only.pos=TRUE表示只保留上调的差异基因

write.csv(scRNA.pos.markers,file = "e:/Rwork/GSE144945/Integrate/DEG analysis/Allmarkers-pos-C20.csv")
#该markers文件用于后边的细胞类型注释
head(scRNA.pos.markers)

Top_pos_genes<-scRNA.pos.markers%>%group_by(cluster)%>%top_n(n=1,wt=avg_log2FC)
#找出各个cluster中的TOP差异基因，按logFC值排序，n=2展示各个cluster的TOP2

# 展示候选marker基因的表达量 --------------------------------------------------------
#从上一部分“输出marker基因”输出的cell_type_markers.csv文件里选出每类p_val_adj排名第一的基因，例如第II类(class为1)的ENSG00000205364。

library("RColorBrewer")
display.brewer.all(type = "qual")

top<-Top_genes$gene%>%unique()

FeaturePlot(scRNA, reduction = "tsne",
            features = top[10:18], 
            ncol  = 3, #画在2列里
            cols = c("#66C2A5","grey","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            #do.return =T, #返回ggplot2对象，便于组图
            #no.legend = F, #显示图例
            pt.size = 1
) 
ggsave("pos_marker_C9-C20.pdf")

# 细胞类型手动注释Cell types annotation ------------------------------------------------------------------

#读入参考marker文件
ref_marker<-fread("E:/Rwork/GSE144945/Integrate/cell type annotation/RGEP_ref_marker.txt")%>%as.data.frame()

#cluster-specific up-markers
pbmc.markers<-fread("E:/Rwork/GSE144945/Integrate/cell type annotation/Allmarkers-pos-C20.csv")

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
                     celltype=paste(ref_marker$CellType[j], pbmc.markers$cluster[i],sep = "_"),
                     phyper_pvalue=p_value)
    p<-rbind(p,df_p)
    
  }
}

annotation_p<-subset(p,p$phyper_pvalue<0.5)

test<-p %>% group_by(cluster) %>% dplyr::mutate(min_rank(phyper_pvalue))
annotation_top<-subset(test,test$`min_rank(phyper_pvalue)`==1)


#####重新定义细胞类型######

new.cluster.ids<-annotation_top$celltype

names(new.cluster.ids)<-levels(scRNA)

Idents(scRNA)
#原来细胞聚类的名称

scRNA<-RenameIdents(scRNA,new.cluster.ids)
Idents(scRNA)
head(scRNA@meta.data)
#更改后的聚类名

DimPlot(scRNA,reduction = "tsne",label = TRUE,pt.size = 0.5)
DimPlot(scRNA,reduction = "umap",label = TRUE,pt.size = 0.5)#+NoLegend()
#可视化
saveRDS(scRNA,file = "E:/Rwork/GSE144945/Integrate/cell type annotation/scRNA_Annotation.rds")
