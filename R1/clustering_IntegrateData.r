library(Seurat)
library(SeuratData)
library(SeuratObject)
library(patchwork)
library(stringr)
library(ggplot2)
library(ggridges)


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

# 整合10个样本的单细胞数据 -----------------------------------------------------------

names(scRNAlist) <- samples_name
#为了整合，选择跨数据集重复可变的特征基因
features <- SelectIntegrationFeatures(object.list = scRNAlist)

#使用FindIntegrationAnchors()函数来识别anchors ，该函数将Seurat对象的列表作为输入，
immune.anchors <- FindIntegrationAnchors(object.list = scRNAlist,anchor.features = features)

saveRDS(immune.anchors,file = "immune.anchors.rds")


# 服务器-整合 ------------------------------------------------------------------


immune.anchors<-readRDS('E:/Rwork/GSE144945/Integrate_Anchor/immune.anchors.rds')

#并且IntegrateData()函数使用这些anchors将两个数据集整合在一起。
immune.combined <- IntegrateData(anchorset = immune.anchors)
DefaultAssay(immune.combined) <- "integrated"
 
saveRDS(immune.combined,file = "scRNAlist_anchor.rds")



###给列表命名，并保存数据
dir.create("Integrate_Anchor")
setwd("./Integrate_Anchor")

DefaultAssay(immune.combined) <- "integrated"


