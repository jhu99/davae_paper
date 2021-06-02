library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

brain <- LoadData("stxBrain", type = "anterior1")
InstallData("stxBrain")
InstallData("/Users/zhongyuanke/data/seurat_data/stxBrain.SeuratData_0.1.1.tar.gz",repo=NULL)
devtools::install_local("/Users/zhongyuanke/data/seurat_data/stxBrain.SeuratData_0.1.1.tar.gz")

install.packages('/Users/zhongyuanke/data/seurat_data/stxBrain.SeuratData_0.1.1.tar.gz', type='source')
load('/Users/zhongyuanke/data/seurat_data/stxBrain.SeuratData/data/anterior1.rda')

#--------ubuntu
library('Seurat')
allen_reference<-readRDS('./data/spatial/allen_cortex.rds')
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE)
save(allen_reference,file='./data/spatial/allen_trans.Rdata')

load(file='/Users/zhongyuanke/data/seurat_data/spatial/allen_trans.rds')
data<-allen_reference@assays[["SCT"]]@scale.data
data <- t(data)
row.names(data)<-allen_reference@assays[["SCT"]]@var.features

DimPlot(allen_reference, group.by = "subclass", label = TRUE)

allen_reference<-readRDS('/Users/zhongyuanke/data/seurat_data/spatial/allen_cortex.rds')
library(dplyr)
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(dims = 1:30)

