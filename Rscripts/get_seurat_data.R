library(Seurat)
library(SeuratData)
# ----------------------------panc8----------------------------------
data("panc8")
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq", "celseq2", "fluidigmc1", "smartseq2", "indrop")]

for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 13000, verbose = TRUE)
#  pancreas.list[[i]] <- ScaleData(pancreas.list[[i]], features =c(VariableFeatures(pancreas.list[[i]]) ))
  pancreas.list[[i]] <- subset(pancreas.list[[i]], features =c(VariableFeatures(pancreas.list[[i]]) ))
}

# panc8_celltype<-panc8@meta.data[["celltype"]]
# cell_type_path = '/Users/zhongyuanke/data/seurat_data/panc8_cell_type.csv'
# panc8_tech<-panc8@meta.data[["tech"]]
# tech_path = '/Users/zhongyuanke/data/seurat_data/panc8_tech.csv'
# a<-GetAssayData(panc8)
# data_path = '/Users/zhongyuanke/data/seurat_data/panc8.csv'
# write.csv(panc8_tech, tech_path)


celseq_scale <- pancreas.list[["celseq"]][["RNA"]]@scale.data
celseq_scale <- t(celseq_scale)
celseq_celltype<-pancreas.list[["celseq"]]@meta.data[["celltype"]]
celseq_data<-GetAssayData(pancreas.list[["celseq"]])
celseq_data <- as.matrix(celseq_data)
celseq_data <- t(celseq_data)
celltype_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq_cell_type.csv'
data_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq_data_norm_13000.csv'
scale_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq_data_scale.csv'
#write.csv(celseq_celltype, celltype_path)
write.csv(celseq_scale, scale_path)
write.csv(celseq_data, data_path)


celseq2_celltype<-pancreas.list[["celseq2"]]@meta.data[["celltype"]]
celseq2_scale <- pancreas.list[["celseq2"]][["RNA"]]@scale.data
celseq2_scale <- t(celseq2_scale)
celseq2_data<-GetAssayData(pancreas.list[["celseq2"]])
celseq2_data <- as.matrix(celseq2_data)
celseq2_data <- t(celseq2_data)
celltype_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq2_cell_type.csv'
data_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq2_data_norm_13000.csv'
scale_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/celseq2_data_scale.csv'
# write.csv(celseq2_celltype, celltype_path)
write.csv(celseq2_data, data_path)
write.csv(celseq2_scale, scale_path)


smartseq2_celltype<-pancreas.list[["smartseq2"]]@meta.data[["celltype"]]
smartseq2_scale <- pancreas.list[["smartseq2"]][["RNA"]]@scale.data
smartseq2_scale <- t(smartseq2_scale)
smartseq2_data<-GetAssayData(pancreas.list[["smartseq2"]])
smartseq2_data <- as.matrix(smartseq2_data)
smartseq2_data <- t(smartseq2_data)
celltype_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/smartseq2_cell_type.csv'
scale_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/smartseq2_data_scale.csv'
data_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/smartseq2_data_norm_13000.csv'
# write.csv(smartseq2_celltype, celltype_path)
write.csv(smartseq2_data, data_path)
write.csv(smartseq2_scale, scale_path)



fluidigmc1_celltype<-pancreas.list[["fluidigmc1"]]@meta.data[["celltype"]]
fluidigmc1_scale <- pancreas.list[["fluidigmc1"]][["RNA"]]@scale.data
fluidigmc1_scale <- t(fluidigmc1_scale)
fluidigmc1_data<-GetAssayData(pancreas.list[["fluidigmc1"]])
fluidigmc1_data <- as.matrix(fluidigmc1_data)
fluidigmc1_data <- t(fluidigmc1_data)
celltype_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/fluidigmc1_cell_type.csv'
scale_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/fluidigmc1_data_scale.csv'
data_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/fluidigmc1_data_norm_13000.csv'
write.csv(fluidigmc1_celltype, celltype_path)
write.csv(fluidigmc1_data, data_path)
write.csv(fluidigmc1_scale, scale_path)


indrop_celltype<-pancreas.list[["indrop"]]@meta.data[["celltype"]]
indrop_scale <- pancreas.list[["indrop"]][["RNA"]]@scale.data
indrop_scale <- t(indrop_scale)
indrop_data<-GetAssayData(pancreas.list[["indrop"]])
indrop_data <- as.matrix(indrop_data)
indrop_data <- t(indrop_data)
celltype_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/indrop_cell_type.csv'
scale_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/indrop_data_scale.csv'
data_path <- '/Users/zhongyuanke/data/seurat_data/panc_8/indrop_data_norm_13000.csv'
# write.csv(smartseq2_celltype, celltype_path)
write.csv(indrop_data, data_path)
write.csv(indrop_scale, scale_path)

# ----------------------------ifnb----------------------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
InstallData('ifnb')

data('ifnb')
ifnb<-UpdateSeuratObject(ifnb)
library(SeuratDisk)
SaveH5Seurat(ifnb, filename = "/Users/zhongyuanke/data/seurat_data/ifnb/ifnb.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_data/ifnb/ifnb.h5Seurat", dest = "h5ad")

ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 5000)
  x <- subset(x, features =c(VariableFeatures(x) ))
})
controll<-GetAssayData(ifnb.list[["CTRL"]])
controll <- as.matrix(controll)
controll <- t(controll)

stim<- GetAssayData(ifnb.list[["STIM"]])
stim <- as.matrix(stim)
stim <- t(stim)

ctrl_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/ctrl_5000.csv'
stim_path <- '/Users/zhongyuanke/data/seurat_data/ifnb/stim_5000.csv'
write.csv(controll, ctrl_path)
write.csv(stim, stim_path)


# ----------------------------panc8----------------------------------
InstallData('panc8')
data('panc8')


