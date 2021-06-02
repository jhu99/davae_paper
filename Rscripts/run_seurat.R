library(Seurat)
library(SeuratData)
data("panc8")
# -----------------save-----------------------------
pancreas.list <- SplitObject(panc8, split.by = "replicate")
library(SeuratDisk)
panc8<-UpdateSeuratObject(panc8)
SaveH5Seurat(panc8, filename = "/Users/zhongyuanke/data/seurat_data/panc_8/panc8_8.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_data/panc_8/panc8_8.h5Seurat", dest = "h5ad")

#-------------------------------------
pancreas.list <- SplitObject(panc8, split.by = "tech")
pancreas.list <- pancreas.list[c("celseq2", "celseq",'fluidigmc1','indrop','smartseq2')]
for (i in 1:length(pancreas.list)) {
  pancreas.list[[i]] <- NormalizeData(pancreas.list[[i]], verbose = FALSE)
  pancreas.list[[i]] <- FindVariableFeatures(pancreas.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}

reference.list <- pancreas.list[c("celseq2", "celseq", 'fluidigmc1','indrop','smartseq2')] 
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

library(SeuratDisk)
DefaultAssay(pancreas.integrated) <- "integrated"
SaveH5Seurat(pancreas.integrated, filename = "/Users/zhongyuanke/data/seurat_result/panc8_all.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/panc8_all.h5Seurat", dest = "h5ad")

integrated_matrix<-GetAssayData(pancreas.integrated)
integrated_matrix <- as.matrix(integrated_matrix)
integrated_matrix <- t(integrated_matrix)
data_path <- '/Users/zhongyuanke/data/seurat_result/pan8_seurat.csv'
write.csv(integrated_matrix, data_path)

# -------------------pbmc3k----------------------------
library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30)
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "tech")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "celltype", label = TRUE, 
              repel = TRUE) + NoLegend()


install.packages("/Users/zhongyuanke/data/seurat_data/pbmc3k.SeuratData_3.1.4.tar.gz", repos = NULL, type = "source")


# -------------------------mcl---------------------
gc()
file1 <- '/Users/zhongyuanke/data/pbmc/zheng/293t.h5ad'
file2 <- '/Users/zhongyuanke/data/pbmc/zheng/jurkat.h5ad'
file3 <- '/Users/zhongyuanke/data/pbmc/zheng/293t_jurkat_50_50.h5ad'
library(SeuratDisk)
adata1 <- ReadH5AD(file1)
adata2 <- ReadH5AD(file2)
adata3 <- ReadH5AD(file3)

Convert(file3, dest = "h5seurat", overwrite = TRUE)
adata3<-LoadH5Seurat('/Users/zhongyuanke/data/pbmc/zheng/293t_jurkat_50_50.h5seurat')

adata1 <- NormalizeData(adata1, verbose = FALSE)
adata1 <- FindVariableFeatures(adata1, selection.method = "vst", 
                              nfeatures = 2000, verbose = TRUE)

adata2 <- NormalizeData(adata2, verbose = FALSE)
adata2 <- FindVariableFeatures(adata2, selection.method = "vst", 
                               nfeatures = 2000, verbose = TRUE)
adata3 <- NormalizeData(adata3, verbose = FALSE)
adata3 <- FindVariableFeatures(adata3, selection.method = "vst", 
                               nfeatures = 2000, verbose = TRUE)

RenameCells(adata1, new.names = adata1@meta.data[["barcodes"]])
RenameCells(adata2, new.names = adata2@meta.data[["barcodes"]])
RenameCells(adata3, new.names = adata3@meta.data[["barcodes"]])

reference.list <- c(pbmc293t=adata1,jurkat=adata2,mix=adata3)
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

integrated_matrix<-GetAssayData(integrated)
integrated_matrix <- as.matrix(integrated_matrix)
integrated_matrix <- t(integrated_matrix)
data_path <- '/Users/zhongyuanke/data/seurat_result/293t.csv'
write.csv(integrated_matrix, data_path)

SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/mcl.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/mcl.h5Seurat", dest = "h5ad")

# --------------ifnb---------------------------

library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)

data("ifnb")
ifnb.list <- SplitObject(ifnb, split.by = "stim")

ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

integrated_matrix<-GetAssayData(immune.combined)
integrated_matrix <- as.matrix(integrated_matrix)
integrated_matrix <- t(integrated_matrix)
data_path <- '/Users/zhongyuanke/data/seurat_result/ifnb.csv'
write.csv(integrated_matrix, data_path)

immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
SaveH5Seurat(immune.combined, filename = "/Users/zhongyuanke/data/seurat_result/ifnb.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/ifnb.h5Seurat", dest = "h5ad")

# -----------------dc------------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
help(package='Seurat')
file1 <- '/Users/zhongyuanke/data/dann_vae/benchmark1/batch1.csv'
file2 <- '/Users/zhongyuanke/data/dann_vae/benchmark1/batch2.csv'

adata1 <- read.csv(file1,header=TRUE,row.names=1)
adata2 <- read.csv(file2,header=TRUE,row.names=1)

data1 <- CreateSeuratObject(counts = adata1)
data2 <- CreateSeuratObject(counts = adata2)

reference.list <- c(data1, data2)
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

integrated_matrix<-GetAssayData(integrated)
integrated_matrix <- as.matrix(integrated_matrix)
integrated_matrix <- t(integrated_matrix)
data_path <- '/Users/zhongyuanke/data/seurat_result/smartseq.csv'
write.csv(integrated_matrix, data_path)


# ------------------------------atac------------------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
library(SeuratDisk)
# help(package='Seurat')
file_atac <- '/Users/zhongyuanke/data/seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix_8000.csv'
file_rna <- '/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3_filtered_feature_bc_matrix.h5ad'
Convert(file_rna,
        dest = "h5seurat", overwrite = TRUE)
data1 <- LoadH5Seurat("/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3_filtered_feature_bc_matrix.h5seurat")
adata2 <-(pbmc@assays[["RNA"]]@data)
data2 <- CreateSeuratObject(counts = adata2)
reference.list <- c(data1, data2)
reference.list <- lapply(X = reference.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

pbmc_rna <- readRDS('/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3.rds')
DefaultAssay(pbmc) <- 'RNA'
anchors <- FindIntegrationAnchors(
  object.list = c(pbmc_rna, pbmc)
)

features <- SelectIntegrationFeatures(object.list = reference.list)
anchors <- FindIntegrationAnchors(object.list = reference.list, anchor.features=features)

integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
integrated <- AddMetaData(object = integrated, metadata = pbmc@active.ident, col.name = 'celltype')
SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/atac.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/atac.h5Seurat", dest = "h5ad")

integrated_matrix<-GetAssayData(integrated)
integrated_matrix <- as.matrix(integrated_matrix)
integrated_matrix <- t(integrated_matrix)
data_path <- '/Users/zhongyuanke/data/seurat_result/atac.csv'
write.csv(integrated_matrix, data_path)


# -------------------------spatial-------------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
help(package='Seurat')

file_Anterior_cortex <- '/Users/zhongyuanke/data/spatial/mouse_brain/10x_mouse_brain_Anterior/Anterior_cortex.h5ad'
file_Posterior_cortex <- '/Users/zhongyuanke/data/spatial/mouse_brain/10x_mouse_brain_Posterior/Posterior_cortex.h5ad'
library(SeuratDisk)
Convert(file_Anterior_cortex, dest = "h5seurat", overwrite = TRUE)
Convert(file_Posterior_cortex, dest = "h5seurat", overwrite = TRUE)

Anterior <- LoadH5Seurat("/Users/zhongyuanke/data/spatial/mouse_brain/10x_mouse_brain_Anterior/Anterior_cortex.h5seurat")
Posterior<- LoadH5Seurat("/Users/zhongyuanke/data/spatial/mouse_brain/10x_mouse_brain_Posterior/Posterior_cortex.h5seurat")

file_rna<-'/Users/zhongyuanke/data/spatial/mouse_brain/cortex_for_seurat.h5ad'
Convert(file_rna, dest = "h5seurat", overwrite = TRUE)
rna<-LoadH5Seurat('/Users/zhongyuanke/data/spatial/mouse_brain/cortex_for_seurat.h5seurat')


reference.list <- c(Posterior, rna)
reference.list <- lapply(X = reference.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

rm(rna)
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)


integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 10, verbose = FALSE)

SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/spatial/posterior_rna.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/spatial/posterior_rna.h5Seurat", dest = "h5ad")


# ----------------multimodal--------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(patchwork)
library(SeuratDisk)

rna<-LoadH5Seurat('/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5Seurat')
atac<-LoadH5Seurat('/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/activaty_matrix.h5Seurat')

reference.list <- c(atac, rna)
reference.list <- lapply(X = reference.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x <- subset(x, features =c(VariableFeatures(x) ))
})

rna<-NormalizeData(rna)
rna<-FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna<-subset(rna, features =c(VariableFeatures(rna) ))
SaveH5Seurat(rna, filename = "/Users/zhongyuanke/Desktop/temp_multimodal/rna.h5seurat")
Convert("/Users/zhongyuanke/Desktop/temp_multimodal/rna.h5seurat", dest = "h5ad")

atac<-NormalizeData(atac)
atac<-FindVariableFeatures(atac, selection.method = "vst", nfeatures = 2000)
atac<-subset(atac, features =c(VariableFeatures(atac)))
SaveH5Seurat(atac, filename = "/Users/zhongyuanke/Desktop/temp_multimodal/atac.h5seurat")
Convert("/Users/zhongyuanke/Desktop/temp_multimodal/atac.h5seurat", dest = "h5ad")




rm(atac)
gc()
anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
integrated <- IntegrateData(anchorset = anchors, dims = 1:50)

integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)

SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/multimodal.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/multimodal.h5Seurat", dest = "h5ad")
