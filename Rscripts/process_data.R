library(Seurat)
library(MixGHD)
library(ggplot2)
adata<-ReadH5AD(file = '/Users/zhongyuanke/data/pbmc/zheng/293t_jurkat_merge.h5ad')
adata <- NormalizeData(adata, verbose = FALSE)
adata <- FindVariableFeatures(adata, selection.method = "vst", 
                                           nfeatures = 3000, verbose = TRUE)
adata <- subset(adata, features =c(VariableFeatures(adata) ))

matrix <- GetAssayData(adata)
matrix <- as.matrix(matrix)
matrix <- t(matrix)

write.csv(matrix, '/Users/zhongyuanke/data/pbmc/zheng/293t_jurkat_merge_filt.csv')

#--------------------HCA---------------------------------------
adata <- Read10X_h5('/Users/zhongyuanke/data/HCA/ica_cord_blood_h5.h5')
adata<-CreateSeuratObject(adata, assay='RNA')
adata <- NormalizeData(adata, verbose = FALSE)
adata <- FindVariableFeatures(adata, selection.method = "vst", 
                              nfeatures = 7000, verbose = TRUE)
adata <- subset(adata, features =c(VariableFeatures(adata) ))

matrix <- GetAssayData(adata)

writeMM(matrix)
matrix <- t(matrix)
write.csv(matrix, '/Users/zhongyuanke/data/pbmc/zheng/ica_cord_blood7000.csv')

#--------------------kbet--293t-------------------------------------
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(MixGHD)

adata<-ReadH5AD(file = '/Users/zhongyuanke/data/dann_vae/k_bet/293t/293t_scan.h5ad')
batch <-adata@meta.data[["batch"]]
cluster <-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/orig_label.csv')

davae_batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/davae_batch.csv')
davae_batch <- t(davae_batch)
davae_batch<-split(davae_batch,cluster)
batch0 <- davae_batch[['0']]
batch1 <- davae_batch[['1']]
# davae_label<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/davae_label.csv')
# davae_label <- t(davae_label)
# 
# orig_label<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/orig_label.csv')
# orig_label <- t(orig_label)

data_orig<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/orig.csv')
data_desc<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/desc.csv')
data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/davae.csv')
data_scan<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/scan.csv')
data_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/seurat.csv')
data_scgen<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/scgen.csv')

batch<-split(batch,cluster)
batch0 <- batch[['0']]
batch1 <- batch[['1']]

scgen_batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/scgen_batch.csv')
scgen_label<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/293t/scgen_label.csv')
scgen_batch <-t(scgen_batch)
scgen_batch <- as.integer(scgen_batch)
batch<-split(scgen_batch,scgen_label)
batch0 <- batch[['0']]
batch1 <- batch[['1']]

data_orig <- umap(data_orig, n_components=6)
data_orig <- data_orig[["layout"]]
 
data_desc <- umap(data_desc, n_components=6)
data_desc <- data_desc[["layout"]]

data_scan <- umap(data_scan, n_components=6)
data_scan <- data_scan[["layout"]]

data_seurat <- umap(data_seurat, n_components=6)
data_seurat <- data_seurat[["layout"]]

data_davae<-split(data_davae,cluster) 
data_desc<-split(data_desc,cluster) 
data_orig<-split(data_orig,cluster) 
data_scan<-split(data_scan,cluster) 
data_seurat<-split(data_seurat,cluster) 
data_scgen<-split(data_scgen,scgen_label)

data_1<-data_davae[["0"]]
data_2<-data_davae[["1"]]

data_1<-data_orig[["0"]]
data_2<-data_orig[["1"]]

data_1<-data_desc[["0"]]
data_2<-data_desc[["1"]]

data_1<-data_scan[["0"]]
data_2<-data_scan[["1"]]

data_1<-data_seurat[["0"]]
data_2<-data_seurat[["1"]]

data_1<-data_scgen[["0"]]
data_2<-data_scgen[["1"]]

subset_size <- 0.05 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch0), size = floor(subset_size * length(batch0)), replace=FALSE)
batch.estimate_1 <- kBET(data_1[subset_id,], batch0[subset_id])

subset_id <- sample.int(n = length(batch1), size = floor(subset_size * length(batch1)), replace=FALSE)
batch.estimate_2 <- kBET(data_2[subset_id,], batch1[subset_id])

signif(mean(batch.estimate_1[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_2[["summary"]][["kBET.observed"]], 3))

subset_size <- 0.2 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
batch.estimate <- kBET(data_orig[subset_id,], batch[subset_id])

write.csv(result, '/Users/zhongyuanke/data/dann_vae/k_bet/plot/293t_davae_1.csv', 
          row.names=FALSE
)
# ----------------kbet smartseq-------------------------------------
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(MixGHD)

cluster <-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/label.csv')
batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/batch.csv')
batch <-t(batch)
batch <- as.integer(batch)
batch<-split(batch,cluster)

data_orig<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/orig.csv')
data_desc<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/desc.csv')
data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/davae.csv')
data_scan<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/scan.csv')
data_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/seurat.csv')
data_scgen<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/smartseq/scgen.csv')

batch0 <- batch[['0']]
batch1 <- batch[['1']]
batch2 <- batch[['2']]
batch3 <- batch[['3']]

data_orig<-split(data_orig,cluster) 
data_davae<-split(data_davae,cluster) 
data_desc<-split(data_desc,cluster) 
data_scan<-split(data_scan,cluster) 
data_seurat<-split(data_seurat,cluster) 
data_scgen<-split(data_scgen,cluster) 

data_0<-data_orig[["0"]]
data_1<-data_orig[["1"]]
data_2<-data_orig[["2"]]
data_3<-data_orig[["3"]]
 
data_0<-data_davae[["0"]]
data_1<-data_davae[["1"]]
data_2<-data_davae[["2"]]
data_3<-data_davae[["3"]]


data_0<-data_desc[["0"]]
data_1<-data_desc[["1"]]
data_2<-data_desc[["2"]]
data_3<-data_desc[["3"]]

data_0<-data_scan[["0"]]
data_1<-data_scan[["1"]]
data_2<-data_scan[["2"]]
data_3<-data_scan[["3"]]

data_0<-data_seurat[["0"]]
data_1<-data_seurat[["1"]]
data_2<-data_seurat[["2"]]
data_3<-data_seurat[["3"]]

data_0<-data_scgen[["0"]]
data_1<-data_scgen[["1"]]
data_2<-data_scgen[["2"]]
data_3<-data_scgen[["3"]]
require('FNN')
knn <- get.knn(data_0, k=2, algorithm = 'cover_tree')
subset_size <- 0.05 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch0), size = floor(subset_size * length(batch0)), replace=FALSE)
batch.estimate_0 <- kBET(data_0[subset_id,], batch0[subset_id],k=2,knn=knn)

subset_id <- sample.int(n = length(batch1), size = floor(subset_size * length(batch1)), replace=FALSE)
batch.estimate_1 <- kBET(data_1[subset_id,], batch1[subset_id],k=2,knn=knn)

subset_id <- sample.int(n = length(batch2), size = floor(subset_size * length(batch2)), replace=FALSE)
batch.estimate_2 <- kBET(data_2[subset_id,], batch2[subset_id],k=2,knn=knn)

subset_id <- sample.int(n = length(batch3), size = floor(subset_size * length(batch3)), replace=FALSE)
batch.estimate_3 <- kBET(data_3[subset_id,], batch3[subset_id],k=2,knn=knn)

signif(mean(batch.estimate_0[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_1[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_2[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_3[["summary"]][["kBET.observed"]], 3))

write.csv(result, '/Users/zhongyuanke/data/dann_vae/k_bet/plot/293t_davae_1.csv', 
          row.names=FALSE
)


# -------------------------------kbet ifnb---------------------------------------------
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(MixGHD)

cluster <-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/label.csv')
batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/batch.csv')

batch <-t(batch)
batch <- as.integer(batch)

batch<-split(batch,cluster)

data_orig<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/orig.csv')
data_desc<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/desc.csv')
data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/davae.csv')
data_scan<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/scan.csv')
data_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/seurat.csv')
data_scgen <- read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/ifnb/scgen.csv')

batch2 <- batch[['2']]
batch3 <- batch[['3']]
batch4 <- batch[['4']]
batch5 <- batch[['5']]

data_davae<-split(data_davae,cluster) 
data_desc<-split(data_desc,cluster) 
data_orig<-split(data_orig,cluster) 
data_scan<-split(data_scan,cluster) 
data_seurat<-split(data_seurat,cluster) 
data_scgen<-split(data_scgen,cluster)

data_2<-data_davae[["2"]]
data_3<-data_davae[["3"]]
data_4<-data_davae[["4"]]
data_5<-data_davae[["5"]]

data_2<-data_orig[["2"]]
data_3<-data_orig[["3"]]
data_4<-data_orig[["4"]]
data_5<-data_orig[["5"]]

data_2<-data_desc[["2"]]
data_3<-data_desc[["3"]]
data_4<-data_desc[["4"]]
data_5<-data_desc[["5"]]

data_2<-data_scan[["2"]]
data_3<-data_scan[["3"]]
data_4<-data_scan[["4"]]
data_5<-data_scan[["5"]]

data_2<-data_seurat[["2"]]
data_3<-data_seurat[["3"]]
data_4<-data_seurat[["4"]]
data_5<-data_seurat[["5"]]

data_2<-data_scgen[["2"]]
data_3<-data_scgen[["3"]]
data_4<-data_scgen[["4"]]
data_5<-data_scgen[["5"]]

names(batch0)

subset_size <- 0.05#subsample to 10% of the data
subset_id <- sample.int(n = length(batch2), size = floor(subset_size * length(batch2)), replace=FALSE)
batch.estimate <- kBET(data_2[subset_id,], batch2[subset_id])
signif(mean(batch.estimate[["summary"]][["kBET.observed"]], 3))

subset_size <- 0.2
subset_id <- sample.int(n = length(batch2), size = floor(subset_size * length(batch2)), replace=FALSE)
batch.estimate_0 <- kBET(data_2[subset_id,], batch2[subset_id])

subset_id <- sample.int(n = length(batch3), size = floor(subset_size * length(batch3)), replace=FALSE)
batch.estimate_1 <- kBET(data_3[subset_id,], batch3[subset_id])

subset_id <- sample.int(n = length(batch4), size = floor(subset_size * length(batch4)), replace=FALSE)
batch.estimate_2 <- kBET(data_4[subset_id,], batch4[subset_id])

subset_id <- sample.int(n = length(batch5), size = floor(subset_size * length(batch5)), replace=FALSE)
batch.estimate_3 <- kBET(data_5[subset_id,], batch5[subset_id])

signif(mean(batch.estimate_0[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_1[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_2[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_3[["summary"]][["kBET.observed"]], 3))

# -------------------------------kbet multimodal---------------------------------------------
gc()
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(MixGHD)

cluster <-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/label.csv')
batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/batch.csv')

batch <-t(batch)
batch<-split(batch,cluster)

data_orig<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/orig.csv')
data_desc<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/desc.csv')
data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/davae.csv')
data_scan<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/scan.csv')
data_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/seurat.csv')
data_scgen <- read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/scgen.csv')

batch2 <- batch[['5']]
batch3 <- batch[['2']]
batch4 <- batch[['3']]
batch5 <- batch[['1']]

data_davae<-split(data_davae,cluster) 
data_desc<-split(data_desc,cluster) 
data_orig<-split(data_orig,cluster) 
data_scan<-split(data_scan,cluster) 
data_seurat<-split(data_seurat,cluster) 
data_scgen<-split(data_scgen,cluster)

data_2<-data_davae[["5"]]
data_3<-data_davae[["2"]]
data_4<-data_davae[["3"]]
data_5<-data_davae[["1"]]

data_2<-data_orig[["5"]]
data_3<-data_orig[["2"]]
data_4<-data_orig[["3"]]
data_5<-data_orig[["1"]]

data_2<-data_desc[["5"]]
data_3<-data_desc[["2"]]
data_4<-data_desc[["3"]]
data_5<-data_desc[["1"]]

data_2<-data_scan[["5"]]
data_3<-data_scan[["2"]]
data_4<-data_scan[["3"]]
data_5<-data_scan[["1"]]

data_2<-data_seurat[["5"]]
data_3<-data_seurat[["2"]]
data_4<-data_seurat[["3"]]
data_5<-data_seurat[["1"]]

data_2<-data_scgen[["5"]]
data_3<-data_scgen[["2"]]
data_4<-data_scgen[["3"]]
data_5<-data_scgen[["1"]]

names(batch0)

subset_size <- 0.25
subset_id <- sample.int(n = length(batch2), size = floor(subset_size * length(batch2)), replace=FALSE)
batch.estimate_0 <- kBET(data_2[subset_id,], batch2[subset_id])

subset_id <- sample.int(n = length(batch3), size = floor(subset_size * length(batch3)), replace=FALSE)
batch.estimate_1 <- kBET(data_3[subset_id,], batch3[subset_id])

subset_id <- sample.int(n = length(batch4), size = floor(subset_size * length(batch4)), replace=FALSE)
batch.estimate_2 <- kBET(data_4[subset_id,], batch4[subset_id])

subset_id <- sample.int(n = length(batch5), size = floor(subset_size * length(batch5)), replace=FALSE)
batch.estimate_3 <- kBET(data_5[subset_id,], batch5[subset_id])

signif(mean(batch.estimate_0[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_1[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_2[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_3[["summary"]][["kBET.observed"]], 3))




# --------------------------------kbet mouse -------------------------------
gc()
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(MixGHD)

cluster <-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/label.csv')
batch<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/batch.csv')

batch <-t(batch)
batch<-split(batch,cluster)

data_orig<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/orig.csv')
data_desc<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/desc.csv')
data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/davae.csv')
data_scan<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/scan.csv')
data_scgen <- read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/mouse/scgen.csv')

batch2 <- batch[['8']]
batch3 <- batch[['9']]
batch4 <- batch[['0']]
batch5 <- batch[['1']]

data_davae<-split(data_davae,cluster) 
data_desc<-split(data_desc,cluster) 
data_orig<-split(data_orig,cluster) 
data_scan<-split(data_scan,cluster) 
data_seurat<-split(data_seurat,cluster) 
data_scgen<-split(data_scgen,cluster)

data_2<-data_davae[["8"]]
data_3<-data_davae[["9"]]
data_4<-data_davae[["0"]]
data_5<-data_davae[["1"]]

data_2<-data_orig[["8"]]
data_3<-data_orig[["9"]]
data_4<-data_orig[["0"]]
data_5<-data_orig[["1"]]

data_2<-data_desc[["8"]]
data_3<-data_desc[["9"]]
data_4<-data_desc[["0"]]
data_5<-data_desc[["1"]]

data_2<-data_scan[["8"]]
data_3<-data_scan[["9"]]
data_4<-data_scan[["0"]]
data_5<-data_scan[["1"]]

data_2<-data_seurat[["8"]]
data_3<-data_seurat[["9"]]
data_4<-data_seurat[["0"]]
data_5<-data_seurat[["1"]]

data_2<-data_scgen[["8"]]
data_3<-data_scgen[["9"]]
data_4<-data_scgen[["0"]]
data_5<-data_scgen[["1"]]

names(batch0)

subset_size <- 0.05
subset_id <- sample.int(n = length(batch2), size = floor(subset_size * length(batch2)), replace=FALSE)
batch.estimate_0 <- kBET(data_2[subset_id,], batch2[subset_id])

subset_id <- sample.int(n = length(batch3), size = floor(subset_size * length(batch3)), replace=FALSE)
batch.estimate_1 <- kBET(data_3[subset_id,], batch3[subset_id])

subset_id <- sample.int(n = length(batch4), size = floor(subset_size * length(batch4)), replace=FALSE)
batch.estimate_2 <- kBET(data_4[subset_id,], batch4[subset_id])

subset_id <- sample.int(n = length(batch5), size = floor(subset_size * length(batch5)), replace=FALSE)
batch.estimate_3 <- kBET(data_5[subset_id,], batch5[subset_id])

signif(mean(batch.estimate_0[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_1[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_2[["summary"]][["kBET.observed"]], 3))
signif(mean(batch.estimate_3[["summary"]][["kBET.observed"]], 3))

# --------------------------------kbet atac -------------------------------
library(kBET)
library(ggplot2)
library(umap)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(MixGHD)

Convert("/Users/zhongyuanke/data/dann_vae/atac/davae_02.h5ad",
        dest = "h5seurat", overwrite = TRUE)
atac_davae <- LoadH5Seurat("/Users/zhongyuanke/data/dann_vae/atac/davae_save03.h5seurat")

data_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/atac/davae.csv')
data_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/atac/seurat.csv')
batch_davae<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/atac/batch_davae.csv')
batch_davae <-t(batch_davae)
batch_seurat<-read.csv('/Users/zhongyuanke/data/dann_vae/k_bet/atac/batch_seurat.csv')
batch_seurat <-t(batch_seurat)

subset_size <- 0.1#subsample to 10% of the data
subset_id <- sample.int(n = length(batch_davae), size = floor(subset_size * length(batch_davae)), replace=FALSE)
batch.estimate <- kBET(data_davae[subset_id,], batch_davae[subset_id])
signif(mean(batch.estimate[["summary"]][["kBET.observed"]]), 3)

subset_id <- sample.int(n = length(batch_seurat), size = floor(subset_size * length(batch_seurat)), replace=FALSE)
batch.estimate <- kBET(data_seurat[subset_id,], batch_seurat[subset_id])
signif(mean(batch.estimate[["summary"]][["kBET.observed"]]), 3)

subset_size <- 0.1#subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
batch.estimate <- kBET(data[subset_id,], batch[subset_id])
signif(batch.estimate[["summary"]][["kBET.observed"]], 3)
signif(mean(batch.estimate[["summary"]][["kBET.observed"]]), 3)

library(umap)
integrated<-LoadH5Seurat("/Users/zhongyuanke/data/seurat_result/atac.h5Seurat")

# -------plot-------------------

SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/atac_umap.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/atac_umap.h5Seurat", dest = "h5ad")

inte_pca<-as.matrix(integrated@reductions[["pca"]]@cell.embeddings)
batch = integrated@meta.data[["orig.ident"]]
# -----------------plot------------------------
p1<-DimPlot(atac, group.by = 'batch')
p2<-DimPlot(atac, group.by = 'cell.type')
p1+p2

atac@reductions[["seurat_umap"]]<-integrated@reductions[["umap"]]

data <- atac@reductions[["davae"]]@cell.embeddings
batch <-atac@meta.data[["batch"]]







