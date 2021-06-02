# ------------------------ATAC----------------------------------------
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)

inputdata.10x <- Read10X_h5("/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"

frag.file <- "/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)
DefaultAssay(pbmc) <- "ATAC"
gene.activities <- GeneActivity(pbmc)
pbmc[['ACTIVATY']] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) 
pbmc <- RunPCA(pbmc, npcs=500)
pbmc <- RunUMAP(pbmc, dims = 1:500, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc, n=500)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:500, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DefaultAssay(pbmc) <- "ACTIVATY"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc, n=500)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:500, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# -------------------wnn--------------------------------------
pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:500, 2:500))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
pbmc <- FindSubCluster(pbmc, cluster = 6, graph.name = "wsnn", algorithm = 3)
Idents(pbmc) <- "sub.cluster"
pbmc <- RenameIdents(pbmc, '19' = 'pDC','20' = 'HSPC','15' = 'cDC')
pbmc <- RenameIdents(pbmc, '0' = 'CD14 Mono', '8' ='CD14 Mono', '5' = 'CD16 Mono')
pbmc <- RenameIdents(pbmc, '17' = 'Naive B', '11' = 'Intermediate B', '10' = 'Memory B', '21' = 'Plasma')
pbmc <- RenameIdents(pbmc, '7' = 'NK')
pbmc <- RenameIdents(pbmc, '4' = 'CD4 TCM', '13'= "CD4 TEM", '3' = "CD4 TCM", '16' ="Treg", '1' ="CD4 Naive", '14' = "CD4 Naive")
pbmc <- RenameIdents(pbmc, '2' = 'CD8 Naive', '9'= "CD8 Naive", '12' = 'CD8 TEM_1', '6_0' = 'CD8 TEM_2', '6_1' ='CD8 TEM_2')
pbmc <- RenameIdents(pbmc, '18' = 'MAIT')
pbmc <- RenameIdents(pbmc, '6_2' ='gdT', '6_3' = 'gdT')
pbmc$celltype <- Idents(pbmc)

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p3
# ---------------------------- get data ----------------------------------
celltype <- pbmc@meta.data[["celltype"]]
celltype_path <- '/Users/zhongyuanke/data/multimodal/pbmc_10k/celltype_filt.csv'
write.csv(celltype, celltype_path)

rna <- GetAssayData(pbmc@assays[["RNA"]])
rna <- as.matrix(rna)
rna <- t(rna)
rna_path <- '/Users/zhongyuanke/data/multimodal/pbmc_10k/rna_filt.csv'
write.csv(rna, rna_path)


gene.activities <- GeneActivity(pbmc)
gene_active <- GetAssayData(pbmc@assays[["SCT"]])
gene.activities <- as.matrix(gene.activities)
gene.activities <- t(gene.activities)
atac_path <- '/Users/zhongyuanke/data/multimodal/pbmc_10k/atac_filt.csv'
write.csv(gene.activities, atac_path)

lsi <- as.matrix(pbmc@reductions[["lsi"]]@cell.embeddings)
lsi_path <- '/Users/zhongyuanke/data/multimodal/pbmc_10k/lsi.csv'
write.csv(lsi, lsi_path)

atac <-as.matrix(pbmc@reductions[["lsi"]]@cell.embeddings)
atac_path <- '/Users/zhongyuanke/data/multimodal/pbmc_10k/atac_filt.csv'
write.csv(atac, atac_path)

save(pbmc,file='/Users/zhongyuanke/data/multimodal/pbmc_10k_filt.Rdata')
load(file='/Users/zhongyuanke/data/multimodal/pbmc_10k_filt.Rdata')
load(file='/Users/zhongyuanke/data/multimodal/pbmc.Rdata')


# ------------------------atac_h5ad_getdata-----------------------------
library(SeuratDisk)
DefaultAssay(pbmc) <- "RNA"
SaveH5Seurat(pbmc, filename = "/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5Seurat")
Convert("/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/rna.h5Seurat", dest = "h5ad")

DefaultAssay(pbmc) <- "ATAC"
SaveH5Seurat(pbmc, filename = "/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/atac.h5Seurat")
Convert("/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/atac.h5Seurat", dest = "h5ad")

DefaultAssay(pbmc) <- "ACTIVATY"
SaveH5Seurat(pbmc, filename = "/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/activaty_matrix.h5Seurat")
Convert("/Users/zhongyuanke/data/multimodal/atac_pbmc_10k/activaty_matrix.h5Seurat", dest = "h5ad")

# ------------------------Protein----------------------------------------
library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)
bm <- LoadData(ds = "bmcite")

# write h5ad
library(SeuratDisk)
DefaultAssay(bm) <- 'RNA'
SaveH5Seurat(bm, filename = "/Users/zhongyuanke/data/multimodal/bmcite/bmcite_rna.h5Seurat")
Convert("/Users/zhongyuanke/data/multimodal/bmcite/bmcite_rna.h5Seurat", dest = "h5ad")


DefaultAssay(bm) <- 'RNA'
bm <- NormalizeData(bm) %>% FindVariableFeatures(nfeatures=5000) %>% ScaleData() %>% RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

bm <- FindMultiModalNeighbors(
  bm, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
)

# ------------------------get data --------------------------------------
rna <- bm@assays[["RNA"]]@scale.data 
rna <- t(rna)
rna_path <- '/Users/zhongyuanke/data/multimodal/bmcite/rna_5k.csv'
write.csv(rna, rna_path)

protein <- bm@assays[["ADT"]]@scale.data
protein <- t(protein)
protein_path <- '/Users/zhongyuanke/data/multimodal/bmcite/protein.csv'
write.csv(protein, protein_path)

celltype <- bm@meta.data[["celltype.l2"]]
celltype_path <- '/Users/zhongyuanke/data/multimodal/bmcite/celltype_l2.csv'
write.csv(celltype, celltype_path)

# -------------------------- 10x protein -----------------------------
inputdata.10x <- Read10X_h5("/Users/zhongyuanke/data/multimodal/protein/vdj_v1_hs_pbmc2_5gex_protein_filtered_feature_bc_matrix.h5")

rna <- as.matrix(inputdata.10x[["Gene Expression"]])
rna <- t(rna)
rna_path <- '/Users/zhongyuanke/data/multimodal/protein/rna.csv'
write.csv(rna, rna_path)

protein <- as.matrix(inputdata.10x[["Antibody Capture"]])
protein <- t(protein)
protein_path <- '/Users/zhongyuanke/data/multimodal/protein/protein.csv'
write.csv(protein,protein_path)


# ---------------------------SHARE-seq----------------------------
