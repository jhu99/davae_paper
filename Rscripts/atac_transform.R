library(Seurat)
library(ggplot2)
library(patchwork)
base_path <- '/Users/zhongyuanke/data/atac/'
atac_path <- paste(base_path, 'p0_BrainCortex/chromatin/',
                   sep = "")
gtf_path <- paste(base_path, 'p0_BrainCortex/chromatin/Mus_musculus.GRCm38.100.gtf',
                  sep = "")
dna_path <- paste(base_path, 'p0_BrainCortex/dna/',
                  sep = "")
barcodes_path <- paste(base_path, 'p0_BrainCortex/chromatin/barcodes.csv',
                       sep = "")
chro<-Read10X(atac_path, gene.column = 1)
geneActMat<-CreateGeneActivityMatrix(peak.matrix = chro,
                                     annotation.file = gtf_path,
								                      seq.levels = c(1:19,"X","Y"))
barcodes = colnames(geneActMat)
write.csv(barcodes, barcodes_path)
chro <- CreateSeuratObject(counts = geneActMat, assay = "ATAC", project = "10x_ATAC")

chro <- FindVariableFeatures(chro, selection.method = "vst", 
                             nfeatures = 5000, verbose = TRUE)
chro <- subset(chro, features =c(VariableFeatures(chro) ))
chro <- NormalizeData(chro)
chro <- ScaleData(chro)

atac_norm <- GetAssayData(chro)
atac_norm <- as.matrix(atac_norm)
atac_norm <- t(atac_norm)
norm_write_path <- paste(base_path, 'p0_BrainCortex/chromatin/orig_norm_count_matrix.csv',
                    sep = "")
write.csv(atac_norm, norm_write_path)

scale_write_path <- paste(base_path, 'p0_BrainCortex/chromatin/scale_count_matrix.csv',
                         sep = "")
atac_scale <- atac[["ATAC"]]@scale.data
atac_scale <- t(atac_scale)
write.csv(atac_scale, scale_write_path)

dna<-Read10X(dna_path, gene.column = 1)
dna_seurat <- CreateSeuratObject(counts = dna, assay = "DNA", project = "10x_ATAC")
dna_seurat <- FindVariableFeatures(dna_seurat, selection.method = "vst", 
                                   nfeatures = 5000, verbose = TRUE)
dna_seurat <- subset(dna_seurat, features =c(VariableFeatures(dna_seurat) ))
dna_seurat <- NormalizeData(dna_seurat)
dna_seurat <- ScaleData(dna_seurat)

dna_norm <- GetAssayData(dna_seurat)
dna_norm <- t(dna_norm)
dna_norm_write_path <- paste(base_path, 'p0_BrainCortex/dna/norm_count_matrix.csv',
                         sep = "")
write.csv(dna_norm, dna_norm_write_path)

dna_scale_write_path <- paste(base_path, 'p0_BrainCortex/dna/scale_count_matrix.csv',
                          sep = "")
dna_scale <- dna_seurat[["DNA"]]@scale.data
dna_scale <- t(dna_scale)
write.csv(dna_scale, dna_scale_write_path)


geneActMat <- as.matrix(geneActMat)					 
geneActMat <- t(geneActMat)

write_path <- paste(base_path, 'p0_BrainCortex/dna/count_matrix.csv',
                   sep = "")
write.csv(peak, write_path)

label_path = '/Users/zhongyuanke/data/atac/p0_BrainCortex/P0BrainCortex_SNAREseq_metadata_full.rds'
celltype_path <- '/Users/zhongyuanke/data/atac/p0_BrainCortex/ref_barcodes.csv'
label <- readRDS(label_path)
ident <- label['Ident']
ident_char <- label['IdentChar']
ref_barcodes <- label['Barcode']
write.csv(ref_barcodes, celltype_path)

# ------------------------ h5ad -----------------------
library(Seurat)
library(ggplot2)
library(patchwork)

peaks <- Read10X_h5("/Users/zhongyuanke/data/seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")

# create a gene activity matrix
activity.matrix <- CreateGeneActivityMatrix(
  peak.matrix = peaks,
  annotation.file = "/Users/zhongyuanke/data/seurat_data/sc_atac/Homo_sapiens.GRCh37.82.gtf",
  seq.levels = c(1:22, "X", "Y"),
  upstream = 2000,
  verbose = TRUE
)

pbmc.atac <- CreateSeuratObject(counts = peaks, assay = "ATAC", project = "10x_ATAC")
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
meta <- read.table(
  "/Users/zhongyuanke/data/seurat_data/sc_atac/atac_v1_pbmc_10k_singlecell.csv",
  sep = ",",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)
meta <- meta[colnames(pbmc.atac), ]
pbmc.atac <- AddMetaData(pbmc.atac, metadata = meta)

# filter
pbmc.atac <- subset(pbmc.atac, subset = nCount_ATAC > 5000)
pbmc.atac$tech <- "atac"