library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
library(zoo)
library(dplyr)
library(gggenes)
set.seed(1234)

base_path <- '/Users/zhongyuanke/data/seurat_data/'
counts_path <- paste(base_path, 'sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5' ,sep = "")
meta_path <- paste(base_path, 'sc_atac/atac_v1_pbmc_10k_singlecell.csv' ,sep = "")
rna_path <- paste(base_path, 'sc_atac/pbmc_10k_v3.rds' ,sep = "")
fragment.path <- paste(base_path, 'sc_atac/atac_v1_pbmc_10k_fragments.tsv.gz' ,sep = "")
# source('./r_script/signac_function.R')

counts <- Read10X_h5(filename = counts_path)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = fragment.path,
  min.cells = 10,
  min.features = 200
)

metadata <- read.csv(
  file = meta_path,
  header = TRUE,
  row.names = 1
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  meta.data = metadata,
)

# ------------------Annotation--------------------------------

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg19"
# add the gene information to the object
Annotation(pbmc) <- annotations
# ------------------   QC--------------------------------

pbmc <- NucleosomeSignal(object = pbmc)
# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)
# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
  peak_region_fragments < 20000 &
  pct_reads_in_peaks > 15 &
  blacklist_ratio < 0.05 &
  nucleosome_signal < 10 &
  TSS.enrichment > 2
)

#--------------------------------------------------------------------------
# pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(
  object = pbmc,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)


pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3)

# ------------- create gene activity ---------------------------------
gene.activities <- GeneActivity(pbmc)

pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)

pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)


#----------------------------------------------
DefaultAssay(pbmc) <- 'RNA'
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 12000)
pbmc <- subset(pbmc, features =c(VariableFeatures(pbmc) ))

data<-GetAssayData(pbmc)
data <- as.matrix(data)
data <- t(data)
data_path <- '/Users/zhongyuanke/data/seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix_norm.csv'
write.csv(data, data_path)


#--------------transfer----------------------------
DefaultAssay(pbmc) <- 'RNA'
pbmc_rna <- readRDS('/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3.rds')

transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc[['lsi']],
  dims = 2:30
)
write.csv(predicted.labels[['predicted.id']], '/Users/zhongyuanke/data/dann_vae/atac/seurat_pred_type.csv')
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

#------------------------------------------
mylabel = read.csv('/Users/zhongyuanke/data/dann_vae/atac/pred_type_save03.csv', header = FALSE)
mylabel <- as.factor(unlist(mylabel))
list1 <-rownames(predicted.labels)
names(mylabel) <- list1
mylabel
pbmc@active.ident<-mylabel


#------------------------------------------

pbmc <- AddMetaData(object = pbmc, metadata = mylabel, col.name = 'active.ident')

pbmc <- subset(pbmc,idents =13, invert = TRUE)

pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'B cell progenitor',
  '1' = 'CD14+ Monocytes',
  '2' = 'CD16+ Monocytes',
  '3' = 'CD4 Memory',
  '4' = 'CD4 Naive',
  '5' = 'CD8 Naive',
  '6' = 'CD8 effector',
  '7' = 'Dendritic cell',
  '8' = 'Double negative T cell',
  '9' = 'NK cell',
  '10' = 'Platelets',
  '11' = 'pDC',
  '12' = 'pre-B cell'
)
pbmc <- subset(pbmc, idents = 14, invert = TRUE)
pbmc <- RenameIdents(
  object = pbmc,
  '0' = 'CD14 Mono',
  '1' = 'CD4 Memory',
  '2' = 'CD8 Effector',
  '3' = 'CD4 Naive',
  '4' = 'CD14 Mono',
  '5' = 'DN T',
  '6' = 'CD8 Naive',
  '7' = 'NK CD56Dim',
  '8' = 'pre-B',
  '9' = 'CD16 Mono',
  '10' = 'pro-B',
  '11' = 'DC',
  '12' = 'NK CD56bright',
  '13' = 'pDC'
)

DefaultAssay(pbmc) <- 'RNA'
p1 <- DimPlot(pbmc, reduction = "umap")


DefaultAssay(pbmc_right) <- 'peaks'
DefaultAssay(pbmc) <- 'peaks'
da_peaks <- FindMarkers(
  object = pbmc,
  ident.1 = "Platelets",
  # ident.2 = "CD14 Mono",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)


open_cd4naive <- rownames(da_peaks[da_peaks$avg_logFC > 0.5, ])
open_cd14mono <- rownames(da_peaks[da_peaks$avg_logFC < -0.5, ])

closest_genes_cd4naive <- ClosestFeature(
  regions = open_cd4naive,
  annotation = gene.ranges,
  sep = c(':', '-')
)
closest_genes_cd14mono <- ClosestFeature(
  regions = open_cd14mono,
  annotation = gene.ranges,
  sep = c(':', '-')
)
rownames(da_peaks)[c(1,5)]

levels(pbmc) <- c('B cell progenitor','CD14+ Monocytes','CD16+ Monocytes','CD4 Memory',
'CD4 Naive', 'CD8 Naive', 'CD8 effector', 'Dendritic cell','Double negative T cell',
'NK cell', 'Platelets', 'pDC', 'pre-B cell')

levels(pbmc) <- c('B cell progenitor','CD14+ Monocytes','CD16+ Monocytes','CD4 Memory',
                  'CD4 Naive', 'CD8 Naive', 'CD8 effector', 
                  'Dendritic cell','Double negative T cell','pDC')

chr1:6896024-6929974

a<- c("chr12:69742121-69748015","chr2:85912298-85925978","chr11:60223225-60238234",
      "chr14:99635624-99737862","chr2:87011729-87035520","chr11:118175260-118186891",
      "chr12:6896024-6929975","chr5:140011313-140013260")
my_marker<-c('chr16:28943260-28950667','chr4:15778275-15853232',
             'chr11:67202981-67205538')
compose_markers<-c("chr12:69742121-69748015","chr2:85912298-85925978",
                   "chr14:99635624-99737862","chr2:87011729-87035520","chr11:118175260-118186891",
                   "chr12:6896024-6929975",'chr16:28943260-28950667',
                   'chr4:15778275-15853232',
                   'chr11:67202981-67205538')
CoveragePlot(
  object = pbmc_right,
  region = a,
  sep = c(":", "-"),
  peaks = StringToGRanges(rownames(pbmc), sep = c(":", "-")),
  annotation = gene.ranges,
  extend.upstream = 0,
  extend.downstream = 0,
  ncol = 9,
)

# marker_index <- which(gene.ranges@ranges@NAMES == "ENSG00000156738")
# start<-gene.ranges@ranges@start[marker_index]
# end<-gene.ranges@ranges@width[marker_index]
# start
# start+end

# 
# lyz: ENSG00000090382:chr12:69742121-69748015
# gnly: ENSG00000115523: chr2:85912298-85925978
# ms4a1: ENSG00000156738:chr11:60223225-60238234
# bcl11b: ENSG00000127152: chr14:99635624-99737862
# cd8a: ENSG00000153563: chr2:87011729-87035520
# cd3e: ENSG00000198851: chr11:118175260-118186891
# cd4: ENSG00000010610: chr12:6896024-6929975
# cd14: ENSG00000170458: chr5:140011313-140013260


# save(pbmc,file='/Users/zhongyuanke/Desktop/r_para/pbmc_seurat.Rdata')
# save(gene.ranges,file='/Users/zhongyuanke/Desktop/r_para/geneRanges.Rdata')
save(predicted.labels,file='/Users/zhongyuanke/Desktop/r_para/predlabels.Rdata')
# 
# load('/Users/zhongyuanke/Desktop/r_para/geneRanges.Rdata')
# 
load('/Users/zhongyuanke/Desktop/r_para/pbmc_seurat.Rdata')
pbmc_right <-pbmc
load(file='/Users/zhongyuanke/Desktop/r_para/pbmc_seurat.Rdata')


load(file='/Users/zhongyuanke/data/seurat_data/sc_atac/win_pbmc.Rdata')

SaveH5Seurat(integrated, filename = "/Users/zhongyuanke/data/seurat_result/atac.h5Seurat")
Convert("/Users/zhongyuanke/data/seurat_result/atac.h5Seurat", dest = "h5ad")

