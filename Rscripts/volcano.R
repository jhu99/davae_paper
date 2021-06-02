#  ------------------------FindMarkers---------------------------

library(Seurat)
library(SeuratDisk)

Convert("/Users/zhongyuanke/data/vipcca/mixed_cell_lines_result/output_save.h5ad",
        dest = "h5seurat", overwrite = TRUE)
mcl <- LoadH5Seurat("/Users/zhongyuanke/data/vipcca/mixed_cell_lines_result/output_save.h5seurat")

mcl <- split(mcl,mcl@meta.data[["celltype"]]) 

Idents(mcl)<-"X_batch"

Idents(mcl) <- 'celltype'
mcl$celltype.cond <- paste(Idents(mcl), mcl@meta.data[['X_batch']], sep = "_")
Idents(mcl) <- "celltype.cond"

br <- FindMarkers(mcl, ident.1 = '293t_293t', ident.2 = '293t_mixed', 
                  slot = "data",
                  logfc.threshold = 0.,
                  test.use = "wilcox", 
                  verbose = FALSE, 
                  min.cells.feature = 1,
                  min.cells.group = 1,
                  min.pct = 0.1)

boxplot(-log(br$p_val_adj [, v]),
        main = "Adjusted P-value for each gene",
        xlab = "293T",
        ylab = "Adjusted P-value")
br$p_val_adj <- as.numeric(br$p_val_adj)
class(br$avg_log2FC)

qtest<-createData(method = "test",pval = br$p_val)
p <- ggplot(qtest,
            aes(x=expected,y=observed,group=Method))

#  ------------------------EnhancedVolcano---------------------------

BiocManager::install('EnhancedVolcano')
devtools::install_github('kevinblighe/EnhancedVolcano')

library(EnhancedVolcano)
EnhancedVolcano(br,
                lab = rownames(br),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Volcano plot for 293t and mixed batch',
                )


# ---------------------spatial------------------------
library(Seurat)
library(SeuratDisk)
Convert("/Users/zhongyuanke/data/seurat_result/spatial/anterior_label_01.h5ad",
        dest = "h5seurat", overwrite = TRUE)
an_seurat<-LoadH5Seurat('/Users/zhongyuanke/data/seurat_result/spatial/anterior_label_01.h5seurat')
Convert("/Users/zhongyuanke/data/dann_vae/spatial/anterior_label_01.h5ad",
        dest = "h5seurat", overwrite = TRUE)
an_davae<-LoadH5Seurat('/Users/zhongyuanke/data/dann_vae/spatial/anterior_label_01.h5seurat')

spatial<-an_seurat
spatial<-an_davae
spatial <- split(spatial,spatial@meta.data[["celltype"]]) 


Idents(an_seurat) <- 'celltype'
br_seurat <- FindMarkers(an_seurat, ident.1 = 'L6 CT',
                  slot = "data",
                  logfc.threshold = 0.,
                  test.use = "wilcox", 
                  verbose = FALSE, 
                  min.cells.feature = 1,
                  min.cells.group = 1,
                  min.pct = 0.1)

Idents(an_davae) <- 'celltype'
br_davae <- FindMarkers(an_davae, ident.1 = 'L6 CT',
                         slot = "data",
                         logfc.threshold = 0.,
                         test.use = "wilcox", 
                         verbose = FALSE, 
                         min.cells.feature = 1,
                         min.cells.group = 1,
                         min.pct = 0.1)

boxplot(-log(br_seurat$p_val_adj),
        main = "Adjusted P-value for each gene",
        xlab = "L6 CT",
        ylab = "Adjusted P-value")

boxplot(-log(br_seurat$p_val_adj), -log(br_davae$p_val_adj),
        main = "Multiple boxplots for comparision",
        at = c(1,2),
        names = c("Seurat L6 CT", "DAVAE L6 CT"),
        las = 2,
        col = c("orange","red"),
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)


br$p_val_adj <- as.numeric(br$p_val_adj)
class(br$avg_log2FC)

qtest<-createData(method = "test",pval = br$p_val)
p <- ggplot(qtest,
            aes(x=expected,y=observed,group=Method))
BiocManager::install('EnhancedVolcano')
devtools::install_github('kevinblighe/EnhancedVolcano')

library(EnhancedVolcano)
EnhancedVolcano(br_davae,
                lab = rownames(br_davae),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Volcano plot for L5 PT',
)
EnhancedVolcano(br_seurat,
                lab = rownames(br_seurat),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Volcano plot for L5 PT',
)



