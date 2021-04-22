import scanpy as sc
import matplotlib
matplotlib.use('TkAgg')

adata=sc.read_10x_h5('/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3_filtered_feature_bc_matrix.h5')
print(adata)
# sc.pl.umap(adata, color='orig.ident', s=3)