import scanpy as sc

file_rna = '/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3_filtered_feature_bc_matrix.h5'
adata = sc.read_10x_h5(file_rna)
adata.var_names_make_unique()
adata.obs_names_make_unique()
print(adata)
adata.write_h5ad('/Users/zhongyuanke/data/seurat_data/sc_atac/pbmc_10k_v3_filtered_feature_bc_matrix.h5ad')
