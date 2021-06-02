import scanpy as sc
import scbean.model.davae as davae
import scbean.tools.utils as tl
import matplotlib
import ttools as tool
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = base_path + '/multimodal/atac_pbmc_10k/rna.h5ad'
file_atac = base_path + 'multimodal/atac_pbmc_10k/activaty_matrix.h5ad'
seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128

# adata1 = sc.read_h5ad(file_rna)
# adata2 = sc.read_h5ad(file_atac)
# adata1.obs_names_make_unique()
# adata2.obs_names_make_unique()

# adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000)
file_seurat = '/Users/zhongyuanke/data/seurat_result/multimodal.h5ad'
adata = sc.read_h5ad(file_seurat)
print(adata)
batch = [0]*10412+[1]*10412
adata.obs['batch_label'] = batch
sc.pl.umap(adata, color='batch_label')