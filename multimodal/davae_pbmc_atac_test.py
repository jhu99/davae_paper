import scanpy as sc
import scbean.model.davae as davae
import scbean.tools.utils as tl
import matplotlib
import ttools as tool
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = '/Users/zhongyuanke/data/dann_vae/multimodal/rna.h5ad'
file_atac = '/Users/zhongyuanke/data/dann_vae/multimodal/atac.h5ad'
seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128

adata1 = sc.read_h5ad(file_atac)
adata2 = sc.read_h5ad(file_rna)
print(adata1)
print(adata2)
# sc.pp.filter_genes(adata1, min_cells=100)
# sc.pp.filter_genes(adata2, min_cells=100)
# sc.pp.log1p(adata1)
# sc.pp.log1p(adata2)
# sc.pp.scale(adata2)
# adata2.obs['celltype'] = adata1.obs['celltype']

# adata2.write_h5ad(base_path + 'multimodal/atac_pbmc_10k/activaty_matrix_label.h5ad')
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000, hvg=False, lognorm=False)
# sc.pp.scale(adata_all)
print(adata_all)
adata_integrate = davae.fit_integration(
    adata_all,
    # mode='DAVAE',
    batch_num=2,
    split_by='batch_label',
    domain_lambda=4.0,
    epochs=60,
    sparse=True,
    hidden_layers=[128, 64, 32, 5]
)
print(adata_integrate)
# adata_integrate.obs['celltyp']=adata
import umap
sc.pp.neighbors(adata_integrate, use_rep='X_davae')
sc.tl.umap(adata_integrate)
# adata_integrate.write_h5ad('/Users/zhongyuanke/data/dann_vae/multimodal/davae_multi_temp.h5ad')
sc.pl.umap(adata_integrate, color=['batch_label'], s=3, cmap='tab20c')
#
# sc.pl.umap(adata_integrate, color=['CD79A', 'CD8B', 'LYZ', 'IL32', 'CD3D', 'CST7', 'LEF1', 'FCER1G', 'PLD4', 'FTL'],
#            s=2, ncols=4)
