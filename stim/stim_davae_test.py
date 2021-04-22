import scbean.model.davae as davae
import scbean.tools.utils as tl
import scanpy as sc
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

adata = tl.read_sc_data('/Users/zhongyuanke/data/seurat_data/ifnb/ifnb.h5ad')
datasets = tl.split_object(adata, by="stim")
print(datasets[0])
print(datasets[1])
adata_all = tl.davae_preprocessing(datasets, n_top_genes=8000)
adata_intagrate = davae.fit_integration(adata_all, epochs=30, hidden_layers=[128, 64, 32, 5], domain_lambda=3.0,)
print(adata_intagrate)
sc.pp.neighbors(adata_intagrate, use_rep='X_davae', n_neighbors=15)
sc.tl.louvain(adata_intagrate)
sc.tl.umap(adata_intagrate)
sc.pl.umap(adata_intagrate, color='louvain', cmap='tab20c')
