import scbean.model.davae as davae
import scbean.tools.utils as tl
import scanpy as sc
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

adata = tl.read_sc_data('/Users/zhongyuanke/data/seurat_data/panc_8/panc8_8.h5ad')
print(adata)
datasets = tl.split_object(adata, by="replicate")

adata_all = tl.davae_preprocessing(datasets, n_top_genes=2000, sparse=True)
print(adata_all)
adata_integrate = davae.fit_integration(
    adata_all,
    batch_num=8,
    split_by='batch_label',
    batch_size=128,
    domain_lambda=100,
    epochs=30,
    sparse=True,
    hidden_layers=[128, 64, 32, 7],
)
print(adata_integrate)
sc.pp.neighbors(adata_integrate, use_rep='X_davae')
sc.tl.umap(adata_integrate)
sc.pl.umap(adata_integrate, color=['batch_label', 'celltype'], cmap='Set2',s=5)
adata_integrate.write_h5ad('/Users/zhongyuanke/data/dann_vae/panc8/davae_8.h5ad')