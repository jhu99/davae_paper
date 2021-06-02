import numpy as np
from scbean.tools import utils as tl
from scbean.model import davae as davae
import anndata
import scanpy as sc
import pandas as pd
import matplotlib
matplotlib.use('TkAgg')

epochs = 40
base_path = '/Users/zhongyuanke/data/'

file1 = base_path+'dann_vae/benchmark1/dc_batch1.h5ad'
file2 = base_path+'dann_vae/benchmark1/dc_batch2.h5ad'
orig_path = base_path+'dann_vae/benchmark1/orig.h5ad'
# -------------train---------------------
adata1 = tl.read_sc_data(file1, fmt='h5ad')
adata2 = tl.read_sc_data(file2, fmt='h5ad')
adata_orig = tl.read_sc_data(orig_path, fmt='h5ad')
# orig_label =adata_orig.obs['label']
print(adata1)
print(adata2)
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=4000, sparse=False)
adata_integrate = davae.fit_integration(adata_all, split_by='batch', epochs=1000,
                                        hidden_layers=[128, 64, 32, 2],
                                        sparse=False,
                                        domain_lambda=6)
adata_integrate.obs['label']=adata_orig.obs['label']
sc.pp.neighbors(adata_integrate, use_rep='X_davae')
sc.tl.umap(adata_integrate)
sc.pl.umap(adata_integrate, color=['batch','label'], s=10, cmap='Dark2')
print(adata_integrate)
adata_integrate.write_h5ad(base_path+'dann_vae/benchmark1/dc_davae_temp.h5ad')







