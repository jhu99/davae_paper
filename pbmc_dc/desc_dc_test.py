import desc
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
# adata_all=adata1.concatenate(adata2)

adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=5000, sparse=False)
# desc.normalize_per_cell(adata_all, counts_per_cell_after=1e4)

adata_out = desc.train(adata_all, dims=[adata_all.shape[1], 32, 16], tol=0.005, n_neighbors=10,
               batch_size=64,
               save_encoder_weights=False)
print(adata_out)
adata_out.write_h5ad(base_path+'desc/desc_dc.h5ad')
