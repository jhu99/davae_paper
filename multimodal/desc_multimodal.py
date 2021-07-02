import desc
import numpy as np
from scbean.tools import utils as tl
from scbean.model import davae as davae
import anndata
import scanpy as sc
import pandas as pd
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
adata_all.X=adata_all.X.A

adata_out = desc.train(adata_all, dims=[adata_all.shape[1], 32, 16], tol=0.005, n_neighbors=20,
               batch_size=64,
               save_encoder_weights=False)
print(adata_out)
adata_out.write_h5ad(base_path+'desc/desc_multimodal.h5ad')
