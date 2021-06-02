import numpy as np
from scbean.tools import utils as tl
from scbean.model import davae as davae
import anndata
import scanpy as sc
import pandas as pd
import scgen
import matplotlib
matplotlib.use('TkAgg')
base_path = '/Users/zhongyuanke/data/'

file1 = base_path+'dann_vae/benchmark1/dc_batch1.h5ad'
file2 = base_path+'dann_vae/benchmark1/dc_batch2.h5ad'

# -------------train---------------------
adata1 = tl.read_sc_data(file1, fmt='h5ad')
adata2 = tl.read_sc_data(file2, fmt='h5ad')
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000, sparse=False)
adata_all = scgen.setup_anndata(adata_all, batch_key="batch_label",  copy=True)
model = scgen.SCGEN(adata_all)
model.train(max_epochs=30,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
    use_gpu=False)
corrected_adata = model.batch_removal()
sc.pp.neighbors(corrected_adata,use_rep='corrected_latent')
sc.tl.umap(corrected_adata)
corrected_adata.write('/Users/zhongyuanke/data/scgen/scgen_dc.h5ad')