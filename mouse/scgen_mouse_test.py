import numpy as np
from scbean.tools import utils as tl
from scbean.model import davae as davae
import anndata
import scanpy as sc
import pandas as pd
import scgen
base_path = '/Users/zhongyuanke/data/'

file1 = 'dropviz/mouse_brain_dropviz_filtered.h5ad'
file2 = 'nuclei/adata_nuclei_filtered.h5ad'
scgen_path = 'results/scgen_mouse.h5ad'

# -------------train---------------------
adata1 = tl.read_sc_data(file1, fmt='h5ad')
adata2 = tl.read_sc_data(file2, fmt='h5ad')
print(adata1)
print(adata2)

adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000)
adata_all.obs_names_make_unique()
adata_all = scgen.setup_anndata(adata_all, batch_key="batch_label", copy=True)
model = scgen.SCGEN(adata_all)
model.train(max_epochs=10,
    batch_size=256,
    early_stopping=True,
    early_stopping_patience=25,
    use_gpu=True)
corrected_adata = model.batch_removal()
sc.pp.neighbors(corrected_adata,use_rep='corrected_latent')
sc.tl.umap(corrected_adata)
# sc.pl.umap(adata_integrate, color=['batch_label', 'class'], s=10, cmap='Dark2')
# print(adata_integrate)
# adata_integrate.write_h5ad(base_path+'dann_vae/benchmark1/dc_davae_temp.h5ad')
corrected_adata.write_h5ad(scgen_path)






