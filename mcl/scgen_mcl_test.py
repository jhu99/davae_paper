import scgen
import scanpy as sc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--result_path", type=str, default='/Users/zhongyuanke/data/scgen/scgen_panc8.h5ad', help="base path")
parser.add_argument("--file_path", type=str, default='/Users/zhongyuanke/data/seurat_data/panc_8/panc8.h5ad', help="base path")

import scbean.model.davae as davae
import scbean.tools.utils as tl
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

base_path = "/Users/zhongyuanke/data/vipcca/mixed_cell_lines/"

adata_b1 = sc.read_h5ad(base_path+"293t.h5ad")
adata_b2 = sc.read_h5ad(base_path+"jurkat.h5ad")
adata_b3 = sc.read_h5ad(base_path+"mixed.h5ad")
# adata_b1.obs_names_make_unique()
# adata_b2.obs_names_make_unique()
# adata_b3.obs_names_make_unique()

adata_all = tl.davae_preprocessing([adata_b1, adata_b2, adata_b3], n_top_genes=3000)
adata_all.obs_names_make_unique()

adata_all = scgen.setup_anndata(adata_all, batch_key="batch_label", copy=True)
model = scgen.SCGEN(adata_all)
model.train(max_epochs=25,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
    use_gpu=False)
corrected_adata = model.batch_removal()
sc.pp.neighbors(corrected_adata,use_rep='corrected_latent')
sc.tl.umap(corrected_adata)
sc.pl.umap(corrected_adata,color='celltype')

# corrected_adata.write('/Users/zhongyuanke/data/scgen/scgen_mcl01.h5ad')
