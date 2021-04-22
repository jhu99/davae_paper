import numpy as np
import scbean.model.davae as davae
from scbean.tools import utils as tl
import scanpy as sc
import argparse
import matplotlib
matplotlib.use('TkAgg')

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/', help="base path")
parser.add_argument("--epoch", type=int, default=10, help="epoch of training")

opt = parser.parse_args()

base_path = opt.base_path
out_path = 'dann_vae/hca/davae_01.h5ad'
file1 = base_path+'HCA/ica_cord_blood_h5.h5'
file2 = base_path+'HCA/ica_bone_marrow_h5.h5'
adata1 = tl.read_sc_data(file1, fmt='10x_h5')
adata2 = tl.read_sc_data(file2, fmt='10x_h5')
adata1.var_names_make_unique()
adata2.var_names_make_unique()
print(adata1)

adata_all = tl.davae_preprocessing([adata1, adata2], hvg=False, lognorm=False)
adata_integrate = davae.fit_integration(adata_all, split_by='batch',
                                  domain_lambda=5,
                                  epochs=1,
                                  hidden_layers=[128, 64, 32, 5],
                                  sparse=True,
                                        )
sc.pp.neighbors(adata_integrate, use_rep='X_davae')
sc.tl.umap(adata_integrate)
sc.pl.umap(adata_integrate, color='batch')
print(adata_integrate)
