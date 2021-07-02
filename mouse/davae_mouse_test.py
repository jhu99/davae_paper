import numpy as np
from scbean.tools import utils as tl
from scbean.model import davae as davae
import anndata
import scanpy as sc
import pandas as pd
import argparse
import tensorflow as tf

parser = argparse.ArgumentParser()
parser.add_argument("--dlambda", type=float, default=4, help="domain adversary")
parser.add_argument("--epoch", type=int, default=40, help="epochs")
parser.add_argument("--outlayer", type=int, default=10, help="outlayer size")
parser.add_argument("--topgene", type=int, default=2000, help="top gene number")
parser.add_argument("--save_path", type=str, default='results/davae_mouse.h5ad', help="save path")
parser.add_argument("--n_neib", type=int, default=10, help="n neighbor")
opt = parser.parse_args()

epoch=opt.epoch
dlambda=opt.dlambda
outlayer=opt.outlayer
topgene=opt.topgene

# gpus = tf.config.experimental.list_physical_devices(device_type='GPU')
# print(gpus)
#
# tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
# tf.config.experimental.set_memory_growth(gpus[0], True)
# for gpu in gpus:
#     tf.config.experimental.set_memory_growth(gpu, True)
# epoch=50
# dlambda=5
# outlayer=6
# topgene=3000
base_path = '/Users/zhongyuanke/data/'

file1 = 'dropviz/mouse_brain_dropviz_filtered.h5ad'
file2 = 'nuclei/adata_nuclei_filtered.h5ad'
davae_path = opt.save_path
print(davae_path)
# file1 = '/Users/zhongyuanke/data/mouse_brain/mouse_brain_dropviz_filtered.h5ad'
# file2 = '/Users/zhongyuanke/data/mouse_brain/adata_nuclei_filtered.h5ad'
# davae_path = '/Users/zhongyuanke/data/dann_vae/mouse/davae_mouse.h5ad'
# -------------train---------------------
adata1 = tl.read_sc_data(file1, fmt='h5ad')
adata2 = tl.read_sc_data(file2, fmt='h5ad')
adata1.write_h5ad(davae_path)
# orig_label =adata_orig.obs['label']
print(adata1)
print(adata2)
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=topgene)
# sc.pp.neighbors(adata_all,use_rep='X')
# sc.tl.umap(adata_all)
# adata_all.write_h5ad('results/davae_orig.h5ad')
# print('finish orig')

adata_integrate = davae.fit_integration(adata_all, split_by='batch', epochs=epoch,
                                        hidden_layers=[128, 64, 32, 16],
                                        domain_lambda=dlambda)

sc.pp.neighbors(adata_integrate, use_rep='X_davae',n_neighbors=opt.n_neib)
sc.tl.umap(adata_integrate)
# sc.pl.umap(adata_integrate, color=['batch_label', 'class'], s=10, cmap='Dark2')
# print(adata_integrate)
# adata_integrate.write_h5ad(base_path+'dann_vae/benchmark1/dc_davae_temp.h5ad')
adata_integrate.write_h5ad(davae_path)






