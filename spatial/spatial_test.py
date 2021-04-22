from scbean.tools import utils as tl
import numpy as np
import scbean.model.davae as davae
import anndata
from sklearn.utils import shuffle
from keras.utils import to_categorical
import scanpy as sc
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')


base_path = '/Users/zhongyuanke/data/'
out_path = 'dann_vae/spatial/davae_save02.h5ad'
file1 = base_path+'spatial/mouse_brain/10x_mouse_brain_Anterior/V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5'
file2 = base_path+'spatial/mouse_brain/10x_mouse_brain_Posterior/V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5'
file1_p = base_path + 'spatial/10x_mouse_brain_Anterior/anterior.h5ad'
file2_p = base_path + 'spatial/10x_mouse_brain_Posterior/posterior.h5ad'
batch_size = 256
epochs = 25

# adata1 = sc.read_10x_h5(file1)
adata1 = tl.read_sc_data(file1, fmt='10x_h5', batch_name='Anterior')
adata2 = tl.read_sc_data(file2, fmt='10x_h5', batch_name='Posterior')
print(adata1)
# print(adata1)
adata1.var_names_make_unique()
adata2.var_names_make_unique()
# len1 = adata1.shape[0]
# len2 = adata2.shape[0]

# sc.pp.filter_genes(adata1, min_cells=30)
# sc.pp.filter_genes(adata2, min_cells=30)
# sc.pp.log1p(adata1)
# sc.pp.log1p(adata2)
# print(adata1)
# print(adata2)
#
# sc.pp.highly_variable_genes(adata1, n_top_genes=6000)
# sc.pp.highly_variable_genes(adata2, n_top_genes=6000)
# adata1 = adata1[:, adata1.var.highly_variable]
# adata2 = adata2[:, adata2.var.highly_variable]
#
# adata1.write_h5ad(file1_p)
# adata2.write_h5ad(file2_p)
# del adata1.var['highly_variable']
# del adata2.var['highly_variable']
# del adata1.var['means']
# del adata2.var['means']
# del adata1.var['dispersions']
# del adata2.var['dispersions']
# del adata1.var['dispersions_norm']
# del adata2.var['dispersions_norm']
print(adata1)
print(adata2)
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=4000)
adata_integrate = davae.fit_integration(
    adata_all,
    epochs=25,
    hidden_layers=[128, 64, 32, 5],
    sparse=True,
    domain_lambda=0.5,
)
# import umap
# adata_integrate.obsm['X_umap']=umap.UMAP().fit_transform(adata_integrate.obsm['X_davae'])
sc.pp.neighbors(adata_integrate, use_rep='X_davae', n_neighbors=8)
sc.tl.umap(adata_integrate)
sc.pl.umap(adata_integrate, color=['_batch'], s=3)
adata_integrate.write_h5ad(base_path+out_path)

