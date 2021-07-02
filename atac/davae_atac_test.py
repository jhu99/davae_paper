from scbean.tools import utils as tl
import numpy as np

from scbean.model import davae as davae
import anndata
from sklearn.utils import shuffle
import scanpy as sc
import matplotlib
matplotlib.use('TkAgg')


base_path = '/Users/zhongyuanke/data/'
out_path = 'dann_vae/atac/davae_02.h5ad'

file1 = base_path+'seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix_8000.csv'
file2 = base_path+'seurat_data/sc_atac/pbmc_10k_v3_8000.csv'


adata1 = sc.read_csv(file1)
adata2 = sc.read_csv(file2)

print(adata2)

adata_list = [adata1, adata2]
# adata_davae = da_cvae.integration(adata_list, latent_size=10, epoch=30)
len1 = adata1.shape[0]
len2 = adata2.shape[0]
#
# print(adata1.shape)
# print(adata2.shape)
# adata = adata1.concatenate(adata2)
# adata.write_h5ad('/Users/zhongyuanke/data/dann_vae/atac/orig.h5ad')
# print('1')

# orig_data = adata.X
# sc.pp.filter_genes(data1, min_cells=3)
# sc.pp.filter_genes(data2, min_cells=3)
# b1 = np.array([0, 1])
# b2 = np.array([1, 0])
# orig_batches = []
# orig_batch_label = []
# loss_weight = []
# for i in range(len1):
#     orig_batches.append(b1)
#     orig_batch_label.append(0)
#     loss_weight.append(0)
# for i in range(len2):
#     orig_batches.append(b2)
#     orig_batch_label.append(1)
#     loss_weight.append(1)
#
# orig_batches = np.array(orig_batches)
# loss_weight = np.array(loss_weight)
# print(adata1.shape)
# print(adata2.shape)
# print(orig_data.shape)

# data, batches, batch_label, loss_weight = shuffle(orig_data, orig_batches, orig_batch_label, loss_weight,
#                                             random_state=0)

# net_x = DACVAE(input_size=data.shape[1], batches=2, latent_size=10)
# net_x.build()
# net_x.compile()
# his = net_x.train(data, batches, loss_weight, epochs=epochs, batch_size=batch_size)
#
# mid = net_x.embedding(orig_data, orig_batches)
# adata_mid = anndata.AnnData(X=mid)
# adata_mid.obs['batch'] = orig_batch_label
# adata.obsm['davae'] = mid
adata_all = tl.davae_preprocessing(adata_list, sparse=False)
adata_out = davae.fit_integration(adata_all, split_by='batch', mode='DACVAE',
                                  domain_lambda=5.0,
                                  epochs=3, hidden_layers=[128, 64, 32, 5],sparse=False)
sc.pp.neighbors(adata_out, use_rep='X_davae')
sc.tl.umap(adata_out)
sc.pl.umap(adata_out, color='celltype')
print(adata_out)
# adata.write_h5ad(base_path+out_path)
