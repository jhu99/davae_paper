import scanpy as sc
import matplotlib.pyplot as plt
import ttools as tl
import numpy as np
import argparse
from sklearn.preprocessing import normalize
import umap
from sklearn.manifold import TSNE
from sklearn.metrics.cluster import adjusted_rand_score, silhouette_score
from sklearn.cluster import KMeans
import collections

base_path = '/Users/zhongyuanke/data/'

orig_path = 'dann_vae/mouse/results/orig_mouse.h5ad'
davae_path = 'dann_vae/mouse/results/davae_mouse_23_16_1k_1.h5ad'
scan_path = 'dann_vae/mouse/results/scan_mouse.h5ad'
seurat_path = 'seurat_result/results/seurat_mouse.h5ad'
desc_path = 'dann_vae/mouse/results/desc_mouse.h5ad'
scgen_path = 'dann_vae/mouse/results/scgen_mouse.h5ad'

adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
# adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)

orig_batches = adata_davae.obs['batch'].values
orig_batches=list(map(int,orig_batches))

label=adata_davae.obs['class']
orig_label=tl.text_label_to_number(label)
print(collections.Counter(orig_label))

# print(adata_seurat)
# sc.pp.neighbors(adata_seurat, use_rep='X_pca')
# sc.tl.umap(adata_seurat)
# adata_seurat.write_h5ad(base_path+seurat_path)


# ---------------------------------------
data_scgen_emb=adata_scgen.obsm['X_umap']
data_orig_emb = adata_orig.obsm['X_umap']
data_davae_emb = adata_davae.obsm['X_umap']
# data_seurat_emb = adata_seurat.obsm['X_umap']
data_scan_emb = adata_scan.obsm['X_umap']
data_desc_emb = adata_desc.obsm['X_umap0.8']

batch_cmap = 'Dark2'
c_map = 'tab20'
size = 1
xy_label_size = 10
title_size = 10
xy_label = 'umap'

#
fig = plt.figure(figsize=(40, 12))

ax = fig.add_subplot(251)
plt.title('Raw', fontsize=title_size)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


ax = fig.add_subplot(256)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ---------------------------------------

ax = fig.add_subplot(252)
plt.title('DAVAE', fontsize=title_size)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


ax = fig.add_subplot(257)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# ---------------------------------------

ax = fig.add_subplot(253)
plt.title('Scanorama', fontsize=title_size)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


ax = fig.add_subplot(258)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


# ---------------------------------------

ax = fig.add_subplot(254)
plt.title('DESC', fontsize=title_size)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)


ax = fig.add_subplot(2,5,9)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ---------------------------------------

# ax = fig.add_subplot(265)
# plt.title('Seurat 3.0', fontsize=title_size)
# ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# # plt.xlabel(xy_label + '1', size=xy_label_size)
# # plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(2, 6, 11)
# ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# -----------------------------------------------

ax = fig.add_subplot(255)
plt.title('scGen', fontsize=title_size)
ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2, 5, 10)
ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

plt.show()
# plt.savefig(figure_merge)
plt.close(fig)

# ---------------scgen--------------
# scgen_batch = adata_scgen.obs['batch_label'].values
# scgen_label = adata_scgen.obs['celltype'].values
# scgen_label=tl.text_label_to_number(scgen_label)
# fig = plt.figure(figsize=(10, 12))
# ax = fig.add_subplot(121)
# plt.title('scGen', fontsize=title_size)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=scgen_batch, cmap=batch_cmap, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(122)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=scgen_label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
#
# plt.show()
# # plt.savefig(figure_merge)
# plt.close(fig)
