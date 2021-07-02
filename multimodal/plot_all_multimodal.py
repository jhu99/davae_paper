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

# ---------panc8---------------------
# cluster_count = 15
# base_path='/Users/zhongyuanke/data/'
# orig_path = '/seurat_data/panc_8/panc8.h5ad'
# desc_path = 'desc/desc_panc8.h5ad'
# davae_path = 'dann_vae/panc8/davae_5.h5ad'
# seurat_path = 'seurat_result/panc8_all.h5ad'
# scan_path = 'scanorama/scan_panc8_all.h5ad'
# scgen_path='scgen/scgen_panc8.h5ad'
# adata_davae = sc.read_h5ad(base_path+davae_path)
# adata_scan = sc.read_h5ad(base_path+scan_path)
# adata_orig = sc.read_h5ad(base_path+orig_path)
# adata_seurat = sc.read_h5ad(base_path+seurat_path)
# adata_scgen=sc.read_h5ad(base_path+scgen_path)
# orig_label =adata_orig.obs['celltype'].values
# orig_label=tl.text_label_to_number(orig_label)
# orig_batches=adata_orig.obs['batch_label'].values
# print(orig_label)
# sc.pp.neighbors(adata_orig)
# sc.tl.umap(adata_orig)
# sc.pp.neighbors(adata_davae, use_rep='X_davae')
# sc.tl.umap(adata_davae)
# sc.pp.neighbors(adata_scan, use_rep='X_scanorama')
# sc.tl.umap(adata_scan)
# sc.pp.neighbors(adata_seurat)
# sc.tl.umap(adata_seurat)
# sc.pp.neighbors(adata_scgen, use_rep='corrected_latent')
# sc.tl.umap(adata_scgen)

orig_path = 'dann_vae/multimodal/multi_orig.h5ad'
base_path = '/Users/zhongyuanke/data/'
davae_path = 'dann_vae/multimodal/davae_multi_save.h5ad'
scan_path = 'scanorama/scan_multimodal.h5ad'
seurat_path = 'seurat_result/seurat_multimodal.h5ad'
desc_path = 'desc/desc_multimodal.h5ad'
scgen_path = 'scgen/scgen_multimodal.h5ad'

seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128
label=tl.get_label_by_txt(seurat_celltype_path)
label=label+label
print(collections.Counter(label))
orig_label=tl.text_label_to_number(label)
print(collections.Counter(orig_label))

adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)
# print(adata_orig)
# adata_orig.obs['celltype']=

# sc.pp.neighbors(adata_orig)
# sc.tl.umap(adata_orig)
# adata_orig.write_h5ad(base_path+orig_path)
import collections
# print(adata_davae.obs['celltype'])
# orig_label =adata_davae.obs['celltype']
# # sc.pl.umap(adata_orig, color=['celltype'], cmap='tab20b')
# orig_label = list(map(int, tl.text_label_to_number(orig_label)))
orig_batches = adata_davae.obs['batch'].values
orig_batches=list(map(int,orig_batches))
# davae_label =adata_davae.obs['label']
# davae_label = list(map(int, tl.text_label_to_number(davae_label)))
# davae_batches = adata_davae.obs['batch_label'].values

# orig_label=davae_label
# sc.pp.neighbors(adata_davae)
# sc.tl.umap(adata_davae)
# adata_davae.write_h5ad(base_path+'dann_vae/benchmark1/davae_save01_uamp.h5ad')
# sc.pp.neighbors(adata_scan,use_rep='X_scanorama')
# sc.tl.umap(adata_scan)
# adata_scan.write_h5ad(base_path+scan_path)
# #
# sc.pp.neighbors(adata_desc, use_rep='X_Embeded_z0.8')
# sc.tl.umap(adata_desc)
# adata_desc.write_h5ad(base_path+desc_path)
# #
# print(adata_seurat)
print(adata_seurat)
sc.pp.neighbors(adata_seurat, use_rep='X_pca')
sc.tl.umap(adata_seurat)
# adata_seurat.write_h5ad(base_path+seurat_path)
#
# sc.pp.neighbors(adata_scgen, use_rep='corrected_latent')
# sc.tl.umap(adata_scgen)
# adata_scgen.write_h5ad(base_path+scgen_path)





# adata_orig.write_h5ad(base_path+'pbmc/zheng/mcl_pre.h5ad')
# print(adata_orig)

# ---------------------------------------
# from matplotlib.pyplot import rc_context
# with rc_context({'figure.figsize': (3, 3)}):
#     sc.pl.umap(adata_orig)
#     sc.pl.umap(adata_davae)
#     sc.pl.umap(adata_scan)
#     sc.pl.umap(adata_seurat)

# adata_desc = sc.read_h5ad(base_path + desc_path)

data_scan = adata_scan.X
# data_davae = adata_davae.obsm['davae']
# data_davae = adata_davae.X
# data_desc = adata_desc.X
# data_desc_emb = adata_desc.obsm['X_umap0.8']
# data_orig_emb = adata_orig.obsm['umap']
data_orig = adata_orig.X
data_seurat = adata_seurat.X
# data_harmony=adata_harmony.X
data_scgen_emb=adata_scgen.obsm['X_umap']

#data_orig_emb = umap.UMAP().fit_transform(data_orig)
# data_davae_emb = umap.UMAP().fit_transform(data_davae)
# data_seurat_emb = umap.UMAP().fit_transform(data_seurat)
# data_scan_emb = umap.UMAP().fit_transform(data_scan)
# data_harmony_emb = umap.UMAP().fit_transform(data_harmony)
data_orig_emb = adata_orig.obsm['X_umap']
data_davae_emb = adata_davae.obsm['X_umap']
data_seurat_emb = adata_seurat.obsm['X_umap']
data_scan_emb = adata_scan.obsm['X_umap']
data_desc_emb = adata_desc.obsm['X_umap']

batch_cmap = 'Dark2'
c_map = 'tab20b'
size = 1
xy_label_size = 10
title_size = 10
xy_label = 'umap'

#
fig = plt.figure(figsize=(40, 12))

ax = fig.add_subplot(261)
plt.title('Raw', fontsize=title_size)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)

# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(267)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(262)
plt.title('DAVAE', fontsize=title_size)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(268)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(263)
plt.title('Scanorama', fontsize=title_size)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(269)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(264)
plt.title('DESC', fontsize=title_size)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2,6,10)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(265)
plt.title('Seurat 3.0', fontsize=title_size)
ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2, 6, 11)
ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# -----------------------------------------------
# ax = fig.add_subplot(266)
# plt.title('Harmony', fontsize=title_size)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(2, 6, 12)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_label, cmap=c_map, s=size)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

# ----------------------------------------
# scgen_batch = adata_scgen.obs['batch_label'].values
# scgen_label = adata_scgen.obs['celltype'].values
# scgen_label = tl.text_label_to_number(scgen_label)
# print(scgen_label)
# print(orig_label)

ax = fig.add_subplot(266)
plt.title('scGen', fontsize=title_size)
ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2, 6, 12)
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
