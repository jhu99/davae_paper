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
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

base_path='/Users/zhongyuanke/data/'
orig_path = 'dann_vae/ifnb/orig.h5ad'
davae_path = 'dann_vae/ifnb/davae_save01.h5ad'
scan_path = 'scanorama/ifnb.h5ad'
seurat_path = 'seurat_result/ifnb.h5ad'
desc_path = 'desc/ifnb.h5ad'
scgen_path = 'scgen/scgen_ifnb.h5ad'


adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen=sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)
sc.pp.filter_genes(adata_orig, min_cells=10)
print(adata_orig)

sc.pp.neighbors(adata_orig)
sc.tl.umap(adata_orig)
adata_orig.write_h5ad(base_path+orig_path)
scgen_batch = adata_scgen.obs['batch_label'].values
orig_label =adata_orig.obs['label']
# orig_label = list(map(int, tl.text_label_to_number(orig_label)))
orig_batches = scgen_batch

# sc.pp.neighbors(adata_orig)
# sc.tl.umap(adata_orig)

# sc.pl.umap(adata_orig, color=['batch', 'celltype'])
print(adata_davae)
# sc.pp.neighbors(adata_davae)
# sc.tl.umap(adata_davae)
# # adata_davae.obs['batch']=adata_orig.obs['batch']
# adata_davae.obs['celltype']=adata_orig.obs['celltype']
# adata_davae.write_h5ad(base_path+'dann_vae/ifnb/davae_save01_umap.h5ad')
#
# sc.pp.neighbors(adata_scan)
# sc.tl.umap(adata_scan)
# adata_scan.obs['batch']=adata_orig.obs['batch']
# adata_scan.obs['celltype']=adata_orig.obs['celltype']
# adata_scan.write_h5ad(base_path+scan_path)
#
# sc.pp.neighbors(adata_desc)
# sc.tl.umap(adata_desc)
# adata_desc.obs['batch']=adata_orig.obs['batch']
# adata_desc.obs['celltype']=adata_orig.obs['celltype']
# adata_desc.write_h5ad(base_path+desc_path)
#
# print(adata_seurat)
# sc.pp.neighbors(adata_seurat)
# sc.tl.umap(adata_seurat)
# adata_seurat.obs['batch']=adata_orig.obs['batch']
# adata_seurat.obs['celltype']=adata_orig.obs['celltype']
# adata_seurat.write_h5ad(base_path+seurat_path)
#
# sc.pp.neighbors(adata_scgen, use_rep='corrected_latent')
# sc.tl.umap(adata_scgen)
# adata_scgen.write_h5ad(base_path+scgen_path)

#
# fig, (ax1, ax2, ax3, ax4, ax5 ,ax6) = plt.subplots(1, 6, figsize=(35,4), gridspec_kw={'wspace':0.9})
# ax1_dict = sc.pl.umap(adata_orig, color=['batch', 'celltype'], show=False)
# ax2_dict = sc.pl.umap(adata_davae, color=['batch', 'celltype'], show=False)
# ax3_dict = sc.pl.umap(adata_scan, color=['batch', 'celltype'], show=False)
# ax4_dict = sc.pl.umap(adata_desc, color=['batch', 'celltype'])



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
data_davae = adata_davae.X
data_desc = adata_desc.X
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

batch_cmap = 'tab20b'
c_map = 'tab20c'
size = 1
xy_label_size = 10
title_size = 5
xy_label = 'umap'


fig = plt.figure(figsize=(42, 12))

ax = fig.add_subplot(251)
plt.title('Raw', fontsize=title_size)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)

plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(256)
ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(252)
plt.title('DAVAE', fontsize=title_size)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(257)
ax.scatter(data_davae_emb[:, 0], data_davae_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(253)
plt.title('Scanorama', fontsize=title_size)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(258)
ax.scatter(data_scan_emb[:, 0], data_scan_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(254)
plt.title('DESC', fontsize=title_size)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2,5,9)
ax.scatter(data_desc_emb[:, 0], data_desc_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

# ---------------------------------------

ax = fig.add_subplot(255)
plt.title('Seurat 3.0', fontsize=title_size)
ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

ax = fig.add_subplot(2, 5, 10)
ax.scatter(data_seurat_emb[:, 0], data_seurat_emb[:, 1], c=orig_label, cmap=c_map, s=size,
           linewidth=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel(xy_label + '1', size=xy_label_size)
plt.ylabel(xy_label + '2', size=xy_label_size)

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
scgen_batch = adata_scgen.obs['batch_label'].values
# scgen_label = adata_scgen.obs['celltype'].values
# scgen_label = tl.text_label_to_number(scgen_label)
# print(scgen_label)
# print(orig_label)
# fig = plt.figure(figsize=(10, 12))
# ax = fig.add_subplot(346)
# plt.title('scGen', fontsize=title_size)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=scgen_batch, cmap=batch_cmap, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(3, 4, 12)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)

plt.show()
# plt.savefig(figure_merge)
plt.close(fig)

# ------------------------
# fig = plt.figure(figsize=(10, 12))
# ax = fig.add_subplot(211)
# plt.title('scGen', fontsize=title_size)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=scgen_batch, cmap=batch_cmap, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(212)
# ax.scatter(data_scgen_emb[:, 0], data_scgen_emb[:, 1], c=orig_label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# plt.xlabel(xy_label + '1', size=xy_label_size)
# plt.ylabel(xy_label + '2', size=xy_label_size)
# plt.show()
# plt.close(fig)
adata_seurat.obs['label']=adata_orig.obs['label']
# sc.pl.umap(adata_orig, color='label',cmap='tab20')
# sc.pl.umap(adata_seurat, color=['seurat_annotations', 'label'])
import collections
l=collections.Counter(adata_orig.obs['label'])
c=collections.Counter(adata_seurat.obs['seurat_annotations'])
print(l)
print(c)