import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import argparse
from sklearn.preprocessing import normalize
import umap
from sklearn.manifold import TSNE
from sklearn.metrics.cluster import adjusted_rand_score, silhouette_score
from sklearn.cluster import KMeans
import ttools as tl


# ---------panc8---------------------
# cluster_count = 15
# base_path='/Users/zhongyuanke/data/'
# orig_path = '/seurat_data/panc_8/panc8.h5ad'
# desc_path = 'desc/desc_panc8.h5ad'
# davae_path = 'dann_vae/panc8/davae_5.h5ad'
# seurat_path = 'seurat_result/panc8_all.h5ad'
# scan_path = 'scanorama/scan_panc8_all.h5ad'
# scgen_path='scgen/scgen_panc8.h5ad'

# ----------------293t--------------------
cluster_count = 4
base_path='/Users/zhongyuanke/data/'
orig_path = 'dann_vae/benchmark1/orig.h5ad'
davae_path = 'dann_vae/benchmark1/dc_davae_temp.h5ad'
scan_path = 'scanorama/smart_seq.h5ad'
seurat_path = 'seurat_result/smartseq.csv'
desc_path = 'desc/desc_dc.h5ad'
scgen_path = 'scgen/scgen_dc.h5ad'

adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
adata_seurat = sc.read_csv(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)
print(adata_orig)

davae_label =adata_davae.obs['label']
davae_label = list(map(int, tl.text_label_to_number(davae_label)))
davae_batches = adata_davae.obs['batch_label'].values
batch_label=adata_davae.obs['batch_label']
sc.pp.neighbors(adata_seurat)
sc.tl.umap(adata_seurat)
# print(adata_scgen)
# sc.pp.neighbors(adata_orig)
# sc.tl.umap(adata_orig)
# sc.pp.neighbors(adata_davae,use_rep='X_davae')
# sc.tl.umap(adata_davae)
# sc.pp.neighbors(adata_scan, use_rep='X_scanorama')
# sc.tl.umap(adata_scan)
# sc.pp.neighbors(adata_desc)
# sc.tl.umap(adata_desc)
# sc.pp.neighbors(adata_seurat)
# sc.tl.umap(adata_seurat)
# sc.pp.neighbors(adata_scgen,use_rep='corrected_latent')
# sc.tl.umap(adata_scgen)
# print(adata_orig)
# adata_desc = sc.read_h5ad(base_path + desc_path)
sc.pp.neighbors(adata_desc, use_rep='X_Embeded_z0.8')
sc.tl.umap(adata_desc)

data_scan = adata_scan.X
# data_davae = adata_davae.obsm['davae']
data_davae = adata_davae.X
# data_desc = adata_desc.X
# data_desc_emb = adata_desc.obsm['X_umap0.8']
# data_orig_emb = adata_orig.obsm['umap']
data_orig = adata_orig.X
data_seurat = adata_seurat.X
data_scgen=adata_scgen.obsm['corrected_latent']

orig_label =adata_orig.obs['label']
# sc.pl.umap(adata_orig, color=['celltype'], cmap='tab20b')
orig_label = list(map(int, tl.text_label_to_number(orig_label)))
orig_batches = adata_orig.obs['batch'].values

# data_orig_emb = umap.UMAP().fit_transform(data_orig)
# data_davae_emb = umap.UMAP().fit_transform(data_davae)
# data_seurat_emb = umap.UMAP().fit_transform(data_seurat)
# data_scan_emb = umap.UMAP().fit_transform(data_scan)
data_desc_emb = adata_desc.obsm['X_umap']
data_orig_emb = adata_orig.obsm['X_umap']
data_davae_emb = adata_davae.obsm['X_umap']
data_seurat_emb = adata_seurat.obsm['X_umap']
data_scan_emb = adata_scan.obsm['X_umap']
data_scgen_emb=adata_scgen.obsm['X_umap']

print(data_scgen_emb.shape)

kmeans_orig = KMeans(n_clusters=cluster_count).fit(data_orig_emb)
# kmeans_orig = KMeans().fit(data_orig_emb)
ari_orig = adjusted_rand_score(orig_label, kmeans_orig.labels_)
sh_orig = silhouette_score(data_orig_emb, kmeans_orig.labels_)
print('original', ari_orig, sh_orig)

kmeans_davae = KMeans(n_clusters=cluster_count).fit(data_davae_emb)
# kmeans_davae = KMeans().fit(data_davae_emb)
ari_davae = adjusted_rand_score(davae_label, kmeans_davae.labels_)
sh_davae = silhouette_score(data_orig_emb, kmeans_davae.labels_)
print('davae', ari_davae, sh_davae)

kmeans_scan = KMeans(n_clusters=cluster_count).fit(data_scan_emb)
# kmeans_scan = KMeans().fit(data_scan_emb)
ari_scan = adjusted_rand_score(davae_label, kmeans_scan.labels_)
sh_scan = silhouette_score(data_scan_emb, kmeans_scan.labels_)
print('scanorama', ari_scan, sh_scan)

kmeans_desc = KMeans(n_clusters=cluster_count).fit(data_desc_emb)
# kmeans_desc = KMeans().fit(data_desc_emb)
ari_desc = adjusted_rand_score(davae_label, kmeans_desc.labels_)
sh_desc = silhouette_score(data_desc_emb, kmeans_desc.labels_)
print('desc', ari_desc, sh_desc)

kmeans_seurat = KMeans(n_clusters=cluster_count).fit(data_seurat_emb)
# kmeans_seurat = KMeans().fit(data_seurat_emb)
ari_seurat = adjusted_rand_score(davae_label, kmeans_seurat.labels_)
sh_seurat = silhouette_score(data_seurat_emb, kmeans_seurat.labels_)
print('seurat', ari_seurat, sh_seurat)

# --------------------------harmony----------------
# harmony_label = tl.get_label_by_txt('/Users/zhongyuanke/data/harmony_result/mcl_celltype.csv')
# print(len(harmony_label))
# kmeans_harmony = KMeans(n_clusters=cluster_count).fit(data_harmony_emb)
# ari_harmony = adjusted_rand_score(harmony_label, kmeans_harmony.labels_)
# sh_harmony = silhouette_score(data_harmony_emb, kmeans_harmony.labels_)
# print('harmony', ari_harmony, sh_harmony)
# -----------------------------------------------------------
print(data_scgen_emb.shape)
kmeans_scgen = KMeans(n_clusters=cluster_count).fit(data_scgen_emb)
# kmeans_scgen = KMeans().fit(data_scgen_emb)
ari_scgen = adjusted_rand_score(davae_label, kmeans_scgen.labels_)
sh_scgen = silhouette_score(data_scgen_emb, kmeans_scgen.labels_)
print('scgen', ari_scgen, sh_scgen)
