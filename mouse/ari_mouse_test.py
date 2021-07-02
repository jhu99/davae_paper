import scanpy as sc
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder
from collections import Counter
import random
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--davae_path", type=str, default='davae_mouse.h5ad', help="save path")
opt = parser.parse_args()

def text_label_to_number(txt_label):
    encoder = LabelEncoder()
    label = encoder.fit_transform(txt_label)
    return label


cluster_count = 9
base_path = 'results/'
orig_path = 'davae_orig.h5ad'
davae_path = opt.davae_path
scan_path = 'scan_mouse.h5ad'
seurat_path = 'seurat_mouse.h5ad'
desc_path = 'desc_mouse.h5ad'
scgen_path = 'scgen_mouse.h5ad'


adata_orig = sc.read_h5ad(base_path+orig_path)
adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
# adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)

orig_batches = adata_davae.obs['batch_label'].values
orig_batches=list(map(int,orig_batches))

orig_label=adata_davae.obs['class'].values
orig_label=text_label_to_number(orig_label)
count=Counter(orig_batches)
print(count[0])
indexes_batch1=list(range(count[0]))
indexes_batch2=list(range(count[1]))
batch1_size=int(count[0]*0.05)
batch2_size=int(count[1]*0.05)
indexes_batch1 = random.sample(indexes_batch1, batch1_size)
indexes_batch2 = random.sample(indexes_batch2, batch2_size)
indexes_batch2=[i +count[0] for i in indexes_batch2]
indexes_sample=indexes_batch1+indexes_batch2

def filt_by_index(nums,idxs):
    res=[]
    for i in idxs:
        res.append(nums[i])
    return res

orig_batches=[0]*batch1_size+[1]*batch2_size

orig_label = filt_by_index(orig_label, indexes_sample)

print(Counter(orig_label))
print(Counter(orig_batches))
data_desc_emb = adata_desc.obsm['X_umap0.8']
data_orig_emb = adata_orig.obsm['X_umap']
data_davae_emb = adata_davae.obsm['X_umap']
#   = adata_seurat.obsm['X_umap']
data_scan_emb = adata_scan.obsm['X_umap']
data_scgen_emb=adata_scgen.obsm['X_umap']

data_orig_emb=filt_by_index(data_orig_emb, indexes_sample)
data_desc_emb=filt_by_index(data_desc_emb, indexes_sample)
data_davae_emb=filt_by_index(data_davae_emb, indexes_sample)
data_scan_emb=filt_by_index(data_scan_emb, indexes_sample)
data_scgen_emb=filt_by_index(data_scgen_emb, indexes_sample)

kmeans_orig = KMeans(n_clusters=cluster_count).fit(data_orig_emb)
# kmeans_orig = KMeans().fit(data_orig_emb)
ari_orig = adjusted_rand_score(orig_label, kmeans_orig.labels_)
ari_batch_orig = adjusted_rand_score(orig_batches, kmeans_orig.labels_)
print('original', ari_orig, ari_batch_orig)

kmeans_davae = KMeans(n_clusters=cluster_count).fit(data_davae_emb)
# kmeans_davae = KMeans().fit(data_davae_emb)
ari_davae = adjusted_rand_score(orig_label, kmeans_davae.labels_)
ari_batch_davae = adjusted_rand_score(orig_batches, kmeans_davae.labels_)
print('davae', ari_davae, ari_batch_davae)

kmeans_scan = KMeans(n_clusters=cluster_count).fit(data_scan_emb)
# kmeans_scan = KMeans().fit(data_scan_emb)
ari_scan = adjusted_rand_score(orig_label, kmeans_scan.labels_)
sh_scan = adjusted_rand_score(orig_batches, kmeans_scan.labels_)
print('scanorama', ari_scan, sh_scan)

kmeans_desc = KMeans(n_clusters=cluster_count).fit(data_desc_emb)
# kmeans_desc = KMeans().fit(data_desc_emb)
ari_desc = adjusted_rand_score(orig_label, kmeans_desc.labels_)
sh_desc = adjusted_rand_score(orig_batches, kmeans_desc.labels_)
print('desc', ari_desc, sh_desc)

# kmeans_seurat = KMeans(n_clusters=cluster_count).fit(data_seurat_emb)
# kmeans_seurat = KMeans().fit(data_seurat_emb)
# ari_seurat = adjusted_rand_score(orig_label, kmeans_seurat.labels_)
# sh_seurat = adjusted_rand_score(orig_batches, kmeans_seurat.labels_)
# print('seurat', ari_seurat, sh_seurat)

kmeans_scgen = KMeans(n_clusters=cluster_count).fit(data_scgen_emb)
# kmeans_scgen = KMeans().fit(data_scgen_emb)
ari_scgen = adjusted_rand_score(orig_label, kmeans_scgen.labels_)
sh_scgen = adjusted_rand_score(orig_batches, kmeans_scgen.labels_)
print('scGen', ari_scgen, sh_scgen)
