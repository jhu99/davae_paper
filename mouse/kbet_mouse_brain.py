import scanpy as sc
import numpy as np
import umap
import ttools as tl
import collections
from collections import Counter
import random

base_path = '/Users/zhongyuanke/data/'

orig_path = 'dann_vae/mouse/results/orig_mouse.h5ad'
davae_path = 'dann_vae/mouse/results/davae_mouse_save01.h5ad'
scan_path = 'dann_vae/mouse/results/scan_mouse.h5ad'
seurat_path = 'seurat_result/results/seurat_mouse.h5ad'
desc_path = 'dann_vae/mouse/results/desc_mouse.h5ad'
scgen_path = 'dann_vae/mouse/results/scgen_mouse.h5ad'

davae_kbet = base_path + 'dann_vae/k_bet/mouse/davae.csv'
orig_kbet = base_path + 'dann_vae/k_bet/mouse/orig.csv'
scan_kbet = base_path + 'dann_vae/k_bet/mouse/scan.csv'
desc_kbet = base_path + 'dann_vae/k_bet/mouse/desc.csv'
seurat_kbet = base_path + 'dann_vae/k_bet/mouse/seurat.csv'
scgen_kbet=base_path + 'dann_vae/k_bet/mouse/scgen.csv'

label_path='/Users/zhongyuanke/data/dann_vae/k_bet/mouse/label.csv'
batch_path='/Users/zhongyuanke/data/dann_vae/k_bet/mouse/batch.csv'

adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
# adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)

orig_batches = adata_davae.obs['batch'].values
orig_batches=list(map(int,orig_batches))

label=adata_davae.obs['class']
print(Counter(label))
print(Counter(orig_batches))
orig_label=tl.text_label_to_number(label)
print(Counter(orig_label))
count=Counter(orig_batches)

indexes_batch1=list(range(count[0]))
indexes_batch2=list(range(count[1]))
batch1_size=int(count[0]*0.2)
batch2_size=int(count[1]*0.1)

indexes_batch1 = random.sample(indexes_batch1, batch1_size)
indexes_batch2 = random.sample(indexes_batch2, batch2_size)
indexes_batch2=[i +count[0] for i in indexes_batch2]
indexes_sample=indexes_batch1+indexes_batch2

sample_path = '/Users/zhongyuanke/data/dann_vae/mouse/results/sample_index.txt'
np.savetxt(sample_path, indexes_sample)

# indexes_sample = np.loadtxt(sample_path,dtype='float')
# indexes_sample=list(map(int,indexes_sample))

def filt_by_index(nums,idxs):
    res=[]
    for i in idxs:
        res.append(nums[i])
    return res

orig_batches=[0]*batch1_size+[1]*batch2_size

orig_label = filt_by_index(orig_label, indexes_sample)
print(orig_label)
print(Counter(orig_label))
print(Counter(orig_batches))

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

data_orig_emb=filt_by_index(data_orig_emb,indexes_sample)
data_scgen_emb=filt_by_index(data_scgen_emb,indexes_sample)
data_davae_emb=filt_by_index(data_davae_emb,indexes_sample)
data_scan_emb=filt_by_index(data_scan_emb,indexes_sample)
data_desc_emb=filt_by_index(data_desc_emb,indexes_sample)

np.savetxt(batch_path, orig_batches, delimiter=',')
np.savetxt(label_path, orig_label, delimiter=',')
#
np.savetxt(davae_kbet, data_davae_emb, delimiter=',')
np.savetxt(orig_kbet, data_orig_emb, delimiter=',')
np.savetxt(scan_kbet, data_scan_emb, delimiter=',')
np.savetxt(desc_kbet, data_desc_emb, delimiter=',')
# np.savetxt(seurat_kbet, data_seurat_emb, delimiter=',')
np.savetxt(scgen_kbet, data_scgen_emb, delimiter=',')


