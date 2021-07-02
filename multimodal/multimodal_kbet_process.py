import scanpy as sc
import numpy as np
import umap
import ttools as tl
import collections

base_path = '/Users/zhongyuanke/data/'
orig_path = base_path + 'dann_vae/multimodal/multi_orig.h5ad'
davae_path = base_path + 'dann_vae/multimodal/davae_multi_save.h5ad'
scan_path = base_path + 'scanorama/scan_multimodal.h5ad'
seurat_path = base_path + 'seurat_result/seurat_multimodal.h5ad'
desc_path = base_path + 'desc/desc_multimodal.h5ad'
scgen_path = base_path + 'scgen/scgen_multimodal.h5ad'

davae_kbet = base_path + 'dann_vae/k_bet/multimodal/davae.csv'
orig_kbet = base_path + 'dann_vae/k_bet/multimodal/orig.csv'
scan_kbet = base_path + 'dann_vae/k_bet/multimodal/scan.csv'
desc_kbet = base_path + 'dann_vae/k_bet/multimodal/desc.csv'
seurat_kbet = base_path + 'dann_vae/k_bet/multimodal/seurat.csv'
scgen_kbet=base_path + 'dann_vae/k_bet/multimodal/scgen.csv'

label_path='/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/label.csv'
batch_path='/Users/zhongyuanke/data/dann_vae/k_bet/multimodal/batch.csv'
adata_davae = sc.read_h5ad(davae_path)
adata_scan = sc.read_h5ad(scan_path)
adata_orig = sc.read_h5ad(orig_path)
adata_seurat = sc.read_h5ad(seurat_path)
adata_desc = sc.read_h5ad(desc_path)
adata_scgen = sc.read_h5ad(scgen_path)

data_scan = adata_scan.obsm['X_scanorama']
data_davae = adata_davae.obsm['X_davae']
data_desc = adata_desc.obsm['X_Embeded_z0.8']
data_orig = adata_orig.X
data_seurat = adata_seurat.obsm['X_pca']
data_scgen = adata_scgen.obsm['corrected_latent']

# data_orig_emb = umap.UMAP(n_components=6).fit_transform(data_orig)
data_seurat_emb = umap.UMAP(n_components=6).fit_transform(data_seurat)
# data_scan_emb = umap.UMAP(n_components=6).fit_transform(data_scan)
# data_desc_emb = umap.UMAP(n_components=6).fit_transform(data_desc)
data_davae_emb = umap.UMAP(n_components=6).fit_transform(data_davae)
data_scgen_emb = umap.UMAP(n_components=6).fit_transform(data_scgen)

seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128
label = tl.get_label_by_txt(seurat_celltype_path)
label = label+label
print(collections.Counter(label))
orig_label = tl.text_label_to_number(label)
print(collections.Counter(orig_label))
orig_batches = adata_davae.obs['batch']
orig_batches=list(map(int, orig_batches))

# np.savetxt(batch_path, orig_batches, delimiter=',')
# np.savetxt(label_path, orig_label, delimiter=',')
#
np.savetxt(davae_kbet, data_davae_emb, delimiter=',')
# np.savetxt(orig_kbet, data_orig_emb, delimiter=',')
# np.savetxt(scan_kbet, data_scan_emb, delimiter=',')
# # np.savetxt(desc_kbet, data_desc_emb, delimiter=',')
# np.savetxt(seurat_kbet, data_seurat_emb, delimiter=',')
# np.savetxt(scgen_kbet, data_scgen_emb, delimiter=',')


