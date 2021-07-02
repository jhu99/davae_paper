import scanpy as sc
import numpy as np
import umap
import ttools as tl
base_path = '/Users/zhongyuanke/data/'

# -----------------------------------293t-------------------------------------
# davae_path = base_path + 'dann_vae/pbmc/293t_save04_label.h5ad'
# # pbmc_path = base_path + 'dann_vae/pbmc/davae_save01.h5ad'
# orig_path = base_path + 'pbmc/zheng/293t_jurkat_merge.h5ad'
# scan_path = base_path + 'scanorama/scan_mcl.h5ad'
# seurat_path = base_path + 'seurat_result/293t.csv'
# desc_path = base_path + 'desc/desc_jurkat.h5ad'
# label_path = base_path + 'pbmc/zheng/293t_jurkat_cluster.txt'
# scgen_path = base_path+'scgen/scgen_mcl.h5ad'
#
# scgen_batch_path = base_path+'dann_vae/k_bet/293t/scgen_batch.csv'
# scgen_label_path = base_path+'dann_vae/k_bet/293t/scgen_label.csv'
#
# davae_kbet = base_path + 'dann_vae/k_bet/293t/davae.csv'
# orig_kbet = base_path + 'dann_vae/k_bet/293t/orig.csv'
# scan_kbet = base_path + 'dann_vae/k_bet/293t/scan.csv'
# desc_kbet = base_path + 'dann_vae/k_bet/293t/desc.csv'
# seurat_kbet = base_path + 'dann_vae/k_bet/293t/seurat.csv'
# scgen_kbet = base_path+'dann_vae/k_bet/293t/scgen.csv'

# -------------------------------smartseq-----------------------------------------
# davae_path = base_path + 'dann_vae/benchmark1/dc_davae_save.h5ad'
# orig_path = base_path + 'dann_vae/benchmark1/orig.h5ad'
# scan_path = base_path + 'scanorama/smart_seq.h5ad'
# seurat_path = base_path + 'seurat_result/smartseq.csv'
# desc_path = base_path + 'desc/desc_dc.h5ad'
# scgen_path = base_path+ 'scgen/scgen_dc.h5ad'
#
# davae_kbet = base_path + 'dann_vae/k_bet/smartseq/davae.csv'
# orig_kbet = base_path + 'dann_vae/k_bet/smartseq/orig.csv'
# scan_kbet = base_path + 'dann_vae/k_bet/smartseq/scan.csv'
# desc_kbet = base_path + 'dann_vae/k_bet/smartseq/desc.csv'
# seurat_kbet = base_path + 'dann_vae/k_bet/smartseq/seurat.csv'
# scgen_kbet=base_path + 'dann_vae/k_bet/smartseq/scgen.csv'
#
# label_path = base_path+'dann_vae/k_bet/smartseq/label.csv'
# batch_path = base_path+'dann_vae/k_bet/smartseq/batch.csv'

# # ----------------------------ifnb--------------------------------------------
# davae_path = base_path + 'dann_vae/ifnb/davae_save01.h5ad'
# orig_path = base_path + 'dann_vae/ifnb/orig.h5ad'
# scan_path = base_path + 'scanorama/ifnb.h5ad'
# seurat_path = base_path + 'seurat_result/ifnb.h5ad'
# desc_path = base_path + 'desc/ifnb.h5ad'
# scgen_path=base_path+'scgen/scgen_ifnb.h5ad'
#
davae_kbet = base_path + 'dann_vae/k_bet/ifnb/davae.csv'
orig_kbet = base_path + 'dann_vae/k_bet/ifnb/orig.csv'
scan_kbet = base_path + 'dann_vae/k_bet/ifnb/scan.csv'
desc_kbet = base_path + 'dann_vae/k_bet/ifnb/desc.csv'
seurat_kbet = base_path + 'dann_vae/k_bet/ifnb/seurat_new.csv'
label_path = base_path+'dann_vae/k_bet/ifnb/label.csv'
batch_path = base_path+'dann_vae/k_bet/ifnb/batch.csv'
scgen_kbet = base_path+'dann_vae/k_bet/ifnb/scgen.csv'

# ----------------------multimodal------------------------
# base_path = '/Users/zhongyuanke/data/'
orig_path = 'dann_vae/multimodal/multi_orig.h5ad'
davae_path = 'dann_vae/multimodal/davae_multi_save.h5ad'
scan_path = 'scanorama/scan_multimodal.h5ad'
seurat_path = 'seurat_result/seurat_multimodal.h5ad'
desc_path = 'desc/desc_multimodal.h5ad'
scgen_path = 'scgen/scgen_multimodal.h5ad'

adata_davae = sc.read_h5ad(davae_path)
adata_scan = sc.read_h5ad(scan_path)
adata_orig = sc.read_h5ad(orig_path)
adata_seurat = sc.read_h5ad(seurat_path)
adata_desc = sc.read_h5ad(desc_path)
adata_scgen = sc.read_h5ad(scgen_path)
# print(adata_seurat.obs['seurat_annotations'])
# print(adata_orig.obs['label'])
# print(adata_orig.obs['celltype'])

data_scan = adata_scan.X
data_davae = adata_davae.obsm['X_davae']
data_desc = adata_desc.X
# data_orig_emb = adata_orig.obsm['umap']
data_orig = adata_orig.X
data_seurat = adata_seurat.obsm['X_pca']
# data_seurat = adata_seurat.X
label = adata_orig.obs['label']
data_scgen = adata_scgen.obsm['corrected_latent']

data_orig_emb = umap.UMAP(n_components=6).fit_transform(data_orig)
data_seurat_emb = umap.UMAP(n_components=6).fit_transform(data_seurat)
data_scan_emb = umap.UMAP(n_components=6).fit_transform(data_scan)
data_desc_emb = umap.UMAP(n_components=6).fit_transform(data_desc)
data_scgen_emb = umap.UMAP(n_components=6).fit_transform(data_scgen)
data_davae_emb = umap.UMAP(n_components=6).fit_transform(data_davae)
# print(adata_scgen.obs['celltype'])
# print(adata_orig.obs['celltype'])
# print(data_orig_emb.shape)
# data_orig_emb = np.array(data_orig_emb)
# =-------batch=================================
# np.savetxt(label_path, label, delimiter=',')
# batch = adata_orig.obs['batch']
# batch = np.array(batch, dtype=float)
# print(batch)

# np.savetxt(batch_path, batch, delimiter=',')
# np.savetxt(scgen_kbet, data_scgen_emb, delimiter=',')
# scgen_batch = adata_scgen.obs['batch_label']
# np.savetxt(scgen_batch_path, scgen_batch, delimiter=',')
# scgen_label = adata_scgen.obs['celltype'].values
# scgen_label = tl.text_label_to_number(scgen_label)
# print(scgen_label)
np.savetxt(scgen_label_path, scgen_label, delimiter=',')
np.savetxt(davae_kbet, data_davae_emb, delimiter=',')
np.savetxt(orig_kbet, data_orig_emb, delimiter=',')
np.savetxt(scan_kbet, data_scan_emb, delimiter=',')
np.savetxt(desc_kbet, data_desc_emb, delimiter=',')
np.savetxt(seurat_kbet, data_seurat_emb, delimiter=',')
