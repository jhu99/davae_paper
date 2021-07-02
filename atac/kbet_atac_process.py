import scanpy as sc
import numpy as np
import umap
import ttools as tl
base_path = '/Users/zhongyuanke/data/'


davae_kbet = base_path + 'dann_vae/k_bet/atac/davae.csv'
seurat_kbet = base_path + 'dann_vae/k_bet/atac/seurat.csv'
batch_seurat_path = base_path+'dann_vae/k_bet/atac/batch_seurat.csv'
batch_davae_path = base_path+'dann_vae/k_bet/atac/batch_davae.csv'

davae_path = '/Users/zhongyuanke/data/dann_vae/atac/davae_save_final.h5ad'
seurat_path = '/Users/zhongyuanke/data/seurat_result/atac.h5ad'

adata_davae = sc.read_h5ad(davae_path)
adata_seurat = sc.read_h5ad(seurat_path)
# print(adata_seurat.obs['seurat_annotations'])
# print(adata_orig.obs['label'])
# print(adata_orig.obs['celltype'])

data_davae = adata_davae.obsm['davae']
# data_orig_emb = adata_orig.obsm['umap']
data_seurat = adata_seurat.obsm['X_pca']
# data_seurat = adata_seurat.X

data_seurat_emb = umap.UMAP(n_components=6).fit_transform(data_seurat)
data_davae_emb = umap.UMAP(n_components=6).fit_transform(data_davae)

# =-------batch=================================
# np.savetxt(label_path, label, delimiter=',')
batch=[0]*7064+[1]*(16496-7064)
np.savetxt(batch_seurat_path, batch, delimiter=',')
davae_batch=list(map(int,adata_davae.obs['batch'].values))
np.savetxt(batch_davae_path, davae_batch, delimiter=',')
# np.savetxt(scgen_kbet, data_scgen_emb, delimiter=',')
# scgen_batch = adata_scgen.obs['batch_label']
# np.savetxt(scgen_batch_path, scgen_batch, delimiter=',')
# scgen_label = adata_scgen.obs['celltype'].values
# scgen_label = tl.text_label_to_number(scgen_label)
# print(scgen_label)

np.savetxt(davae_kbet, data_davae_emb, delimiter=',')
np.savetxt(seurat_kbet, data_seurat_emb, delimiter=',')
