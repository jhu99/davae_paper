import scanpy as sc
import matplotlib.pyplot as plt

base_path = '/Users/zhongyuanke/data/'
file_davae = base_path+'dann_vae/atac/davae_save_final.h5ad'
adata_davae = sc.read_h5ad(file_davae)
adata_seurat=sc.read_h5ad('/Users/zhongyuanke/data/seurat_result/atac.h5ad')
print(adata_seurat)
print(adata_davae)
batch=[0]*7064+[1]*(16496-7064)
sc.pp.neighbors(adata_seurat, use_rep='X_pca')
sc.tl.umap(adata_seurat)
adata_seurat.obs['batch']=batch
sc.pl.umap(adata_seurat, color=['batch'],cmap='Dark2')

# data_davae_emb = adata_davae.obsm['X_umap']
# data_seurat_emb = adata_seurat.obsm['X_umap']
#
#
# batch_cmap = 'Dark2'
# c_map = 'tab20b'
# size = 1
# xy_label_size = 10
# title_size = 10
# xy_label = 'umap'
#
# #
# fig = plt.figure(figsize=(40, 20))
#
# ax = fig.add_subplot(121)
# plt.title('Raw', fontsize=title_size)
# ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_batches, cmap=batch_cmap, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# # ax.spines['bottom'].set_visible(False)
# # ax.spines['left'].set_visible(False)
# #
# # plt.xlabel(xy_label + '1', size=xy_label_size)
# # plt.ylabel(xy_label + '2', size=xy_label_size)
#
# ax = fig.add_subplot(122)
# ax.scatter(data_orig_emb[:, 0], data_orig_emb[:, 1], c=orig_label, cmap=c_map, s=size,
#            linewidth=0)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)