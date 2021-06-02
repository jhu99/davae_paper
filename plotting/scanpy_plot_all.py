import scanpy as sc
import matplotlib.pyplot as plt
import ttools as tl
base_path='/Users/zhongyuanke/data/'
orig_path = 'pbmc/zheng/mcl_pre.h5ad'
desc_path = 'desc/desc_jurkat.h5ad'
davae_path = 'dann_vae/pbmc/293t_save04_label.h5ad'
seurat_path = 'seurat_result/mcl.h5ad'
scan_path = 'scanorama/scan_mcl.h5ad'
scgen_path = 'scgen/scgen_mcl.h5ad'

adata_davae = sc.read_h5ad(base_path+davae_path)
adata_scan = sc.read_h5ad(base_path+scan_path)
adata_orig = sc.read_h5ad(base_path+orig_path)
adata_seurat = sc.read_h5ad(base_path+seurat_path)
adata_scgen = sc.read_h5ad(base_path+scgen_path)
adata_desc = sc.read_h5ad(base_path+desc_path)
sc.pl.umap(adata_orig, color=['batch', 'celltype'])

# fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(35,4), gridspec_kw={'wspace':0.9})
# ax1_dict = sc.pl.umap(adata_orig, color=['batch'], show=False)
# ax2_dict = sc.pl.umap(adata_davae, color=['batch'], show=False)
# ax3_dict = sc.pl.umap(adata_scan, color=['batch'], show=False)
# ax4_dict = sc.pl.umap(adata_desc, color=['batch'])