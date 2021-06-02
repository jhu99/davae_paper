import scanpy as sc
import scbean.model.davae as davae
import matplotlib
import ttools as tool
matplotlib.use('TkAgg')
import ttools as tl

base_path = '/Users/zhongyuanke/data/'
file_rna = '/Users/zhongyuanke/data/dann_vae/multimodal/rna.h5ad'
file_atac = '/Users/zhongyuanke/data/dann_vae/multimodal/atac.h5ad'
file_davae='/Users/zhongyuanke/data/dann_vae/multimodal/davae_multi_save.h5ad'
seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128
label=tl.get_label_by_txt(seurat_celltype_path)
label=label+label
adata = sc.read_h5ad(file_davae)
adata.obs['celltype'] = label
batch = ['0']*10412+['1']*10412
adata.obs['batch']=batch
sc.pl.umap(adata, color=['batch'],frameon=False, s=2)
sc.pl.umap(adata, color=['celltype'],frameon=False, s=2)
# sc.tl.leiden(adata)
# sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
# sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False)

# sc.pl.umap(adata, color=['LYN', 'ARHGAP26', 'PAX5', 'CD8A', 'ZEB2','LEF1','SULF2','CST7'],
#            s=1, ncols=4,frameon=False, cmap='viridis')
