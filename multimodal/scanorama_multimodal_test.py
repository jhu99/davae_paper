import scanpy as sc
import scbean.model.davae as davae
import scbean.tools.utils as tl
import scanorama
import matplotlib
import ttools as tool
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
file_rna = '/Users/zhongyuanke/data/dann_vae/multimodal/rna.h5ad'
file_atac = '/Users/zhongyuanke/data/dann_vae/multimodal/atac.h5ad'
seurat_celltype_path = base_path + 'multimodal/atac_pbmc_10k/celltype_filt.csv'
batch_size = 128

adata1 = sc.read_h5ad(file_atac)
adata2 = sc.read_h5ad(file_rna)
# adata_orig=adata1.concatenate(adata2)
# adata_orig.write_h5ad('/Users/zhongyuanke/data/dann_vae/multimodal/multi_orig.h5ad')
print(adata1)
print(adata2)
# sc.pp.filter_genes(adata1, min_cells=100)
# sc.pp.filter_genes(adata2, min_cells=100)
# sc.pp.log1p(adata1)
# sc.pp.log1p(adata2)
# sc.pp.scale(adata2)
# adata2.obs['celltype'] = adata1.obs['celltype']

# adata2.write_h5ad(base_path + 'multimodal/atac_pbmc_10k/activaty_matrix_label.h5ad')
adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000, hvg=False, lognorm=False)
# sc.pp.scale(adata_all)
datas = [adata1, adata2]

corrected = scanorama.correct_scanpy(datas, return_dimred=True, dimred=16)
print(corrected[0])

adata_corrected=corrected[0].concatenate(corrected[1])

print(adata_corrected)
sc.pp.neighbors(adata_corrected,use_rep='X_scanorama')
sc.tl.umap(adata_corrected)
# sc.pl.umap(adata_corrected,color='batch')
adata_corrected.write_h5ad('/Users/zhongyuanke/data/scanorama/scan_multimodal.h5ad')