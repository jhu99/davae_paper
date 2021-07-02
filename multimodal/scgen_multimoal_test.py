import scgen
import scanpy as sc
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
print(adata1)
print(adata2)
# adata_b1.obs_names_make_unique()
# adata_b2.obs_names_make_unique()
# adata_b3.obs_names_make_unique()

adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000, hvg=False, lognorm=False)
adata_all.obs_names_make_unique()

adata_all = scgen.setup_anndata(adata_all, batch_key="batch_label", copy=True)
model = scgen.SCGEN(adata_all)
model.train(max_epochs=15,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
    use_gpu=False)
corrected_adata = model.batch_removal()
sc.pp.neighbors(corrected_adata,use_rep='corrected_latent')
sc.tl.umap(corrected_adata)
sc.pl.umap(corrected_adata,color='batch')

corrected_adata.write('/Users/zhongyuanke/data/scgen/scgen_multimodal.h5ad')
