import scbean.model.davae as davae
import scbean.tools.utils as tl
import scanpy as sc
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

base_path = "/Users/zhongyuanke/data/vipcca/mixed_cell_lines/"

adata_b1 = tl.read_sc_data(base_path+"293t.h5ad", batch_name="293t")
adata_b2 = tl.read_sc_data(base_path+"jurkat.h5ad", batch_name="jurkat")
adata_b3 = tl.read_sc_data(base_path+"mixed.h5ad", batch_name="mixed")

adata_all = tl.davae_preprocessing([adata_b1, adata_b2, adata_b3], n_top_genes=3000)
print(adata_all)
print(adata_all)
adata_integrate = davae.fit_integration(
    adata_all,
    batch_num=3,
    split_by='batch_label',
    domain_lambda=3.0,
    epochs=25,
    sparse=True,
    hidden_layers=[128, 64, 32, 5]
)
# sc.pp.neighbors(adata_integrate, use_rep='X_davae', n_neighbors=10)
# sc.tl.umap(adata_integrate)
import umap
adata_integrate.obsm['X_umap']=umap.UMAP().fit_transform(adata_integrate.obsm['X_davae'])
sc.pl.umap(adata_integrate, color=['_batch', 'celltype'], s=3)

# adata_integrate.write_h5ad('/Users/zhongyuanke/data/dann_vae/pbmc/davae_save02.h5ad')
