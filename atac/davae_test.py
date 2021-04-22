from scbean.tools import utils as tl
from scbean.model import davae as davae
import scanpy as sc
import matplotlib
matplotlib.use('TkAgg')

base_path = '/Users/zhongyuanke/data/'
out_path = 'dann_vae/atac/davae_02.h5ad'
file1 = base_path+'seurat_data/sc_atac/atac_v1_pbmc_10k_filtered_peak_bc_matrix_8000.csv'
file2 = base_path+'seurat_data/sc_atac/pbmc_10k_v3_8000.csv'

adata1 = sc.read_csv(file1)
adata2 = sc.read_csv(file2)

adata_all = tl.davae_preprocessing([adata1, adata2], hvg=False, lognorm=False, sparse=False)
adata_integrate = davae.fit_integration(adata_all, split_by='batch',
                                  domain_lambda=5,
                                  epochs=50,
                                  hidden_layers=[128, 64, 32, 5],
                                  sparse=False,
                                  mode='DACVAE'
                                        )
sc.pp.neighbors(adata_integrate, use_rep='X_davae')
sc.tl.umap(adata_integrate)
sc.pl.umap(adata_integrate, color='batch')
print(adata_integrate)
# adata_integrate.write_h5ad(base_path+out_path)
