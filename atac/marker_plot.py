import scanpy as sc

base_path = '/Users/zhongyuanke/data/'
file_davae = base_path+'dann_vae/atac/davae_save_final.h5ad'
adata_davae = sc.read_h5ad(file_davae)

sc.pl.umap(adata_davae, color=['LYN', 'VCAN', 'CCL5', 'BANK1', 'SULF2', 'LDHB'],
            s=1, frameon=False, ncols=3, vmax='p99')