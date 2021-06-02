import ttools as tools
import anndata
import scbean.model.davae as davae
import scbean.tools.utils as tl
import scanpy as sc
import matplotlib
import numpy as np
import scanorama
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')

base_path = "/Users/zhongyuanke/data/vipcca/mixed_cell_lines/"
adata1 = sc.read_h5ad(base_path+"293t.h5ad")
adata2 = sc.read_h5ad(base_path+"jurkat.h5ad")
adata3 = sc.read_h5ad(base_path+"mixed.h5ad")
# adata1.obs_names_make_unique()
# adata2.obs_names_make_unique()
# adata3.obs_names_make_unique()
datas = [adata1, adata2, adata3]

corrected = scanorama.correct_scanpy(datas, return_dimred=True, dimred=16)
print(corrected)

adata_corrected=corrected[0]
adata_corrected.obs=datas[0].obs
for i in np.arange(1,len(corrected)):
    adata_i=corrected[i]
    adata_i.obs=datas[i].obs
    adata_corrected=adata_corrected.concatenate(adata_i,index_unique=None)




# adata_davae = sc.read_h5ad('/Users/zhongyuanke/data/dann_vae/panc8/davae_5.h5ad')
# adata_corrected = sc.read_h5ad('/Users/zhongyuanke/data/scanorama/scan_panc8.h5ad')
#
# adata_corrected.obs['batch_label']=adata_davae.obs['batch_label']
print(adata_corrected)
# print(adata_davae)
# sc.pp.neighbors(adata_corrected, use_rep='X_scanorama')
# sc.tl.umap(adata_corrected)
# sc.pl.umap(adata_corrected, color=['celltype', 'batch_label'], s=5)
adata_corrected.write_h5ad('/Users/zhongyuanke/data/scanorama/scan_mcl.h5ad')