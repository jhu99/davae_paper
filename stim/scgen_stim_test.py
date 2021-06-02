import scbean.tools.utils as tl
import scanpy as sc
import matplotlib
from numpy.random import seed
import scgen
seed(2021)
matplotlib.use('TkAgg')

adata = tl.read_sc_data('/Users/zhongyuanke/data/seurat_data/ifnb/ifnb.h5ad')
datasets = tl.split_object(adata, by="stim")
print(datasets[0])
print(datasets[1])
adata_all = tl.davae_preprocessing(datasets, n_top_genes=5000)
adata_all = scgen.setup_anndata(adata_all, batch_key="batch_label",  copy=True)
model = scgen.SCGEN(adata_all)
model.train(max_epochs=10,
    batch_size=32,
    early_stopping=True,
    early_stopping_patience=25,
    use_gpu=False)
corrected_adata = model.batch_removal()
corrected_adata.write('/Users/zhongyuanke/data/scgen/scgen_ifnb.h5ad')
