import numpy as np
from scbean.tools import utils as tl
import scanpy as sc
import pandas as pd
import scanorama
import argparse

base_path = '/Users/zhongyuanke/data/'

file1 = 'dropviz/mouse_brain_dropviz_filtered.h5ad'
file2 = 'nuclei/adata_nuclei_filtered.h5ad'
scan_path = 'results/scan_mouse.h5ad'

# -------------train---------------------
adata1 = tl.read_sc_data(file1, fmt='h5ad')
adata2 = tl.read_sc_data(file2, fmt='h5ad')
# orig_label =adata_orig.obs['label']
print(adata1)
print(adata2)
datas = [adata1, adata2]
corrected = scanorama.correct_scanpy(datas, return_dimred=True, dimred=16)
adata_corrected=corrected[0].concatenate(corrected[1])

print(adata_corrected)
sc.pp.neighbors(adata_corrected,use_rep='X_scanorama')
sc.tl.umap(adata_corrected)
adata_corrected.write_h5ad(scan_path)






