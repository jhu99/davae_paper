import desc
import scanpy as sc
from scbean.tools import utils as tl


file1 = 'dropviz/mouse_brain_dropviz_filtered.h5ad'
file2 = 'nuclei/adata_nuclei_filtered.h5ad'
desc_path = 'results/desc_mouse.h5ad'
adata1 = sc.read_h5ad(file1)
adata2 = sc.read_h5ad(file2)

adata_all = tl.davae_preprocessing([adata1, adata2], n_top_genes=2000)

adata_all.X=adata_all.X.A
adata_out = desc.train(adata_all, dims=[adata_all.shape[1], 32, 16], tol=0.03, n_neighbors=10,
               batch_size=256, save_dir="result", do_tsne=False, learning_rate=300,
               do_umap=True,pretrain_epochs=50, louvain_resolution=0.8,
               save_encoder_weights=False,use_GPU=True)
adata_out.write_h5ad(desc_path)
