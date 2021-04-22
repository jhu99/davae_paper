import numpy as np
import time
import scanpy as sc
import desc
import argparse
import anndata

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/dann_vae/scalability/', help="base path")

opt = parser.parse_args()
base_path = opt.base_path

time_list = []
adata1 = sc.read_h5ad(base_path+'blood_5w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_5w.h5ad')
print(adata1)
print(adata2)

t0 = time.time()
adata = adata1.concatenate(adata2)
adata = anndata.AnnData(X=adata.X.A, obs=adata.obs, var=adata.var)
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)

adata_out = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
                   batch_size=256, louvain_resolution=[0.8],
                   save_dir="result", do_tsne=False, learning_rate=300,
                   do_umap=False,
                   save_encoder_weights=False)
t1 = time.time()
print("Total time running DAVAE 10w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)


adata1 = sc.read_h5ad(base_path+'blood_10w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_10w.h5ad')


t0 = time.time()
adata = adata1.concatenate(adata2)
adata = anndata.AnnData(X=adata.X.A,obs=adata.obs, var=adata.var)
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
adata_out = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
                   batch_size=256, louvain_resolution=[0.8],
                   save_dir="result", do_tsne=False, learning_rate=300,
                   do_umap=False,
                   save_encoder_weights=False)
t1 = time.time()
print("Total time running DAVAE 20w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_20w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_20w.h5ad')

t0 = time.time()
adata = adata1.concatenate(adata2)
adata = anndata.AnnData(X=adata.X.A,obs=adata.obs, var=adata.var)
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
adata_out = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
                   batch_size=256, louvain_resolution=[0.8],
                   save_dir="result", do_tsne=False, learning_rate=300,
                   do_umap=False,
                   save_encoder_weights=False)
t1 = time.time()
print("Total time running DAVAE 40w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_30w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_30w.h5ad')

t0 = time.time()
adata = adata1.concatenate(adata2)
adata = anndata.AnnData(X=adata.X.A, obs=adata.obs, var=adata.var)
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)
adata_out = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
                   batch_size=256, louvain_resolution=[0.8],
                   save_dir="result", do_tsne=False, learning_rate=300,
                   do_umap=False,
                   save_encoder_weights=False)
t1 = time.time()
print("Total time running DAVAE 60w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

time_list = np.array(time_list)
print('--------------DESC-------------------')
print(time_list)
np.savetxt(base_path+'desc_runtime.csv', time_list, delimiter=',')