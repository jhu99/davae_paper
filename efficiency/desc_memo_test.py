import os
import scanpy as sc
import desc
import argparse
import anndata
import numpy as np
import time



parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/hca/', help="base path")
parser.add_argument("--type", type=str, default='1w', help="cell counts")

opt = parser.parse_args()

base_path = opt.base_path
file1 = base_path+'blood_' + opt.type + '.h5ad'
file2 = base_path+'bone_' + opt.type + '.h5ad'

adata1 = sc.read_h5ad(file1)
adata2 = sc.read_h5ad(file2)
print(adata1)
print(adata2)

adata = adata1.concatenate(adata2)
adata = anndata.AnnData(X=adata.X.A, obs=adata.obs, var=adata.var)
t0 = time.time()
# from concurrent.futures import ThreadPoolExecutor
#
# with ThreadPoolExecutor() as executor:
#     monitor = MemoryMonitor()
#     mem_thread = executor.submit(monitor.measure_usage)
#     try:
#         fn_thread = executor.submit(desc.train, adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
#                    batch_size=256, louvain_resolution=[0.8],
#                    save_dir="result", do_tsne=False, learning_rate=300,
#                    do_umap=False,
#                    save_encoder_weights=False)
#         result = fn_thread.result()
#
#     finally:
#         monitor.keep_measuring = False
#         max_usage = mem_thread.result()
#
#     print(f"Peak memory usage: {max_usage/1024/1024/1024} GB")
desc.normalize_per_cell(adata, counts_per_cell_after=1e4)

adata_out = desc.train(adata, dims=[adata.shape[1], 32, 16], tol=0.005, n_neighbors=10,
               batch_size=256, louvain_resolution=[0.8],
               save_dir="result", do_tsne=False, learning_rate=300,
               do_umap=False,
               save_encoder_weights=False)
t1 = time.time()
# info = psutil.virtual_memory()
# print('内存使用：', psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024, 'GB')
# print('总内存：', info.total/1024/1024/1024, 'GB')
# print('内存占比：', info.percent)
# print('cpu个数：', psutil.cpu_count())
t = (t1-t0)/60
print("Total time running: %s min" % (str(t)))
#
# result = [psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024, t1-t0]
# result = np.array(result)
# np.savetxt(base_path+'desc_memo.txt', result)