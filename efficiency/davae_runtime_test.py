import scbean.model.davae as davae
import scbean.tools.utils as tl
import numpy as np
import time
import scanpy as sc
import os
import psutil
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/hca/', help="base path")
parser.add_argument("--epoch", type=int, default=15, help="epochs")
opt = parser.parse_args()

base_path = opt.base_path
epoch=opt.epoch
#
time_list = []
adata1 = sc.read_h5ad(base_path+'blood_5w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_5w.h5ad')
print(adata1)
print(adata2)
adata_all = tl.davae_preprocessing([adata1,adata2])
t0 = time.time()
adata_out = davae.fit_integration(adata_all, batch_size=256, epochs=epoch, sparse=True)
t1 = time.time()
print("Total time running DAVAE 10w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)
info = psutil.virtual_memory()
print('内存使用：', psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024, 'GB')
print('总内存：', info.total/1024/1024/1024, 'GB')
print('内存占比：', info.percent)
print('cpu个数：', psutil.cpu_count())


adata1 = sc.read_h5ad(base_path+'blood_10w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_10w.h5ad')
print(adata1)
print(adata2)
adata_all = tl.davae_preprocessing([adata1,adata2])
t0 = time.time()
adata_out = davae.fit_integration(adata_all, batch_size=3000, epochs=opt.epoch, sparse=True)
t1 = time.time()
print("Total time running DAVAE 20w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_20w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_20w.h5ad')
data_list = [adata1, adata2]
t0 = time.time()
adata_out = davae.fit_integration(data_list, batch_size=3000, epochs=opt.epoch, sparse=True)
t1 = time.time()
print("Total time running DAVAE 40w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_30w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_30w.h5ad')
data_list = [adata1, adata2]
t0 = time.time()
adata_out = davae.fit_integration(data_list, batch_size=3000, epochs=opt.epoch, sparse=True)
t1 = time.time()
print("Total time running: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)
time_list = np.array(time_list)
print('---------------DAVAE-------------------')
print(time_list)
np.savetxt(base_path+'davae_runtime.csv', time_list, delimiter=',')