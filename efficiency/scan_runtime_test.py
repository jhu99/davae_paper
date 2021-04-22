import numpy as np
import time
import scanpy as sc
import scanorama
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/dann_vae/scalability/', help="base path")

opt = parser.parse_args()

base_path = opt.base_path
file1 = base_path+'blood_5w.h5ad'
file2 = base_path+'bone_5w.h5ad'

time_list = []
adata1 = sc.read_h5ad(base_path+'blood_5w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_5w.h5ad')
print(adata1)
print(adata2)
data_list = [adata1, adata2]
t0 = time.time()
integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True)
t1 = time.time()
print("Total time running DAVAE 10w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)


adata1 = sc.read_h5ad(base_path+'blood_10w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_10w.h5ad')

data_list = [adata1, adata2]
t0 = time.time()
integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True)
t1 = time.time()
print("Total time running DAVAE 20w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_20w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_20w.h5ad')
data_list = [adata1, adata2]
t0 = time.time()
integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True)
t1 = time.time()
print("Total time running DAVAE 40w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

adata1 = sc.read_h5ad(base_path+'blood_30w.h5ad')
adata2 = sc.read_h5ad(base_path+'bone_30w.h5ad')
data_list = [adata1, adata2]
t0 = time.time()
integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True)
t1 = time.time()
print("Total time running DAVAE 60w cells: %s seconds" % (str(t1-t0)))
time_list.append(t1-t0)

time_list = np.array(time_list)
print('--------------Scanorama-------------------')
print(time_list)
np.savetxt(base_path+'scan_runtime.csv', time_list, delimiter=',')