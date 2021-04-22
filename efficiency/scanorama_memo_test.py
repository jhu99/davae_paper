import scanorama
import scanpy as sc
import argparse
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
# adata = anndata.AnnData(X=adata.X.A,obs=adata.obs, var=adata.var)
print(adata1)
print(adata2)
t0 = time.time()
data_list = [adata1, adata2]

integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True, approx=False)
# integrated, corrected = scanorama.correct_scanpy(data_list, return_dimred=True)

t1 = time.time()
#
# info = psutil.virtual_memory()
# print('内存使用：', psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024, 'GB')
# print('总内存：', info.total/1024/1024/1024, 'GB')
# print('内存占比：', info.percent)
# print('cpu个数：', psutil.cpu_count())
t = (t1-t0)/60
print("Total time running: %s min" % (str(t)))

# result = [psutil.Process(os.getpid()).memory_info().rss/1024/1024/1024, t1-t0]
# result = np.array(result)
# np.savetxt(base_path+'scan_memo.txt', result)