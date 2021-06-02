from numpy.random import seed
seed(2021)
import scanpy as sc
import ttools as tll
import scbean.tools.utils as tl
import matplotlib
matplotlib.use('TkAgg')
base_path='/Users/zhongyuanke/data/'
orig_path = 'dann_vae/ifnb/orig.h5ad'
davae_path = base_path+'dann_vae/ifnb/davae_save01.h5ad'
# davae_path = base_path+'dann_vae/ifnb/davae_01.h5ad'

adata = sc.read_h5ad(davae_path)
print(adata)
tech1 = 'ctrl'
tech2 = 'stim'
type = '_8000'
# sc.pp.neighbors(adata,use_rep='X_davae')
# sc.tl.umap(adata)
data1, celltype1 = tll.get_ifnb(tech1, type)
data2, celltype2 = tll.get_ifnb(tech2, type)
orig_len_x = data1.shape[0]
orig_len_y = data2.shape[0]
adata1 = adata[0:orig_len_x, ]
adata2 = adata[orig_len_x:, ]
adata1.obs['celltype']=celltype1
adata2.obs['celltype']=celltype2

# sc.pp.scale(adata1)
# sc.pp.scale(adata2)

# adata = tl.read_sc_data('/Users/zhongyuanke/data/seurat_data/ifnb/ifnb.h5ad')
# datasets = tl.split_object(adata, by="stim")
# print(datasets[0])
# print(datasets[1])
# adata1=datasets[0]
# adata2=datasets[1]
# print(adata1.var['features'])
# sc.pl.umap(adata, color=['CD74', 'FTL', 'BTG1', 'NKG7', 'FCGR3A', 'PABPC1',
#                          'SEC61B', 'GIMAP7', 'HBB', 'CXCR4', 'HLA-DPB1', 'NPM1', 'CCL2'],
#            s=2, frameon=False, ncols=4, vmax='p99')
# sc.pl.umap(adata, color=['batch', 'celltype'],
#            s=2, frameon=False, ncols=4, vmax='p99')
# sc.pl.umap(adata, color=['CD14'],
#            s=2, frameon=False, ncols=4, vmax='p99')




# sc.pl.umap(adata2, color=["CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", "CCL2", "PPBP"],
#           s=3, frameon=False, ncols=3, vmax='p99')
# sc.pl.umap(adata2, color=['CD69', 'TNFSF10', 'APOBEC3A', "CCL8",  "SELL", 'CXCL10',  'ISG15', 'LAG3'],
#           s=3, frameon=False, ncols=3, vmax='p99')
# sc.pl.umap(adata1, color=['CD69', 'TNFSF10', 'APOBEC3A', "CCL8",  "SELL", 'CXCL10','ISG15', 'LAG3'],
#           s=3, frameon=False, ncols=3, vmax='p99')
# sc.pl.stacked_violin(adata1, ['CD69', 'TNFSF10', 'APOBEC3A', "CCL8",  "SELL", 'CXCL10',  'ISG15', 'LAG3'],
#                      groupby='cell type2', swap_axes=True)
# sc.pl.stacked_violin(adata2, ['CD69', 'TNFSF10', 'APOBEC3A', "CCL8",  "SELL", 'CXCL10',  'ISG15', 'LAG3'],
#                      groupby='cell type2', swap_axes=True)

sc.pl.umap(adata2, color=['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10',  'ISG15', 'LAG3'],
          s=2, frameon=False, ncols=6, vmax=4, vmin=0, cmap='viridis',
           legend_fontsize=4)
sc.pl.umap(adata1, color=['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10', 'ISG15', 'LAG3'],
          s=2, frameon=False, ncols=6, vmax=4, vmin=0, cmap='viridis',
           legend_fontsize=4)
sc.pl.stacked_violin(adata1, ['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10',  'ISG15', 'LAG3'],
                     groupby='celltype', swap_axes=True, show=True, standard_scale='var',
                     vmax=1, vmin=0)
sc.pl.stacked_violin(adata2, ['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10',  'ISG15', 'LAG3'],
                     groupby='celltype', swap_axes=True, show=True, standard_scale='var',
                     vmax=1, vmin=0)

# marker_genes_dict = {
#     'shared': ['TNFSF10'],
#     'T cell':['LAG3'],
# }
# sc.pl.matrixplot(adata1, marker_genes_dict, 'seurat_annotations', dendrogram=True,)
# # sc.pl.matrixplot(adata1, ['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10',  'ISG15', 'LAG3'], groupby='seurat_annotations',)
#
# sc.pl.matrixplot(adata2, ['TNFSF10', 'APOBEC3A', "CCL8", 'CXCL10',  'ISG15', 'LAG3'], groupby='seurat_annotations',
#                    cmap='viridis', dendrogram=False)