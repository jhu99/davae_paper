import scanpy as sc
import scbean.tools.utils as tl
import numpy as np
from scbean.model import davae as davae
from sklearn.metrics.pairwise import cosine_distances
import anndata
from sklearn.utils import shuffle
from keras.utils import to_categorical
import pandas as pd
import matplotlib
from numpy.random import seed
seed(2021)
matplotlib.use('TkAgg')


base_path = '/Users/zhongyuanke/data/'
anterior_out_path = 'dann_vae/spatial/rna_anterior_davae_01.h5ad'
posterior_out_path = 'dann_vae/spatial/rna_posterior_davae_01.h5ad'
file_rna = base_path+'spatial/mouse_brain/adata_processed_sc.h5ad'
rna_anterior_orig = base_path+'dann_vae/spatial/rna_anterior_orig.h5ad'

file1_spatial = base_path+'spatial/mouse_brain/10x_mouse_brain_Anterior/'
file2_spatial = base_path+'spatial/mouse_brain/10x_mouse_brain_Posterior/'
file1 = base_path+'spatial/mouse_brain/10x_mouse_brain_Anterior/V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5'
file2 = base_path+'spatial/mouse_brain/10x_mouse_brain_Posterior/V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5'
figure_umap = base_path+'dann_vae/spatial/umap.png'

adata_spatial_anterior = sc.read_visium(file1_spatial, count_file=file1)
adata_spatial_posterior = sc.read_visium(file2_spatial, count_file=file2)
adata_spatial_anterior.var_names_make_unique()
adata_spatial_posterior.var_names_make_unique()
adata_rna = sc.read_h5ad(file_rna)
# sc.pp.filter_genes(adata_rna, min_cells=500)
# sc.pp.highly_variable_genes(adata_rna, n_top_genes=5000)
# features = adata_rna.var_names[adata_rna.var['highly_variable']]
# adata_rna = adata_rna[:,features]
# adata_rna.write_h5ad(base_path+'spatial/mouse_brain/cortex_for_seurat.h5ad')
print(adata_rna)
print(adata_spatial_anterior)
print(adata_spatial_posterior)
adata_spatial_anterior = adata_spatial_anterior[
    adata_spatial_anterior.obsm["spatial"][:, 1] < 6000, :
]
adata_spatial_posterior = adata_spatial_posterior[
    (adata_spatial_posterior.obsm["spatial"][:, 1] < 4000)
    & (adata_spatial_posterior.obsm["spatial"][:, 0] < 6000),
    :,
]
adata_all = tl.spatial_rna_preprocessing(
        adata_spatial_anterior,
        adata_rna,
    )




def label_transfer(dist, labels):
    lab = pd.get_dummies(labels).to_numpy().T
    print(lab.shape)
    print(lab)
    class_prob = lab @ dist
    print(class_prob.shape)
    print(dist.shape)
    norm = np.linalg.norm(class_prob, 2, axis=0)
    print(norm.shape)
    class_prob = class_prob / norm
    class_prob = (class_prob.T - class_prob.min(1)) / class_prob.ptp(1)
    return class_prob


def integrate_spatial_rna(adata_spatial, adata_integrate, type='anterior'):
    sc.pp.neighbors(adata_integrate, use_rep='X_pca')
    # sc.tl.umap(adata_integrate)
    # sc.pl.umap(adata_integrate,color='batch')
    len_spatial = adata_spatial.shape[0]
    len_rna = adata_rna.shape[0]
    seurat_embedding = adata_integrate.X
    adata_spatial.obsm["seurat_embedding"] = seurat_embedding[0:len_spatial, :]
    adata_rna.obsm['seurat_embedding'] = seurat_embedding[len_spatial:len_rna + len_spatial, :]
    distances = 1 - cosine_distances(
        adata_rna.obsm["seurat_embedding"],
        adata_spatial.obsm['seurat_embedding']
    )
    class_prob_anterior = label_transfer(distances, adata_rna.obs.cell_subclass)
    cp_spatial_df = pd.DataFrame(
        class_prob_anterior, columns=np.sort(adata_rna.obs.cell_subclass.unique())
    )
    label = cp_spatial_df.idxmax(axis='columns').values
    cp_spatial_df.index = adata_spatial.obs.index
    adata_transfer = adata_spatial.copy()
    adata_transfer.obs = pd.concat(
        [adata_spatial.obs, cp_spatial_df], axis=1
    )
    sc.pl.spatial(
        adata_transfer,
        img_key="hires",
        color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
        size=1.5,
        color_map='Blues',
        ncols=2,
        legend_fontsize='xx-small'
    )

    adata_spatial.obs['celltype'] = label
    sc.pl.spatial(
        adata_spatial,
        img_key="hires",
        color='celltype',
        size=1.5,
        color_map='Set2'
    )
    label = list(label)
    from collections import Counter
    print(Counter(label))
    adata_spatial.write_h5ad('/Users/zhongyuanke/data/seurat_result/spatial/'+type+'_label_01.h5ad')

adata=sc.read_h5ad('/Users/zhongyuanke/data/seurat_result/spatial/posterior_rna.h5ad')
adata.obs['batch']=adata_all.obs['batch']
print(adata)
integrate_spatial_rna(adata_spatial_posterior, adata, type='posterior')
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(adata, color='batch', alpha=0.5)
#
# sc.pp.neighbors(adata_all)
# sc.tl.umap(adata_all)
# sc.pl.umap(adata_all, color='batch', alpha=0.5)