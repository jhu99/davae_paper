import numpy as np
import scbean.tools.utils as tl
from scbean.model import davae as davae
from keras.utils import to_categorical
import scanpy as sc
from sklearn.preprocessing import LabelEncoder
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from collections import Counter
import davae_test.spatial.deep_classifier as dc
import matplotlib
matplotlib.use('TkAgg')
from numpy.random import seed
seed(2021)

def deep_label_transfer(adata_spatial, adata_rna, type='anterior'):
    adata_all = tl.spatial_rna_preprocessing(adata_spatial, adata_rna)
    adata_integrate = davae.fit_integration(
        adata_all,
        epochs=45,
        hidden_layers=[128, 64, 32, 5],
        sparse=True,
        domain_lambda=3.0,
    )
    sc.pp.neighbors(adata_integrate, use_rep='X_davae')
    sc.tl.umap(adata_integrate)
    sc.pl.umap(adata_integrate, color='batch')
    rna_celltype = adata_rna.obs.cell_subclass
    print(rna_celltype)
    encoder = LabelEncoder()
    orig_label = encoder.fit_transform(rna_celltype)
    print(orig_label)
    orig_label.dtype = 'int64'

    davae_emb = adata_integrate.obsm['X_davae']
    len_spatial = adata_spatial.shape[0]
    len_rna = adata_rna.shape[0]
    test_set = davae_emb[0:len_spatial]
    train_set = davae_emb[len_spatial:len_spatial + len_rna]

    label = to_categorical(orig_label)
    print(label)
    class_num = label.shape[1]

    net_x = dc.CLASSIFIER(input_size=train_set.shape[1], class_num=class_num)
    net_x.build()
    net_x.compile()
    net_x.train(x=train_set, label=label, epochs=25, batch_size=128)
    pred_label = net_x.prediction(test_set)
    pred_label.dtype = 'int64'
    pred_type = encoder.inverse_transform(pred_label)

    # df = pd.DataFrame(pred_type)
    # df.to_csv('/Users/zhongyuanke/data/dann_vae/atac/pred_type_save03.csv')
    # np.savetxt('/Users/zhongyuanke/data/dann_vae/atac/pred_label_save03.csv', pred_label, delimiter=',')
    #
    # all_label = np.concatenate([pred_label, orig_label])
    # all_type = encoder.inverse_transform(all_label)

    print(pred_type)
    type_list = list(pred_type)
    print(Counter(type_list))

    adata_spatial.obs['celltype'] = pred_type
    # adata_davae.obs['cell type'] = all_type
    adata_spatial.write_h5ad(base_path + 'dann_vae/spatial/'+type+'_label_02.h5ad')


base_path = '/Users/zhongyuanke/data/'
file1 = base_path + 'spatial/mouse_brain/10x_mouse_brain_Anterior/' \
                    'V1_Mouse_Brain_Sagittal_Anterior_filtered_feature_bc_matrix.h5'
file2 = base_path + 'spatial/mouse_brain/10x_mouse_brain_Posterior/' \
                    'V1_Mouse_Brain_Sagittal_Posterior_filtered_feature_bc_matrix.h5'
file1_spatial = base_path+'spatial/mouse_brain/10x_mouse_brain_Anterior/'
file2_spatial = base_path+'spatial/mouse_brain/10x_mouse_brain_Posterior/'
rna_path = base_path+'spatial/mouse_brain/adata_processed_sc.h5ad'

adata1 = sc.read_visium(file1_spatial, count_file=file1)
adata2 = sc.read_visium(file2_spatial, count_file=file2)
adata_rna = sc.read_h5ad(rna_path)
adata1.var_names_make_unique()
adata2.var_names_make_unique()
adata1 = adata1[
    adata1.obsm["spatial"][:, 1] < 6000, :
]
adata2 = adata2[
    (adata2.obsm["spatial"][:, 1] < 4000)
    & (adata2.obsm["spatial"][:, 0] < 6000),
    :,
]

deep_label_transfer(adata2, adata_rna, type='posterior')
