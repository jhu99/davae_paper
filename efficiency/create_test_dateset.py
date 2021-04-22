import scanpy as sc
import argparse
import anndata
parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/hca/', help="base path")

opt = parser.parse_args()

base_path = opt.base_path
file1 = base_path+'ica_cord_blood_h5.h5'
file2 = base_path+'ica_bone_marrow_h5.h5'

adata_blood = sc.read_10x_h5(file1)
sc.pp.filter_genes(adata_blood, min_cells=3000)
adata_bone = sc.read_10x_h5(file2)
sc.pp.filter_genes(adata_bone, min_cells=3000)
adata_blood.var_names_make_unique()
adata_bone.var_names_make_unique()
print(adata_blood.shape)
print(adata_bone.shape)

# adata_blood = anndata.AnnData(X=adata_blood.X.A, obs=adata_blood.obs, var=adata_blood.var)
# adata_bone = anndata.AnnData(X=adata_bone.X.A, obs=adata_bone.obs, var=adata_bone.var)
# adata_blood.write_h5ad(base_path+'blood_all.h5ad')
# adata_bone.write_h5ad(base_path+'bone_all.h5ad')
adata1 = adata_blood[0:10000, ]
adata2 = adata_bone[0:10000, ]
adata1.write_h5ad(base_path + 'blood_1w.h5ad')
adata2.write_h5ad(base_path + 'bone_1w.h5ad')


adata1 = adata_blood[0:50000, ]
adata2 = adata_bone[0:50000, ]
adata1.write_h5ad(base_path + 'blood_5w.h5ad')
adata2.write_h5ad(base_path + 'bone_5w.h5ad')

adata1 = adata_blood[0:100000, ]
adata2 = adata_bone[0:100000, ]
adata1.write_h5ad(base_path + 'blood_10w.h5ad')
adata2.write_h5ad(base_path + 'bone_10w.h5ad')
#
#
adata1 = adata_blood[0:200000, ]
adata2 = adata_bone[0:200000, ]
adata1.write_h5ad(base_path + 'blood_20w.h5ad')
adata2.write_h5ad(base_path + 'bone_20w.h5ad')
#
#
adata1 = adata_blood[0:300000, ]
adata2 = adata_bone[0:300000, ]
adata1.write_h5ad(base_path + 'blood_30w.h5ad')
adata2.write_h5ad(base_path + 'bone_30w.h5ad')
print('finish write')

