import scanpy as sc
import argparse
import scbean.model.davae as davae
import scbean.tools.utils as tl

parser = argparse.ArgumentParser()
parser.add_argument("--base_path", type=str, default='/Users/zhongyuanke/data/hca/', help="base path")
parser.add_argument("--type", type=str, default='5w', help="cell counts")
parser.add_argument("--epoch", type=int, default=2, help="epochs")
opt = parser.parse_args()


base_path = opt.base_path
file1 = base_path+'blood_' + opt.type + '.h5ad'
file2 = base_path+'bone_' + opt.type + '.h5ad'


adata1 = sc.read_h5ad(file1)
adata2 = sc.read_h5ad(file2)
print(adata1)
print(adata2)

adata_all = tl.davae_preprocessing([adata1, adata2],lognorm=False, hvg=False)
adata_integrate=davae.fit_integration(
    adata_all,
    batch_num=2,
    split_by='batch_label',
    domain_lambda=2.0,
    epochs=opt.epoch,
    sparse=True,
    hidden_layers=[64, 32, 6]
)