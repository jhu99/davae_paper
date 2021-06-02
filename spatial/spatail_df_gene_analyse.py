import scanpy as sc


# davae_path = '/Users/zhongyuanke/data/seurat_result/spatial/anterior_label_01.h5ad'
davae_path = '/Users/zhongyuanke/data/dann_vae/spatial/posterior_label_01.h5ad'

adata_davae = sc.read_h5ad(davae_path)
label = adata_davae.obs['celltype']
label = list(label)
from collections import Counter
print(Counter(label))
label = ['nan' if i == 'Endo' or i == 'L6b' or i == 'CR' else i for i in label] # for davae
# label = ['nan' if i == 'Macrophage' or i == 'Lamp5' else i for i in label] #for seurat
adata_davae.obs['celltype'] = label
print(Counter(label))
# sc.pl.spatial(
#     adata_davae,
#     img_key="hires",
#     color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
#     size=1.5,
#     color_map='Blues'
# )

# sc.pl.spatial(
#         adata_davae,
#         img_key="hires",
#         color='celltype',
#         size=1.5,
#         color_map='Set2'
#     )
# sc.pp.filter_genes(adata_davae, min_cells=80)
# sc.pp.m
sc.pp.neighbors(adata_davae)
sc.tl.louvain(adata_davae)
# sc.tl.rank_genes_groups(adata_davae, groupby='celltype', n_genes=adata_davae.shape[1], method='wilcoxon')
# sc.pl.rank_genes_groups_dotplot(adata_davae,
#                         groups=['L2/3 IT', 'L4', 'L5 IT', 'L6 CT',
#                                 'Vip', 'Lamp5', 'Pvalb', 'Astro', 'Sst',
#                                 'Oligo'],
#                                 n_genes=5)
#
# sc.pl.rank_genes_groups_dotplot(adata_davae,
#                         groups=['L2/3 IT', 'L4', 'L5 IT', 'L6 CT',
#                                 'Vip', 'Pvalb', ],
#                                 n_genes=10)
# sc.tl.umap(adata_davae)
# sc.pl.umap(adata_davae, color=['celltype'])
marker_genes_dict_hca = {
    'Astro': ['Ptgds', 'Mgp', 'Fn1'],
    'L2/3 IT': ['Lamp5', 'Rasgrf2', 'Epop'],
    'L4': ['Lamp5', 'Nptxr', 'Cux2', 'Camk2a'],
    'L5 IT': ['Snap25', 'Scn1b', 'Cck'],
    'L5 PT': ['Nxph3', 'Ttc9b'],
    'L6 CT': ['Ccn2', 'Cplx3'],
    'Lamp5': ['Sez6', 'Gad2', 'Pcp4'],
    'Meis2': ['Tmem243', 'Pak4'],

    'Oligo': ['Plp1', 'Trf', 'Mal'],
    'Pvalb': ['Doc2g', 'S100a5','Mobp'],
    'Sst': ['Apod', 'Atp1a2', 'Apoe'],
    'Vip': ['Penk', 'Pde10a'],
    'test':['Gad1','Gad2','Slc32a1','Slc17a7','Lamp5','Ndnf','Sncg',
            'Vip','Sst','Pvalb','Cux2','Rorb','Fezf2','Sulf1','Foxp2','Nxph4','Aqp4']
                         }

# sc.pl.heatmap(adata_davae, marker_genes_dict_hca, groupby='celltype', dendrogram=False)
geneset=['Necab2','lah1','Nrgn','Cck','Sst','Pvalb','Pcp4','Penk','Plp1','Mbp','Pcbp3','Gng4']
# sc.pl.heatmap(adata_davae, geneset, groupby='celltype', dendrogram=False)
# sc.pl.heatmap(adata_davae, geneset, groupby='celltype', dendrogram=False)
base_path = '/Users/zhongyuanke/data/'
rna_path = base_path+'spatial/mouse_brain/adata_processed_sc.h5ad'
adata_rna = sc.read_h5ad(rna_path)
# sc.tl.rank_genes_groups(adata_rna, groupby='cell_subclass', n_genes=adata_rna.shape[1], method='wilcoxon')
# sc.pl.rank_genes_groups(adata_rna,
#                         groups=['L2/3 IT','L4','L5 IT','L5 PT','L6 CT'],
#                                 n_genes=10)
print(adata_rna.var)
# seurat marker-------------------------------
sc.pl.spatial(
    adata_davae,
    img_key="hires",
    # color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
    color=['Cd4','Ido1',],
    size=1.5,
    color_map='Blues',
    ncols=2,
    legend_fontsize='xx-small'
)

# ---------------------DAVAE marker-------------------------------
sc.pl.spatial(
    adata_davae,
    img_key="hires",
    # color=["L2/3 IT", "L4", "L5 PT", "L6 CT"],
    color=['Ighm','Rprm'],
    size=1.5,
    color_map='Blues',
    ncols=2,
    legend_fontsize='xx-small'
)
#sc.pl.heatmap(adata_rna, marker_genes_dict_hca, groupby='cell_subclass', dendrogram=False)
# file_rna = base_path+'spatial/mouse_brain/adata_processed_sc.h5ad'
# adata_rna = sc.read_h5ad(file_rna)
# print(adata_rna.obs['cell_subclass'])
# sc.pp.filter_genes(adata_rna, min_cells=1)
# sc.tl.rank_genes_groups(adata_rna, groupby='cell_subclass', n_genes=adata_davae.shape[1], method='wilcoxon')
# sc.pl.rank_genes_groups_dotplot(adata_rna,
#                         groups=['L2/3 IT', 'L4', 'L5 IT', 'L6 CT',
#                                 'Vip', 'Lamp5', 'Pvalb', 'Astro', 'Sst',
#                                 'Oligo'],
#                                 n_genes=5)