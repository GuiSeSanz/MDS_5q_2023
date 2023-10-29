import scanpy as sc
import plotnine as p9
import anndata as ad

import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway analysis

import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd



def harmonyze_gene_names(dataset1, dataset2):
    """Harmonyze gene names of two datasets"""
    # Get the intersection of genes
    genes1 = set(dataset1.var_names)
    genes2 = set(dataset2.var_names)
    genes = genes1.intersection(genes2)
    print(f'Number of genes in intersection: {len(genes)}')
    print(f'Number of genes in dataset1: {len(genes1)} and dataset2: {len(genes2)}')
    genes_1_sub = [gene.replace('.', '-') for gene in dataset1.var_names]
    genes_2_sub = [gene.replace('.', '-') for gene in dataset2.var_names]
    print(f'Number of genes in intersection after harmonization: {len(set(genes_1_sub).intersection(set(genes_2_sub)))}')
    # rename the gene names
    dataset1.var_names = genes_1_sub
    dataset1.var['features'] = genes_1_sub
    dataset2.var_names = genes_2_sub
    dataset2.var['features'] = genes_2_sub
    return dataset1, dataset2


plt.rcParams['figure.figsize'] = [10, 10]
plt.rcParams['figure.dpi'] = 100



adata_elder = sc.read(
    "/home/tereshkova/data/gserranos/MDS/Data/Elder_data_integrated_COUNTS.h5ad"
)
adata_elder.layers['counts'] = adata_elder.raw.X.copy()
adata_elder.layers['log1p'] = adata_elder.X.copy()
adata = sc.read(
    "/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_COUNTS.h5ad"
)
adata.layers['counts'] = adata.raw.X.copy()
adata.layers['log1p'] = adata.X.copy()





adata, adata_elder =harmonyze_gene_names(adata, adata_elder)

[i for i in adata.var['features'] if 'HAND2' in i]
[i for i in adata_elder.var_names if 'HAND2' in i]


all_data = adata.concatenate(adata_elder, batch_key='Project', batch_categories=['MDS', 'Elder'], join='inner')
all_data.obs['cell_5q'] = [ 'WT' if pd.isna(cat) else cat  for cat in all_data.obs['cell_5q'] ]


# all_data.obs["Comparison"] = ( all_data.obs["Sample"] + all_data.obs["cell_5q"])
all_data.obs["Comparison"] = [cell_type + '&' + cell for cell_type, cell in zip(all_data.obs["Cluster_names"], all_data.obs["cell_5q"])]
all_data.obs['Comparison'] = pd.Series(all_data.obs['Comparison'], dtype="category")
sample_key = 'Sample'
condition_key = 'Comparison'
groupby = 'Comparison'



# Plot cell-types for reference
plt.style.use('dark_background')
sc.pl.scatter(all_data, basis='umap', color='Sample', frameon=False, show=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/UMAP_pyhton_Sample.pdf", bbox_inches="tight")

plt.style.use('dark_background')
sc.pl.scatter(all_data, basis='umap', color='Cluster_names', frameon=False, show=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/UMAP_pyhton_Cluster_names.pdf", bbox_inches="tight")


sc.pp.filter_cells(all_data, min_genes=200)
sc.pp.filter_genes(all_data, min_cells=3)
# log1p normalize the data
sc.pp.normalize_total(all_data)
sc.pp.log1p(all_data)


plt.style.use('dark_background')
sc.pl.umap(all_data, color=[condition_key, groupby], frameon=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/UMAP_pyhton_Cluster_Test.pdf", bbox_inches="tight")



li.mt.rank_aggregate.by_sample(
    all_data,
    groupby=groupby,
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    verbose=True, # use 'full' to show all information
    n_perms=1000, # reduce number of permutations for speed
    return_all_lrs=True, # return all LR values
    )

all_data.uns["liana_res"].sort_values("magnitude_rank").head(10)



all_data.uns["liana_res"]

plt.rcParams['figure.figsize'] = [30, 25]

(li.pl.dotplot_by_sample(all_data, sample_key=sample_key,
                         colour="magnitude_rank", size="specificity_rank",
                         source_labels=list(set(all_data.obs['Cluster_names'])),
                         target_labels=list(set(all_data.obs['Cluster_names'])),
                         inverse_colour=True,
                         inverse_size=True,
                         ligand_complex=["HMGB1", "APP"],
                         receptor_complex=["CXCR4", "CD74", "RPSA"],
                         figure_size=(25, 20),
                         size_range=(0.5, 5),
                         ) +
 # rotate facet labels
 p9.theme(strip_text=p9.element_text(size=10, colour="black", angle=90))
 )
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Result_bySample.pdf", bbox_inches="tight")

all_data.uns["liana_res"].to_csv('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/ResultsPerSample_CT_and_Condition.csv', sep='\t', mode='a')

tensor = li.multi.to_tensor_c2c(all_data,
                                sample_key=sample_key,
                                score_key='magnitude_rank', # can be any score from liana
                                how='outer_cells' # how to join the samples
                                )


c2c.io.export_variable_with_pickle(tensor, "/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/tensor_MDS.pkl")


context_dict = all_data.obs[[sample_key, condition_key]].drop_duplicates()
context_dict = dict(zip(context_dict[sample_key], context_dict[condition_key]))
context_dict = defaultdict(lambda: 'Unknown', context_dict)

tensor_meta = c2c.tensor.generate_tensor_metadata(interaction_tensor=tensor,
                                                  metadata_dicts=[context_dict, None, None, None],
                                                  fill_with_order_elements=True
                                                  )
tensor = c2c.analysis.run_tensor_cell2cell_pipeline(tensor,
                                                    tensor_meta,
                                                    copy_tensor=True, # Whether to output a new tensor or modifying the original
                                                    rank=None, # Number of factors to perform the factorization. If None, it is automatically determined by an elbow analysis. Here, it was precomuputed.
                                                    tf_optimization='regular', # To define how robust we want the analysis to be.
                                                    random_state=0, # Random seed for reproducibility
                                                    device='cpu', # Device to use. If using GPU and PyTorch, use 'cuda'. For CPU use 'cpu'
                                                    elbow_metric='error', # Metric to use in the elbow analysis.
                                                    smooth_elbow=False, # Whether smoothing the metric of the elbow analysis.
                                                    upper_rank=20, # Max number of factors to try in the elbow analysis
                                                    tf_init='random', # Initialization method of the tensor factorization
                                                    tf_svd='numpy_svd', # Type of SVD to use if the initialization is 'svd'
                                                    cmaps=None, # Color palettes to use in color each of the dimensions. Must be a list of palettes.
                                                    sample_col='Element', # Columns containing the elements in the tensor metadata
                                                    group_col='Category', # Columns containing the major groups in the tensor metadata
                                                    output_fig=True, # Whether to output the figures. If False, figures won't be saved a files if a folder was passed in output_folder.
                                                    )



plt.style.use('ggplot')

factors, axes = c2c.plotting.tensor_factors_plot(interaction_tensor=tensor,
                                                 metadata = tensor_meta, # This is the metadata for each dimension
                                                 sample_col='Element',
                                                 group_col='Category',
                                                 meta_cmaps = ['viridis', 'Dark2_r', 'tab20', 'tab20'],
                                                 fontsize=10, # Font size of the figures generated
                                                 )
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Result_Tensors.pdf", bbox_inches="tight")


factors = tensor.factors
factors.keys()
factors['Contexts']
factors['Sender Cells']
factors['Ligand-Receptor Pairs']
factors['Sender Cells']
factors['Receiver Cells']

lr_loadings = factors['Ligand-Receptor Pairs']
lr_loadings.sort_values("Factor 5", ascending=False).head(10)