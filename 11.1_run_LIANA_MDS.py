import scanpy as sc
import plotnine as p9

import liana as li
import cell2cell as c2c
import decoupler as dc # needed for pathway analysis

import warnings
warnings.filterwarnings('ignore')
from collections import defaultdict
import matplotlib.pyplot as plt
from collections import Counter
import pandas as pd


plt.rcParams['figure.figsize'] = [10, 10]
plt.rcParams['figure.dpi'] = 100



adata = sc.read(
    "/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_COUNTS.h5ad"
)


# meta = pd.read_csv("/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_meta.csv", sep=',', index_col=0)
# umap = pd.read_csv("/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_umap.csv", sep='\t', index_col=0)

adata.layers['counts'] = adata.raw.X.copy()
adata.layers['log1p'] = adata.X.copy()
# adata.obsm['umap'] = umap.values


adata.obs["Comparison"] = ( adata.obs["Sample"] + "&" + adata.obs["cell_5q"] )

sample_key = 'Comparison'
condition_key = 'cell_5q'
groupby = 'Cluster_names'



# Plot cell-types for reference
plt.style.use('dark_background')
sc.pl.scatter(adata, basis='umap', color='Sample', frameon=False, show=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/UMAP_pyhton_Sample.pdf", bbox_inches="tight")

plt.style.use('dark_background')
sc.pl.scatter(adata, basis='umap', color='Cluster_names', frameon=False, show=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/UMAP_pyhton_Cluster_names.pdf", bbox_inches="tight")


sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
# log1p normalize the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)


plt.style.use('dark_background')
sc.pl.umap(adata, color=[condition_key, groupby], frameon=False)
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/UMAP_pyhton_Cluster_Test.pdf", bbox_inches="tight")



li.mt.rank_aggregate.by_sample(
    adata,
    groupby=groupby,
    sample_key=sample_key, # sample key by which we which to loop
    use_raw=False,
    verbose=True, # use 'full' to show all information
    n_perms=1000, # reduce number of permutations for speed
    return_all_lrs=True, # return all LR values
    )

adata.uns["liana_res"].sort_values("magnitude_rank").head(10)



adata.uns["liana_res"]

plt.rcParams['figure.figsize'] = [30, 25]

(li.pl.dotplot_by_sample(adata, sample_key=sample_key,
                         colour="magnitude_rank", size="specificity_rank",
                         source_labels=list(set(adata.obs['Cluster_names'])),
                         target_labels=list(set(adata.obs['Cluster_names'])),
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
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/Result_bySample.pdf", bbox_inches="tight")

tensor = li.multi.to_tensor_c2c(adata,
                                sample_key=sample_key,
                                score_key='magnitude_rank', # can be any score from liana
                                how='outer_cells' # how to join the samples
                                )

adata.uns["liana_res"].to_csv('/home/tereshkova/data/gserranos/MDS/Plots/Liana/ResultsPerSample.csv', sep='\t', mode='a')

c2c.io.export_variable_with_pickle(tensor, "/home/tereshkova/data/gserranos/MDS/Plots/Liana/tensor_MDS.pkl")


context_dict = adata.obs[[sample_key, condition_key]].drop_duplicates()
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
plt.savefig("/home/tereshkova/data/gserranos/MDS/Plots/Liana/Result_Tensors.pdf", bbox_inches="tight")


factors = tensor.factors
factors.keys()
factors['Contexts']
factors['Sender Cells']
factors['Ligand-Receptor Pairs']
factors['Sender Cells']
factors['Receiver Cells']

lr_loadings = factors['Ligand-Receptor Pairs']
lr_loadings.sort_values("Factor 5", ascending=False).head(10)