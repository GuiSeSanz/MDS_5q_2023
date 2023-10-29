
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os


def uncompress_file(file_path):
    import gzip
    import shutil
    # Open the compressed file and uncompress it
    with gzip.open(file_path, 'rb') as f_in:
        # Extract the filename without the .gz extension
        uncompressed_file_path = file_path[:-3]
        # Open a new file to write the uncompressed data
        with open(uncompressed_file_path, 'wb') as f_out:
            # Copy the uncompressed data from the input file to the output file
            shutil.copyfileobj(f_in, f_out)

def saveresults(score, prediction, barcodes, path):
    import pandas as pd
    df = pd.DataFrame({
        'doublet_score': score,
        'predicted_doublet': prediction
    })
    df.to_csv(f'{path}/scrublet_output_table.csv', index=False)
    pd.DataFrame(list(zip(barcodes, score, prediction)), columns=['barcode', 'score', 'predicted_doublet']).to_csv(f'{path}/scrublet_output_table_cellID.csv', index=False)


# Set the graphical parameters
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = "Arial"
plt.rc("font", size=14)
plt.rcParams["pdf.fonttype"] = 42

sample_name = 'FS-0634-post'
# sample_name = 'FS-0406-post'
input_dir = f'/home/tereshkova/data/gserranos/MDS/Data/{sample_name}/outs/filtered_feature_bc_matrix'
PLOTS = True

print(f"Seeking doublets int: {input_dir}")

counts_matrix = scipy.io.mmread(os.path.join(input_dir, "matrix.mtx.gz")).T.tocsc()
if not os.path.exists(os.path.join(input_dir, "features.tsv")):
    uncompress_file(os.path.join(input_dir, "features.tsv.gz"))
genes = np.array(
    scr.load_genes(os.path.join(input_dir, "features.tsv"), delimiter="\t", column=1)
)

if not os.path.exists(os.path.join(input_dir, "barcodes.tsv")):
    uncompress_file(os.path.join(input_dir, "barcodes.tsv.gz"))

with open(os.path.join(input_dir, "barcodes.tsv"), "r") as f:
    barcodes = [x.strip() for x in f.readlines()]


print(
    "Counts matrix shape: {} rows, {} columns".format(
        counts_matrix.shape[0], counts_matrix.shape[1]
    )
)

print("Number of genes in gene list: {}".format(len(genes)))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

doublet_scores, predicted_doublets = scrub.scrub_doublets(
    min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
)

predicted_doublet_mask = scrub.call_doublets(threshold=0.2)
if 'predicted_doublet_mask' in locals() :
    saveresults(doublet_scores, predicted_doublet_mask, barcodes, f'/home/tereshkova/data/gserranos/MDS/Data/{sample_name}/outs/filtered_feature_bc_matrix')
else:
    saveresults(doublet_scores, predicted_doublets, barcodes, f'/home/tereshkova/data/gserranos/MDS/Data/{sample_name}/outs/filtered_feature_bc_matrix')



if PLOTS:
    PLOT_PATH = os.path.join(os.getcwd(), "Plots", "Doublets", sample_name)
    if not os.path.exists(PLOT_PATH):
        os.makedirs(PLOT_PATH)
    scrub.plot_histogram()
    plt.title("Histogram scores")
    plt.savefig(os.path.join(PLOT_PATH, "Hist.pdf"), bbox_inches="tight")
    plt.close()
    print("Running UMAP...")
    scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding("UMAP", order_points=True)
    plt.title("UMAP")
    plt.savefig(os.path.join(PLOT_PATH, "UMAP.png"), bbox_inches="tight")
    plt.close()
    # print("Running tSNE...")
    # scrub.set_embedding("tSNE", scr.get_tsne(scrub.manifold_obs_, angle=0.9))
    # scrub.plot_embedding("tSNE", order_points=True)
    # plt.title("tSNE")
    # plt.savefig(os.path.join(PLOT_PATH, "tSNE.png"), bbox_inches="tight")
    # plt.close()
    # print("Running ForceAtlas2...")
    # scrub.set_embedding(
    #     "FA", scr.get_force_layout(scrub.manifold_obs_, n_neighbors=5, n_iter=1000)
    # )
    # scrub.plot_embedding("FA", order_points=True)
    # plt.title("ForceAtlas2")
    # plt.savefig(os.path.join(PLOT_PATH, "FA.pdf"), bbox_inches="tight")
    # plt.close()
    print("Done.")
