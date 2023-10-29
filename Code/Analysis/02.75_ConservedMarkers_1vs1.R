library(Seurat)
library(argparse)
library(ggplot2)
library(future)
plan("multicore", workers = 64)



parser <- ArgumentParser(description='Get the DE genes between clusters, 1Vs1 on the resolution selected')
parser$add_argument('-R', '--Resolution', default="0.2", type="double",
                    help='Resolution of the clustering')
args <- parser$parse_args()
RESOLUTION <- args$Resolution


PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/')
dir.create(PLOT_PATH, showWarnings = FALSE)
all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/all_seurat_integrated_sct.rds'))
DefaultAssay(all_seurat_integrated_sct) <- 'SCT'
res <- paste0('integrated_snn_res.', RESOLUTION)
message(res)
Idents(object = all_seurat_integrated_sct) <- res
Results_sp <- list()
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster.1=NULL, cluster.2=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster1 in  levels(all_seurat_integrated_sct[[res]][[1]])){
	for (cluster2 in  levels(all_seurat_integrated_sct[[res]][[1]])){
		if(cluster1!=cluster2 & 
			all(table(all_seurat_integrated_sct$Sample, all_seurat_integrated_sct$integrated_snn_res.0.2)[,cluster2]>3) &
			all(table(all_seurat_integrated_sct$Sample, all_seurat_integrated_sct$integrated_snn_res.0.2)[,cluster1]>3)){
			message(paste0('cluster1: ', cluster1, ' cluster2: ', cluster2))
			tmp <- FindConservedMarkers(all_seurat_integrated_sct,
							ident.1 = cluster1,
							ident.2 = cluster2,
							grouping.var = "Sample")º
			tmp$gene_name <- rownames(tmp)
			tmp$Comparison <- paste0(as.character(cluster1), 'Vs', as.character(cluster2))
			Results <- rbind(Results, tmp)
		}
	}
}
saveRDS(Results, paste0(PLOT_PATH,'Conserved_Markers_1Vs1_', res, '_parallel.rds'))
Results_sp <- split(Results, Results$Comparison)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'Conserved_Markers_1Vs1_', res, '_parallel.xlsx'), SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

