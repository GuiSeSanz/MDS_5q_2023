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
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster.1=NULL, cluster.2=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster1 in  levels(all_seurat_integrated_sct[[res]][[1]])){
	for (cluster2 in  levels(all_seurat_integrated_sct[[res]][[1]])){
		if(cluster1!=cluster2){
			tmp <- FindMarkers(all_seurat_integrated_sct,
							ident.1 = cluster1,
							ident.2 = cluster2)
			tmp <- tmp[order(tmp$p_val_adj),]
			tmp$gene_name <- rownames(tmp)
			tmp$cluster.1 <- as.character(cluster1)
			tmp$cluster.2 <- as.character(cluster2)
			tmp <- tmp[, c('gene_name', 'pct.1', 'pct.2', 'cluster.1', 'cluster.2',  'p_val', 'avg_logFC', 'p_val_adj')]
			# tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_logFC) >2,]
			Results <- rbind(Results, tmp)
		}
	}
}
saveRDS(Results, paste0(PLOT_PATH,'Markers_1Vs1_', res, '_parallel.rds'))
Results_sp <- split(Results, Results$cluster.1)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'Markers_1Vs1_', res, '_parallel.xlsx'), SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


