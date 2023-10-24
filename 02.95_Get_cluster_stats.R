library(Seurat)
library(argparse)
library(ggplot2)
library(future)
plan("multicore", workers = 64)





plot_Stat_umap <- function(coords, resolution, stat){
	resolution <- setNames(as.data.frame(all_seurat_integrated_sct[[resolution]]), 'Cluster')
	stat_df <- setNames(as.data.frame(all_seurat_integrated_sct[[stat]]), c('Value'))
	plotter <- merge(coords, resolution, by=0)
	plotter <- merge(plotter, stat_df, by.x= 'Row.names', by.y=0)
	p <- ggplot(plotter, aes(x=UMAP_1, y=UMAP_2, fill = Value)) + 
		geom_point(data = coords[sample(nrow(coords), 1e3),], fill='grey', color='grey', alpha=0.2, size=1) + 
		geom_point(size=1, alpha=1, pch=21) + theme_classic() + 
		ggtitle(stat) + 
		theme(legend.position = "bottom", axis.ticks=element_blank(), axis.text=element_blank()) +
		facet_wrap(~Cluster) 
	if (is.numeric(plotter$Value[1])){
		p <- p + viridis::scale_fill_viridis()
	}else{
		colors <- ggthemes::tableau_color_pal('Classic 10')(length(unique(plotter$Value)))
		p <- p + scale_fill_manual(values = c(colors))
	}
	return(p)
}

get_stats_table <- function(resolution, stat_list){
	resolution <- setNames(as.data.frame(all_seurat_integrated_sct[[resolution]]), 'Cluster')
	table_df <- data.frame(Cluster = unique(resolution$Cluster))
	for(stat in stat_list){
		tmp <- all_seurat_integrated_sct[[stat]]
		tmp <- merge(resolution, tmp, by=0)
		tmp <- setNames(aggregate(tmp[, c(stat)], list(tmp$Cluster), mean), c('Cluster', stat))
		table_df <- merge(table_df, tmp , by='Cluster')
	}
	ncells <- setNames(as.data.frame(table(resolution$Cluster)), c('Cluster', 'Ncells'))
	table_df <- merge(table_df, ncells, by='Cluster')
	return(table_df)
}


PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/Stats_per_cluster/')
dir.create(PLOT_PATH, showWarnings = FALSE)
all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/all_seurat_integrated_sct.rds'))

DefaultAssay(all_seurat_integrated_sct) <- "integrated"
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)



for (res in (paste0('integrated_snn_res.', seq(0.2,1.2, by=0.2)))){
	message(res)
	pdf(paste0(PLOT_PATH, 'Stats_per_cluster',res,'.pdf'))
		print(plot_Stat_umap(coords, res, 'Phase'))
		print(plot_Stat_umap(coords, res, 'S.Score'))
		print(plot_Stat_umap(coords, res, 'G2M.Score'))
		print(plot_Stat_umap(coords, res, 'RPSRatio'))
		print(plot_Stat_umap(coords, res, 'mitoRatio'))
		print(plot_Stat_umap(coords, res, 'nUMI'))
	dev.off()
	stats <- get_stats_table(res, c('S.Score', 'G2M.Score', 'RPSRatio', 'mitoRatio', 'nUMI'))
	WriteXLS::WriteXLS(stats, ExcelFileName=paste0(PLOT_PATH, 'Stats_per_cluster',res,'.xlsx'),SheetNames = c('Stats'), col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
}