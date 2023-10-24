
library(Seurat)
library(ggplot2)

get_histogram <- function(data){
	hist <- ggplot(data, aes(x=Scores))+
  	geom_histogram(color="darkblue", fill="lightblue") + theme_classic()
	return(hist)
}


all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))


doublet_folder <- '/home/sevastopol/data/gserranos/MDS/Data/Doublet_prediction/'

plot_list <- list()
all_doublets <- data.frame(Scores=NULL, Prediction=NULL, Sample=NULL)
all_scores <- data.frame(Scores=NULL, Prediction=NULL, Sample=NULL)
for (sample in unique(all_seurat_integrated_sct$Sample)){
	doublets <- read.table(paste0(doublet_folder, '/', sample, '_doublets_results.txt'), sep='\t', header=TRUE)
	print(paste0('===== ', sample, ' ====='))
	print(table(doublets$Prediction))
	rownames(doublets) <- doublets$X
	doublets <- doublets[,-1]
	doublets$Sample <- sample
	all_scores <- rbind(all_scores, doublets)
	plot_list[[sample]] <- get_histogram(doublets)
	doublets <- doublets[doublets$Prediction == 'True', ]
	all_doublets <- rbind(all_doublets, doublets)
}

pdf('/home/sevastopol/data/gserranos/MDS/Plots/Integrated/Doublet_plots.pdf')
cowplot::plot_grid(plotlist=plot_list, ncol=2, nrow=2)
dev.off()


coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
rownames(all_scores) <- paste0( all_scores$Sample , '_',rownames(all_scores))
coords <- merge(coords, all_scores, by=0)
pdf('/home/sevastopol/data/gserranos/MDS/Plots/Integrated/Doublet_plots_UMAP.pdf')
ggplot()+
	geom_point(data= coords[coords$Prediction == 'False',], mapping= aes(x=UMAP_1, y=UMAP_2, color=Prediction) , size=0.3, alpha = 0.8)+
	geom_point(data= coords[coords$Prediction == 'True',],  mapping= aes(x=UMAP_1, y=UMAP_2, color=Prediction) , size=0.5, alpha = 1)+
	viridis::scale_color_viridis(discrete=TRUE, direction = -1)+
	theme_classic()
ggplot(coords, aes(x=UMAP_1, y=UMAP_2, color=Scores))+
	geom_point(size=0.5, alpha = 0.8)+
	viridis::scale_color_viridis(direction=1, option = 'inferno')+
	theme_classic()
dev.off()
