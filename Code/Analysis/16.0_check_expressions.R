library(Seurat)
library(ggplot2)

library(ggpubr)
library(rstatix)


SAMPLE_NAME <- '5qSamples'
sc_data <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', SAMPLE_NAME, '_Annotated_final.rds'))

DefaultAssay(sc_data) <- 'RNA'
genes_2_check <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/genes_in_the_deleted_region/13-33-Table_1.tsv', sep='\t', header=T)
genes_2_check <- genes_2_check$Gene.name
genes_2_check <- genes_2_check[genes_2_check!= '']

genes_2_check <- genes_2_check[genes_2_check %in% rownames(sc_data)]

gene <- genes_2_check[1]


plotter <- FetchData(sc_data, vars = c(genes_2_check, 'Sample', 'cell_5q'))

get_boxplot <- function(gene, tmp = plotter){
	p <- ggplot(tmp, aes(x=Sample, y=.data[[gene]],  fill=cell_5q)) + geom_boxplot(outlier.alpha=0.5) + 
	theme_classic() + scale_fill_manual(values = c('#e63946', '#457b9d')) + 
	# ggtitle(gene) + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1),
	legend.position = 'top')
	return(p)
}




pdf('./Plots/expression5qCheck.pdf')
plotter_list <- list()
counter <- 1
for (gene in genes_2_check){
	print(gene)
	plotter_list[[gene]] <- get_boxplot(gene)
	if(counter %% 4 == 0){
		print(cowplot::plot_grid(plotlist=plotter_list, ncol=2))
		plotter_list <- list()
		counter <- 1
		# break
	}else{
		counter <- counter + 1
	}
}
# ggplot(plotter, aes(x=cell_5q, y=.data[[gene]],  fill=Sample)) + geom_boxplot() + 
# theme_classic()
dev.off()