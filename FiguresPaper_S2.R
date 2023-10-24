library(Seurat)
library(ggplot2)
library(ComplexHeatmap)


get_gvenn <- function(data, title, scale_Y=TRUE){
	p <-  ggvenn::ggvenn(data,
	fill_color = destiny::cube_helix(length(data)),
	stroke_size = 0.4,
	show_percentage = TRUE,
	fill_alpha = 0.4,
	stroke_color = 'white',
	stroke_alpha = 1,
	stroke_linetype = 'solid',
	text_color = 'black',
	set_name_size = 4,
	text_size = 2,
	label_sep = ','
	)+ theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 7),
		plot.margin = unit(c(0, 0, 0, 0), "cm"),
		axis.text.y=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.y=element_blank()) + ggtitle(title)
	if(scale_Y){
	p <- p + scale_y_continuous(limits = c(-1, 1.5))
	}
	return(p)
}


get_umap_continous <- function(data, color_val){
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.8, size = 0.6) + theme_classic() + 
		viridis::scale_color_viridis( name=color_val) + 
		theme(legend.position='bottom', text = element_text(family = "Helvetica", size = 7), 
				line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
				plot.title = element_text(hjust = 0.5)) 
		return(p)
}

get_umap <- function(data, color_val, mapped_colors, legend_row=3){
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.7, size = 0.3) + theme_classic() + 
		scale_color_manual(values=mapped_colors) + 
		theme(legend.position='bottom', 
		text = element_text(family = "Helvetica", size = 7), 
		line = element_blank(),
		title = element_blank(), 
		axis.text.x =element_blank(), 
		axis.text.y=element_blank(),  
		legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm')) + 
		guides(color = guide_legend(nrow=legend_row, byrow=TRUE, override.aes = list(size=3, alpha=0.9)))
	return(p)
}

get_dotplot <- function(data, features){
	p <- DotPlot(data, features = features, col.min=0) + coord_flip() +
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	# viridis::scale_colour_viridis(option="plasma") +
	# scale_colour_gradient(low = 'blue', high='red') +
    scale_size_continuous(range = c(0.1,4))+
	scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted('red'), mid ='white') +
	theme(legend.position='right', axis.title = element_text(family = "Helvetica", size = 7), 
	 	axis.text = element_text(family = "Helvetica", size = 5), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 7, family = "Helvetica"),
		legend.text = element_text(size = 5, family = "Helvetica")) +
	guides(color = guide_colourbar(title = 'Avg\nExpression', barwidth = 0.5, barheight = 4), 
			size= guide_legend(title= 'Percent\nExpressed'))
	return(p)
}

get_copykat_HM <- function(Sample, neighbours, anno_title=TRUE){
	print(Sample)
	if(grepl('^GSM', Sample)){

		MDSclust_results <-  readRDS(paste0('./Data/CopyKat/Control_',Sample,'_copykat_CNA',neighbours,'_results_chr5_MDS.rds'))
	}else{
		
		MDSclust_results <-  readRDS(paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',neighbours,'_results_chr5_MDS.rds'))
	}
	ELDERclust_results <- readRDS(paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',neighbours,'_results_chr5_MDS.rds'))
	clust_results_mds   <- MDSclust_results[['cluster_results']]
	clust_results_elder <- ELDERclust_results[['cluster_results']]

	cluster_results <- cbind(clust_results_mds, clust_results_elder[, !colnames(clust_results_elder) %in% c('Band')])
	Annotation_rows <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(cluster_results$Band, '^p'), 'p', 'q'))), 'chr_arm')
	Annotation_rows$locus <- as.factor(stringr::str_extract(cluster_results$Band, '(?<=[pq]{1})[\\d]+'))

	rownames(cluster_results) <- paste0(rownames(cluster_results),  '_', cluster_results$Band)
	cluster_results <- cluster_results[, -1]

	rownames(Annotation_rows) <- paste(rownames(cluster_results))
	if(grepl('^GSM', Sample)){
		Annotation_columns <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(colnames(cluster_results), '_Elder_'), 'Control', 'Elder'))), c('Sample'))
	}else{
		Annotation_columns <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(colnames(cluster_results), '_Elder_'), 'Control', 'MDS'))), c('Sample'))
	}
	rownames(Annotation_columns) <- colnames(cluster_results)


	annotation_cols <- list()
	annotation_cols[['chr_arm']] <- viridis::viridis(2)
	names(annotation_cols[['chr_arm']]) <- levels(Annotation_rows$chr_arm)
	annotation_cols[['locus']]   <- viridis::viridis(nlevels(Annotation_rows$locus), option='turbo')
	names(annotation_cols[['locus']]) <- levels(Annotation_rows$locus)
	annotation_cols[['Sample']] <- Control_Elder_MDS_colors

	sample_distribution <- data.frame(Cell = NULL, Cluster= NULL, Sample = NULL)
	tmp <-  setNames(as.data.frame(MDSclust_results[['clustered_cells']]), 'Cluster')
	tmp$Cell <- rownames(tmp)
	tmp$Sample <- ifelse(grepl('^GSM', Sample), 'Elder', 'MDS')
	sample_distribution <- rbind(sample_distribution, tmp)
	tmp <-  setNames(as.data.frame(ELDERclust_results[['clustered_cells']]), 'Cluster')
	tmp$Cell <- rownames(tmp)
	tmp$Sample <- 'Control'
	sample_distribution <- rbind(sample_distribution, tmp)

	sample_distribution$Cluster_name <- paste0('Cluster_', sample_distribution$Sample, '_', sample_distribution$Cluster)
	# HM <- pheatmap::pheatmap(cluster_results, cluster_rows=FALSE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_row= 5, cellheight=0.5, annotation_row=Annotation_rows, annotation_col =  Annotation_columns, annotation_colors = annotation_cols, main='Heatmap from p15.33 - q35.3', silent=TRUE)
	# HM$tree_col$labels[HM$tree_col$order]

	plotter <- setNames(as.data.frame(table(sample_distribution$Cluster_name)), c('Cluster','Number_of_cells'))
	plotter$Phenotype <- stringr::str_extract(plotter$Cluster, '(?<=_)[\\w]+(?=_)')
	# Pearson_Rank  <- readRDS('./Data/Pearson_Rank_signature.rds')
	# # Pearson_Rank  <- readRDS('./Data/Pearson_Rank_Relu.rds')
	# # Zscore_LogFC  <- readRDS('./Data/Zscore_LogFC_signature.rds')
	# # Signature     <- readRDS('./Data/Signature_signature.rds')
	# # SignatureRelu <- readRDS('./Data/SignatureRelu_signature.rds')
	# # Signature_cor <- readRDS('./Data/Signature_cor_signature.rds')
	# Pearson_Rank$Cell_id <- gsub('-', '.', rownames(Pearson_Rank))
	# Pearson_Rank$Sample <- stringr::str_extract(Pearson_Rank$Cell_id, '^[A-Z0-9]+')

	# Pearson_Rank_mds <- Pearson_Rank[Pearson_Rank$Sample == Sample,]
	# Pearson_Rank_mds$Cluster <- apply(Pearson_Rank_mds, 1, FUN=function(x){MDSclust_results[['clustered_cells']][x['Cell_id']]})
	# Pearson_Rank_mds$Cluster <- paste0('Cluster_MDS_', Pearson_Rank_mds$Cluster)

	# Pearson_Rank_elder <- Pearson_Rank[Pearson_Rank$Cell_id %in% names(ELDERclust_results[['clustered_cells']]),]
	# Pearson_Rank_elder$Cluster <- apply(Pearson_Rank_elder, 1, FUN=function(x){ELDERclust_results[['clustered_cells']][x['Cell_id']]})
	# Pearson_Rank_elder$Cluster <- paste0('Cluster_Elder_', Pearson_Rank_elder$Cluster)

	# Pearson_Rank_samples <- rbind(Pearson_Rank_mds, Pearson_Rank_elder)

	# Pearson_Rank_samples2 <- reshape2::dcast(Pearson_Rank_samples[,c('Pearson_Rank', 'Cluster', 'Cell_id') ], Cell_id~Cluster, value.var = 'Pearson_Rank')
	# Pearson_Rank_samples2 <- Pearson_Rank_samples2[, !colnames(Pearson_Rank_samples2) %in% c('Cluster_MDS_NA', 'Cell_id')]

	# Pearson_Rank_samples2[is.na(Pearson_Rank_samples2)] <- 0
	row_ann1 <- ComplexHeatmap::rowAnnotation(chr_arm = Annotation_rows$chr_arm, 
											locus=Annotation_rows$locus, 
											col = annotation_cols, 
											show_legend = c(TRUE, FALSE))
	# cols <- stringr::str_extract(colnames(Pearson_Rank_samples2), '(?<=Cluster_)[A-Za-z]+')
	# if(grepl('^GSM', Sample)){
	# 	cols<- ifelse(cols == 'Elder', Control_Elder_MDS_colors['Control'] , Control_Elder_MDS_colors['Elder'])
	# }else{
	# 	cols<- ifelse(cols == 'Elder', Control_Elder_MDS_colors['Elder'] , Control_Elder_MDS_colors['MDS'])
	# }
	cols <- Control_Elder_MDS_colors[paste(Annotation_columns$Sample)]
	top_ann2 <- ComplexHeatmap::HeatmapAnnotation(Sample = Annotation_columns$Sample, 
													N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
																which='column', height = unit(1, "cm"), gp = grid::gpar(fill = cols), ylim =c(0,1500)),
													# Signature = ComplexHeatmap::anno_boxplot(as.matrix(Pearson_Rank_samples2), 
													# 			which='column', height = unit(1.5, "cm"), gp = grid::gpar(fill = cols), outline = FALSE),
													col = annotation_cols,
													show_annotation_name= anno_title,
													annotation_legend_param = list(
														Sample = list(direction = "horizontal", nrow = 1)),
													annotation_name_gp= grid::gpar(fontsize = 7, family = "Helvetica")
													)

	ht <- ComplexHeatmap::draw(ComplexHeatmap::Heatmap(as.matrix(cluster_results),column_km= 4))
	# r.dend <- ComplexHeatmap::column_dend(ht)
	column_clusters <- ComplexHeatmap::column_order(ht)
	column_names <- colnames(cluster_results)
	p <- grid::grid.grabExpr(ComplexHeatmap::draw( 
					ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
					cluster_rows = FALSE, 
					column_km=4, 
					col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
					show_column_names = FALSE,
					show_row_names = FALSE,
					top_annotation = top_ann2,
					heatmap_legend_param = list(direction = "horizontal")),
					heatmap_legend_side = "bottom", 
				annotation_legend_side = "bottom"))
	results <- list()
	results$HM <- p
	results$raw <- as.matrix(cluster_results)
	results$top_anno <- top_ann2
	results$column_clusters <- column_clusters
	results$column_names <- colnames(cluster_results)
	results$sample_distribution <- sample_distribution
	return(results)
}

get_cluster_deletion <- function(coordinates, data){
	plotter <- merge(coordinates, data, by.x='Row.names', by.y='Var1')
	p <- ggplot(plotter, aes(x = Cluster, fill = value2)) + geom_bar(position='fill') + theme_classic() +
		scale_fill_manual(values=c('neutral'="#999999", 'amplification'="#E69F00", 'deletion'="#56B4E9"))+ 
		theme(axis.text.x = element_text(angle = 45, hjust=1, size=6)) + 
		ggtitle('CNV per cluster, resolution 0.4')
	return(p)
}



get_deletions_chromosomes <- function(Data){
	control <- 'Elder'
	sample <-  setdiff(unique(Data[,'phenotype']), control)
	# names5q <- ifelse(levels(data$Var2) %in% c('5q', '5p'), 'red', 'black')
	Data$phenotype <- factor(Data$phenotype, levels = c(sample, control))
	Data$Alpha <- ifelse(Data$Var2	%in% c('5q', '5p'), 1, 0.9)
	p <- ggplot(Data, aes(x=Var2, fill = value2, alpha=Alpha, color=as.factor(Alpha))) + geom_bar(position='fill') + theme_classic()  +
		scale_fill_manual(values=c('neutral'="#999999", 'amplification'="#E69F00", 'deletion'="#56B4E9"))+ 
		theme(legend.position='right', axis.title = element_text(family = "Helvetica", size = 7), 
	 	axis.text = element_text(family = "Helvetica", size = 5), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 7, family = "Helvetica"),
		legend.text = element_text(size = 5, family = "Helvetica")) +
		scale_alpha(range = c(0.8, 1), guide='none') + labs(x = "Locus", y='Percent of cells') + 
		scale_color_manual(values = c('grey', 'black'), guide='none')#+ facet_wrap(.~phenotype, nrow=2)
	return(p)
}



SAMPLE_NAME <- '5qSamples'
sc_data <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', SAMPLE_NAME, '_Annotated_final.rds'))

Idents(sc_data) <- 'integrated_snn_res.0.8_sub'

clusters <- setNames(as.data.frame(sc_data$integrated_snn_res.0.8_sub), c('Cluster'))


Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722')
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog')
Control_Elder_MDS_colors <- setNames(c('#002642','#840032', '#e59500'), c('Control', 'Elder', 'MDS'))

K_neighors <- 80
CopyKat_MDS_1 <- get_copykat_HM('SMD34459', K_neighors, anno_title=FALSE)
CopyKat_MDS_2 <- get_copykat_HM('SMD35109', K_neighors, anno_title=FALSE)
CopyKat_MDS_3 <- get_copykat_HM('SMD37209', K_neighors, anno_title=FALSE)
CopyKat_MDS_4 <- get_copykat_HM('SMD35303', K_neighors, anno_title=FALSE)

dev.off()
CopyKat_control <- get_copykat_HM('GSM5460411', K_neighors)
dev.off()
min_val <- round(min(range(CopyKat_MDS$raw), range(CopyKat_control$raw)), 2)
max_val <- round(max(range(CopyKat_MDS$raw), range(CopyKat_control$raw)), 2)
min_mds <- round(min(CopyKat_MDS$raw), 2)
max_mds <- round(max(CopyKat_MDS$raw), 2)
min_ctr <- round(min(CopyKat_control$raw), 2)
max_ctr <- round(max(CopyKat_control$raw), 2)

bands_2_show <- c('p15.33', 'p15.2', 'p15.1', 'p14', 'p13.3', 'p13.1', 'p12', 'q11.1', 'q12', 'q13.1',
					'q13.2', 'q13.3', 'q14', 'q15', 'q21', 'q22', 'q23.1', 'q23.2', 'q23.3', 'q31.1', 'q31.2', 
					'q31.3', 'q32', 'q33.1', 'q33.2', 'q33.3', 'q34', 'q35.1', 'q35.2', 'q35.2', 'q35.3')
positions <- c()
for (band in bands_2_show){
	all_bands <- stringr::str_extract(rownames(CopyKat_MDS_1$raw), '(?<=_)[\\w\\d\\.]+')
	position <- grep(band , all_bands)[1]
	positions <- c(positions, position)
}

col_fun = circlize::colorRamp2(seq(-0.3, 0.3, by=(0.3-(-0.3))/(100-1)),
 , colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100))

ha = ComplexHeatmap::rowAnnotation(Locus = ComplexHeatmap::anno_mark(at = positions, labels = bands_2_show, side = "left",
labels_gp = grid::gpar(col = "black", fontsize = 7.5), padding = 0.7))
ht_list =   ComplexHeatmap::Heatmap(CopyKat_MDS_1$raw,cluster_rows = FALSE, name = 'CNA Score',
			column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
			col =col_fun,
			show_column_names = FALSE,
			show_row_names = FALSE, 
			top_annotation = CopyKat_MDS_1$top_anno,
			left_annotation = ha,
			show_heatmap_legend = TRUE) + 
			ComplexHeatmap::Heatmap(CopyKat_MDS_2$raw,cluster_rows = FALSE, name = 'CNA Score SMD35109',
			column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
			col =col_fun,
			show_column_names = FALSE,
			show_row_names = FALSE, 
			top_annotation = CopyKat_MDS_2$top_anno,
			show_heatmap_legend = FALSE) + 
			ComplexHeatmap::Heatmap(CopyKat_MDS_3$raw,cluster_rows = FALSE, name = 'CNA Score SMD37209',
			column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
			col =col_fun,
			show_column_names = FALSE,
			show_row_names = FALSE, 
			top_annotation = CopyKat_MDS_3$top_anno,
			show_heatmap_legend = FALSE) + 
			ComplexHeatmap::Heatmap(CopyKat_MDS_4$raw,cluster_rows = FALSE, name = 'CNA Score SMD35303',
			column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
			col =col_fun,
			show_column_names = FALSE,
			show_row_names = FALSE, 
			top_annotation = CopyKat_MDS_4$top_anno,
			show_heatmap_legend = FALSE) + 
			ComplexHeatmap::Heatmap(CopyKat_control$raw,cluster_rows = FALSE, 
			column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
			col = col_fun,
			show_column_names = FALSE,
			show_row_names = FALSE,
			top_annotation = CopyKat_control$top_anno,
			show_heatmap_legend = FALSE)


all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds')
results_CK <- readRDS('./Data/CopyKat/results_CK.rds')
plotter_list <- readRDS('./Data/CopyKat/plotter_list.rds')
cells_selected_CASPER <- readRDS('./Data/cells_selected_CASPER.rds')

plot_list <- list()
for(Sample in c("SMD34459", "SMD35109", "SMD35303", "SMD37209")){
plot_list[[Sample]] <- get_gvenn(list(
						CopyKat = all_5q_selected_cells[all_5q_selected_cells$Sample == Sample, 'Cell_id'],
						CASPER = cells_selected_CASPER[cells_selected_CASPER$Sample == Sample, 'Cell_id']), 
						title = Sample)
}



pdf('./Plots/PaperFigures/FigS2.pdf')
cowplot::plot_grid(
	cowplot::plot_grid(grid::grid.grabExpr(
		ComplexHeatmap::draw(ht_list, merge_legend = TRUE, 
		heatmap_legend_side = "bottom", 
		annotation_legend_side = "bottom", ht_gap = unit(c(6, 6, 6, 6), "mm")))),
	cowplot::plot_grid(plotlist=plot_list),
ncol=1, rel_heights = c(1, 0.5))
dev.off()