


library(ggplot2)


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



get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}


annotate_band <- function(chr, pos, cytobandAnnotatiojn = cyto_ann){
	tmp <- cytobandAnnotatiojn[cytobandAnnotatiojn$V1 == paste0('chr', chr)  ,]
	band <- as.character(tmp[tmp$V2 <= pos & tmp$V3 >= pos,'V4'])
	if (length(band) ==0 ){
		# band <- chr
		band <- as.character(tmp[nrow(tmp), 'V4'])
	}
	return(band)
}

get_boxplots <- function(data, title, hline = NULL){
	ncolors <- length(unique(data$phenotype))
	p <- ggplot(data, aes(x =band, y = value, fill = phenotype)) + geom_boxplot(alpha = 0.8, outlier.alpha = 0.2) + 
	theme_classic() + scale_fill_manual(values=destiny::cube_helix(ncolors)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position='bottom') + 
	ggtitle(title) + ylab('Estimated copy numbers') + xlab('Chr5 genomic band') + 
	if (!is.null(hline)){
		p <- p + geom_hline(yintercept=hline, linetype="dashed")
	}
	return(p)
}

get_5q_cells <- function(cell, threshold){
	values_cna   <- values_5q_2_plot[values_5q_2_plot$variable == cell & values_5q_2_plot$band %in% bands_2_filter, 'value']
	if (sum(values_cna < threshold) > 0.5 * length(values_cna)){
			return(cell)
		}
}


cluster_and_median <- function(data, k, tag){
	results <- list()
	print('Wow! that will take a while!')
	data_no_band <- data[, !colnames(data) %in% c('band')]
	clust_obj <- hclust(dist(t(data_no_band)))
	clustered_cells <- cutree(clust_obj, k)
	cluster_results <- data.frame(Band= data$band)
	for (cluster in unique(clustered_cells)){
		tmp <- data_no_band[, names(which(clustered_cells == cluster)), drop=FALSE]
		tmp <- setNames(as.data.frame(matrixStats::rowMedians(as.matrix(tmp))), c(paste0('Cluster_',tag,'_', cluster)))
		cluster_results <- cbind(cluster_results, tmp)
	}
	results[['cluster_results']] <- cluster_results
	results[['clustered_cells']] <- clustered_cells
	return(results)
}




# SMD34459
# SMD35109
# SMD35303
# SMD37209

all_samples <- c('FS-0634-post',	'FS-0406-post')
for (Sample in all_samples){
	print(Sample)
	if (!file.exists(paste0('/home/tereshkova/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_CNA_results_chr5.txt'))){
		# predictions <- read.table(paste0('/home/tereshkova/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_prediction.txt'), header=TRUE)
		values      <- read.table(paste0('/home/tereshkova/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_CNA_results.txt'), header=TRUE)
		write.table(values[values$chrom=="5",], paste0('/home/tereshkova/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_CNA_results_chr5.txt'), quote=FALSE, sep='\t', row.names=FALSE)
	}
}

all_results <- list()
for (Sample in all_samples){
	print(Sample)
	tmp <- read.table(paste0('/home/tereshkova/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_CNA_results_chr5.txt'), header=TRUE)
	print(dim(tmp))
	all_results[[Sample]] <- tmp
}


# annotate_band(23, values[1,2])

cyto_ann <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/cytoBand.txt')
for (Sample in all_samples){
	print(Sample)
	tmp <- all_results[[Sample]]
	bands <- apply(tmp[, 1:3], 1, function(x){annotate_band(x[1], x[2])})
	tmp$band <- bands
	all_results[[Sample]] <- tmp
}


bands_2_keep <- paste0('q', seq(13,33))
bands_2_filter <- paste0('q', seq(15,31))



K_neighors = 80
for(Sample in all_samples){
	print(Sample)
	if(Sample == 'GSM5460411'){
		values <- all_results[[Sample]]
		values_mds   <-  values[ , c('band', grep(paste0('^', Sample), colnames(values), value=TRUE))]
		values_elder <-  values[ , c('band', grep(paste0('^', Sample), colnames(values), value=TRUE, invert = TRUE))]

		MDSclust_results   <- cluster_and_median(values_mds, K_neighors, 'Control')
		ELDERclust_results <- cluster_and_median(values_elder, K_neighors, 'Elder')
		saveRDS(MDSclust_results,   paste0('./Data/CopyKat/Control_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
		saveRDS(ELDERclust_results, paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
	}else{
		values <- all_results[[Sample]]
		values_mds   <-  values[ , c('band', grep('^FS', colnames(values), value=TRUE))]
		values_elder <-  values[ , c('band', grep('^GSM', colnames(values), value=TRUE))]

		MDSclust_results   <- cluster_and_median(values_mds, K_neighors, 'MDS')
		ELDERclust_results <- cluster_and_median(values_elder, K_neighors, 'Elder')
		saveRDS(MDSclust_results,   paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
		saveRDS(ELDERclust_results, paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
	}
}



Control_Elder_MDS_colors <- setNames(c('#002642','#840032', '#e59500'), c('Control', 'Elder', 'MDS'))
all_clusters <- data.frame(Cell_id=NULL, Cluster=NULL, Cluster_ID = NULL)
cluster_results <- NULL
for(Sample in all_samples){
	message(Sample)
	MDSclust_results   <- readRDS(paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
	ELDERclust_results <- readRDS(paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
	clust_results_mds   <- MDSclust_results[['cluster_results']]
	colnames(clust_results_mds) <- paste0(names(clust_results_mds), '_', Sample)
	if (Sample == all_samples[1]){
		cluster_results <- clust_results_mds
	}else{
		cluster_results <- cbind(cluster_results, clust_results_mds[, !colnames(clust_results_mds) %in% c(paste0('Band', '_', Sample))])
	}
	tmp <- setNames(as.data.frame(MDSclust_results[['clustered_cells']]), 'Cluster')
	tmp$Cell_id <- rownames(tmp)
	tmp$Cluster_ID <- paste0('Cluster_MDS_', tmp$Cluster, '_', Sample)
	all_clusters <- rbind(all_clusters, tmp)


	tmp <- setNames(as.data.frame(ELDERclust_results[['clustered_cells']]), 'Cluster')
	tmp$Cell_id <- rownames(tmp)
	tmp$Cluster_ID <- paste0('Cluster_Elder_', tmp$Cluster, '_Elder')
	all_clusters <- rbind(all_clusters, tmp)

	clust_results_elder <- ELDERclust_results[['cluster_results']]
	Annotation_rows <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(clust_results_elder$Band, '^p'), 'p', 'q'))), 'chr_arm')
	rownames(cluster_results) <- paste0(rownames(cluster_results),  '_', clust_results_elder$Band)
	Annotation_rows$locus <- as.factor(stringr::str_extract(clust_results_elder$Band, '(?<=[pq]{1})[\\d]+'))
	colnames(clust_results_elder) <- paste0(colnames(clust_results_elder), '_Elder')
	cluster_results <- cbind(cluster_results, clust_results_elder[, !colnames(clust_results_elder) %in% c('Band_Elder')])

	cluster_results <- cluster_results[, -1]

	rownames(Annotation_rows) <- paste(rownames(cluster_results))
	Annotation_columns <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(colnames(cluster_results), '_Elder_'), 'Elder', 'MDS'))), c('Phenotype'))
	rownames(Annotation_columns) <- colnames(cluster_results)

	annotation_cols <- list()
	annotation_cols[['chr_arm']] <- viridis::viridis(2)
	names(annotation_cols[['chr_arm']]) <- levels(Annotation_rows$chr_arm)
	annotation_cols[['locus']]   <- viridis::viridis(nlevels(Annotation_rows$locus), option='turbo')
	names(annotation_cols[['locus']]) <- levels(Annotation_rows$locus)
	annotation_cols[['Phenotype']] <- viridis::viridis(nlevels(Annotation_columns$Phenotype), option='turbo')
	names(annotation_cols[['Phenotype']]) <- levels(Annotation_columns$Phenotype)
	Annotation_columns$Sample <- stringr::str_extract(rownames(Annotation_columns), 'FS-\\d+-\\w+')
	Annotation_columns$Sample[is.na(Annotation_columns$Sample)] <- 'Elder'
	annotation_cols[['Sample']] <- c('#495867', '#577399')
	names(annotation_cols[['Sample']]) <-  all_samples

	plotter <- setNames(as.data.frame(table(all_clusters$Cluster_ID)), c('Cluster','Number_of_cells'))
	# plotter <- plotter[stringr::str_detect(plotter$Cluster , Sample) | stringr::str_detect(plotter$Cluster , 'Elder'),]
	plotter$Phenotype <- ifelse(stringr::str_detect(plotter$Cluster, '_FS_'), annotation_cols[['Phenotype']][['MDS']], annotation_cols[['Phenotype']][['Elder']])
	rownames(plotter) <- paste(plotter$Cluster)


	top_ann2 <-ComplexHeatmap::HeatmapAnnotation(  Phenotype = Annotation_columns$Phenotype, 
													Sample= Annotation_columns$Sample,
													N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
														which='column', height = unit(2, "cm"), gp = grid::gpar(fill = plotter$Phenotype))
													)
	pdf(paste0('./Plots/CopyKat/Test_POST',Sample,'.pdf'), width=8.25, height=11.75)
	print(ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
					name ='CNVA', 
					cluster_rows = FALSE, 
					column_km=8, 
					col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
					show_column_names = FALSE,
					show_row_names = FALSE,
					show_column_dend = FALSE))
	dev.off()
}



all_5q_selected_cells <- data.frame(Cell_id=NULL, Sample=NULL)
for (Sample in all_samples){
	if(Sample %in% c('FS-0406-post')){cluster_2_keep <-  c(1,2)}else{cluster_2_keep <-  c(1)}
	selected_cells <- c()
	for (clust in cluster_2_keep){
		selected_cells <- c(selected_cells, results_CK[[Sample]]$column_clusters[[paste(clust)]])
	}
	selected_cells <- results_CK[[Sample]]$column_names[selected_cells]
	selected_cells <- selected_cells[stringr::str_detect(selected_cells, 'MDS')]
	all_cells <- results_CK[[Sample]]$sample_distribution
	selected_cells <- all_cells[all_cells$Cluster_name %in% selected_cells, 'Cell']
	selected_cells <- sub('\\.', '-',selected_cells)
	all_5q_selected_cells <- rbind(all_5q_selected_cells, data.frame(Cell_id=selected_cells, Sample=Sample))
}
saveRDS(all_5q_selected_cells,  './Data/CopyKat/all_5q_CopyKat_depleted_cells_POST.rds')

