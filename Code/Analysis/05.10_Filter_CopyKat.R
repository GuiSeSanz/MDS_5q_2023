


library(ggplot2)



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

all_samples <- c('SMD34459', 'SMD35109', 'SMD35303', 'SMD37209', 'GSM5460411')
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



for( K_neighors in c(10, 20, 40, 80)){

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
			values_mds   <-  values[ , c('band', grep('^SMD', colnames(values), value=TRUE))]
			values_elder <-  values[ , c('band', grep('^GSM', colnames(values), value=TRUE))]

			MDSclust_results   <- cluster_and_median(values_mds, K_neighors, 'MDS')
			ELDERclust_results <- cluster_and_median(values_elder, K_neighors, 'Elder')
			saveRDS(MDSclust_results,   paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
			saveRDS(ELDERclust_results, paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
		}
	}
}
quit(save='no')	

	# clust_results_mds   <- MDSclust_results[['cluster_results']]
	# clust_results_elder <- ELDERclust_results[['cluster_results']]

	# cluster_results <- cbind(clust_results_mds, clust_results_elder[, !colnames(clust_results_elder) %in% c('Band')])

	# rownames(cluster_results) <- paste0(rownames(cluster_results),  '_', cluster_results$Band)
	# cluster_results <- cluster_results[, -1]

	# Annotation_rows <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(values$band, '^p'), 'p', 'q'))), 'chr_arm')
	# Annotation_rows$locus <- as.factor(stringr::str_extract(values$band, '(?<=[pq]{1})[\\d]+'))
	# rownames(Annotation_rows) <- paste(rownames(cluster_results))
	# Annotation_columns <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(colnames(cluster_results), '_Elder_'), 'Elder', 'MDS'))), c('Sample'))
	# rownames(Annotation_columns) <- colnames(cluster_results)


	# annotation_cols <- list()
	# annotation_cols[['chr_arm']] <- viridis::viridis(2)
	# names(annotation_cols[['chr_arm']]) <- levels(Annotation_rows$chr_arm)
	# annotation_cols[['locus']]   <- viridis::viridis(nlevels(Annotation_rows$locus), option='turbo')
	# names(annotation_cols[['locus']]) <- levels(Annotation_rows$locus)
	# annotation_cols[['Sample']] <- viridis::viridis(nlevels(Annotation_columns$Sample), option='turbo')
	# names(annotation_cols[['Sample']]) <- levels(Annotation_columns$Sample)

	# sample_distribution <- data.frame(Cell = NULL, Cluster= NULL, Sample = NULL)
	# tmp <-  setNames(as.data.frame(MDSclust_results[['clustered_cells']]), 'Cluster')
	# tmp$Cell <- rownames(tmp)
	# tmp$Sample <- 'MDS'
	# sample_distribution <- rbind(sample_distribution, tmp)
	# tmp <-  setNames(as.data.frame(ELDERclust_results[['clustered_cells']]), 'Cluster')
	# tmp$Cell <- rownames(tmp)
	# tmp$Sample <- 'Elder'
	# sample_distribution <- rbind(sample_distribution, tmp)

	# sample_distribution$Cluster_name <- paste0('Cluster_', sample_distribution$Sample, '_', sample_distribution$Cluster)
	# # HM <- pheatmap::pheatmap(cluster_results, cluster_rows=FALSE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_row= 5, cellheight=0.5, annotation_row=Annotation_rows, annotation_col =  Annotation_columns, annotation_colors = annotation_cols, main='Heatmap from p15.33 - q35.3', silent=TRUE)
	# # HM$tree_col$labels[HM$tree_col$order]

	# plotter <- setNames(as.data.frame(table(sample_distribution$Cluster_name)), c('Cluster','Number_of_cells'))
	# plotter$Phenotype <- stringr::str_extract(plotter$Cluster, '(?<=_)[\\w]+(?=_)')
	# # plotter$Cluster <- as.factor(plotter$Cluster)
	# # # plotter <- plotter[match(HM$tree_col$labels[HM$tree_col$order], plotter$Cluster),]
	# # plotter$Cluster <- factor(plotter$Cluster, ordered = TRUE, levels=c(HM$tree_col$labels[HM$tree_col$order]))


# 	top_ann1 <- ComplexHeatmap::HeatmapAnnotation(Sample = Annotation_columns$Sample, 
# 												N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
# 															which='column', height = unit(2, "cm")),
# 												show_legend = c(FALSE, FALSE), annotation_name_side = "left",
# 												col = annotation_cols)
# 	top_ann2 <- ComplexHeatmap::HeatmapAnnotation(Sample = Annotation_columns$Sample, 
# 												N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
# 															which='column', height = unit(2, "cm")),
# 												col = annotation_cols)

# 	row_ann1 <- ComplexHeatmap::rowAnnotation(chr_arm = Annotation_rows$chr_arm, 
# 												locus=Annotation_rows$locus, 
# 												col = annotation_cols, 
# 												show_legend = c(FALSE, FALSE))

# 	row_ann2 <- ComplexHeatmap::rowAnnotation(chr_arm = Annotation_rows$chr_arm, 
# 												locus=Annotation_rows$locus, 
# 												col = annotation_cols)



# 	pdf(paste0('./Plots/CopyKat/',Sample,'_HM.pdf'))
# 	print(ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
# 				cluster_rows = FALSE, 
# 				column_km=2, 
# 				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				show_column_names = FALSE,
# 				show_row_names = FALSE,
# 				top_annotation = top_ann1,
# 				show_heatmap_legend = FALSE) + 
# 	row_ann2 + 
# 	ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
# 				cluster_rows = FALSE, 
# 				column_km=3, 
# 				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				show_column_names = FALSE,
# 				show_row_names = FALSE,
# 				top_annotation = top_ann2))
	
# 	print(ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
# 				cluster_rows = FALSE, 
# 				column_km=4, 
# 				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				show_column_names = FALSE,
# 				show_row_names = FALSE,
# 				top_annotation = top_ann1,
# 				show_heatmap_legend = FALSE) + 
# 	row_ann2 + 
# 	ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
# 				cluster_rows = FALSE, 
# 				column_km=5, 
# 				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				show_column_names = FALSE,
# 				show_row_names = FALSE,
# 				top_annotation = top_ann2))
# 	dev.off()

# }
# }

# Pearson_Rank  <- readRDS('./Data/Pearson_Rank_signature.rds')
# Pearson_Rank  <- readRDS('./Data/Pearson_Rank_Relu.rds')
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

# # Pearson_Rank_samples2[is.na(Pearson_Rank_samples2)] <- 0

# cols <- ifelse(stringr::str_detect(colnames(Pearson_Rank_samples2), '_MDS_'), '#30123BFF', '#7A0403FF')

# pdf('./Plots/CopyKat/Test.pdf')
# top_ann2 <- ComplexHeatmap::HeatmapAnnotation(Sample = Annotation_columns$Sample, 
# 												N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
# 												which='column', height = unit(2, "cm")),
# 												Signature = ComplexHeatmap::anno_boxplot(as.matrix(Pearson_Rank_samples2), 
# 												which='column', height = unit(2, "cm"), gp = grid::gpar(fill = cols)),
# 												col = annotation_cols)
# ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
# 				cluster_rows = FALSE, 
# 				column_km=4, 
# 				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
# 				show_column_names = FALSE,
# 				show_row_names = FALSE,
# 				top_annotation = top_ann2)
# dev.off()









all_samples <- c('SMD34459', 'SMD35109', 'SMD35303', 'SMD37209')
cyto_ann <- read.table('/home/sevastopol/data/gserranos/MDS/Data/Annotation/cytoBand.txt')
all_results <- list()
for (Sample in all_samples){
	print(Sample)
	tmp <- read.table(paste0('/home/sevastopol/data/gserranos/MDS/Data/CopyKat/', Sample , '/',Sample,'_copykat_CNA_results_chr5.txt'), header=TRUE)
	print(dim(tmp))
	bands <- apply(tmp[, 1:3], 1, function(x){annotate_band(x[1], x[2])})
	tmp$band <- bands
	all_results[[Sample]] <- tmp
}


bands_2_keep <- paste0('q', seq(13,33))
bands_2_filter <- paste0('q', seq(15,31))



for(Sample in all_samples){
	values <- all_results[[Sample]]
	for( K_neighors in c(10, 20, 40, 80)){
		print(Sample)
		MDSclust_results   <- readRDS(paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
		ELDERclust_results <- readRDS(paste0('./Data/CopyKat/elder_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
		clust_results_mds   <- MDSclust_results[['cluster_results']]
		clust_results_elder <- ELDERclust_results[['cluster_results']]

		cluster_results <- cbind(clust_results_mds, clust_results_elder[, !colnames(clust_results_elder) %in% c('Band')])

		rownames(cluster_results) <- paste0(rownames(cluster_results),  '_', cluster_results$Band)
		cluster_results <- cluster_results[, -1]

		Annotation_rows <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(values$band, '^p'), 'p', 'q'))), 'chr_arm')
		Annotation_rows$locus <- as.factor(stringr::str_extract(values$band, '(?<=[pq]{1})[\\d]+'))
		rownames(Annotation_rows) <- paste(rownames(cluster_results))
		Annotation_columns <- setNames(as.data.frame(as.factor(ifelse(stringr::str_detect(colnames(cluster_results), '_Elder_'), 'Elder', 'MDS'))), c('Sample'))
		rownames(Annotation_columns) <- colnames(cluster_results)


		annotation_cols <- list()
		annotation_cols[['chr_arm']] <- viridis::viridis(2)
		names(annotation_cols[['chr_arm']]) <- levels(Annotation_rows$chr_arm)
		annotation_cols[['locus']]   <- viridis::viridis(nlevels(Annotation_rows$locus), option='turbo')
		names(annotation_cols[['locus']]) <- levels(Annotation_rows$locus)
		annotation_cols[['Sample']] <- viridis::viridis(nlevels(Annotation_columns$Sample), option='turbo')
		names(annotation_cols[['Sample']]) <- levels(Annotation_columns$Sample)

		sample_distribution <- data.frame(Cell = NULL, Cluster= NULL, Sample = NULL)
		tmp <-  setNames(as.data.frame(MDSclust_results[['clustered_cells']]), 'Cluster')
		tmp$Cell <- rownames(tmp)
		tmp$Sample <- 'MDS'
		sample_distribution <- rbind(sample_distribution, tmp)
		tmp <-  setNames(as.data.frame(ELDERclust_results[['clustered_cells']]), 'Cluster')
		tmp$Cell <- rownames(tmp)
		tmp$Sample <- 'Elder'
		sample_distribution <- rbind(sample_distribution, tmp)

		sample_distribution$Cluster_name <- paste0('Cluster_', sample_distribution$Sample, '_', sample_distribution$Cluster)
		# HM <- pheatmap::pheatmap(cluster_results, cluster_rows=FALSE, cluster_cols=TRUE, show_rownames=FALSE, show_colnames=TRUE, fontsize_row= 5, cellheight=0.5, annotation_row=Annotation_rows, annotation_col =  Annotation_columns, annotation_colors = annotation_cols, main='Heatmap from p15.33 - q35.3', silent=TRUE)
		# HM$tree_col$labels[HM$tree_col$order]

		plotter <- setNames(as.data.frame(table(sample_distribution$Cluster_name)), c('Cluster','Number_of_cells'))
		plotter$Phenotype <- stringr::str_extract(plotter$Cluster, '(?<=_)[\\w]+(?=_)')
		

		Pearson_Rank  <- readRDS('./Data/Pearson_Rank_signature.rds')
		# Pearson_Rank  <- readRDS('./Data/Pearson_Rank_Relu.rds')
		# Zscore_LogFC  <- readRDS('./Data/Zscore_LogFC_signature.rds')
		# Signature     <- readRDS('./Data/Signature_signature.rds')
		# SignatureRelu <- readRDS('./Data/SignatureRelu_signature.rds')
		# Signature_cor <- readRDS('./Data/Signature_cor_signature.rds')
		Pearson_Rank$Cell_id <- gsub('-', '.', rownames(Pearson_Rank))
		Pearson_Rank$Sample <- stringr::str_extract(Pearson_Rank$Cell_id, '^[A-Z0-9]+')
		print(range(Pearson_Rank[Pearson_Rank$Sample == 'GSM5460411', 'Pearson_Rank']))

		Pearson_Rank_mds <- Pearson_Rank[Pearson_Rank$Sample == Sample,]
		Pearson_Rank_mds$Cluster <- apply(Pearson_Rank_mds, 1, FUN=function(x){MDSclust_results[['clustered_cells']][x['Cell_id']]})
		Pearson_Rank_mds$Cluster <- paste0('Cluster_MDS_', Pearson_Rank_mds$Cluster)

		Pearson_Rank_elder <- Pearson_Rank[Pearson_Rank$Cell_id %in% names(ELDERclust_results[['clustered_cells']]),]
		Pearson_Rank_elder$Cluster <- apply(Pearson_Rank_elder, 1, FUN=function(x){ELDERclust_results[['clustered_cells']][x['Cell_id']]})
		Pearson_Rank_elder$Cluster <- paste0('Cluster_Elder_', Pearson_Rank_elder$Cluster)

		Pearson_Rank_samples <- rbind(Pearson_Rank_mds, Pearson_Rank_elder)

		Pearson_Rank_samples2 <- reshape2::dcast(Pearson_Rank_samples[,c('Pearson_Rank', 'Cluster', 'Cell_id') ], Cell_id~Cluster, value.var = 'Pearson_Rank')
		Pearson_Rank_samples2 <- Pearson_Rank_samples2[, !colnames(Pearson_Rank_samples2) %in% c('Cluster_MDS_NA', 'Cell_id')]

		# Pearson_Rank_samples2[is.na(Pearson_Rank_samples2)] <- 0
		row_ann1 <- ComplexHeatmap::rowAnnotation(chr_arm = Annotation_rows$chr_arm, 
												locus=Annotation_rows$locus, 
												col = annotation_cols, 
												show_legend = c(TRUE, FALSE))
		cols <- ifelse(stringr::str_detect(colnames(Pearson_Rank_samples2), '_MDS_'), '#30123BFF', '#7A0403FF')

		pdf(paste0('./Plots/CopyKat/HM_Clustering_',Sample,'_',K_neighors,'.pdf'))
		top_ann2 <- ComplexHeatmap::HeatmapAnnotation(Sample = Annotation_columns$Sample, 
														N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
														which='column', height = unit(2, "cm"), gp = grid::gpar(fill = cols)),
														Signature = ComplexHeatmap::anno_boxplot(as.matrix(Pearson_Rank_samples2), 
														which='column', height = unit(1.5, "cm"), gp = grid::gpar(fill = cols), outline = FALSE),
														col = annotation_cols)
		print(ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
						cluster_rows = FALSE, 
						column_km=4, 
						col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
						show_column_names = FALSE,
						show_row_names = FALSE,
						top_annotation = top_ann2))
		dev.off()
	}
}


# try all samples at once


K_neighors = 80

# The elder clustering is always the same
ELDERclust_results <- readRDS(paste0('./Data/CopyKat/elder_',all_samples[1],'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))

all_clusters <- data.frame(Cell_id=NULL, Cluster=NULL, Cluster_ID = NULL)
cluster_results <- NULL
for(Sample in all_samples){
	# values <- all_results[[Sample]]
	print(Sample)
	MDSclust_results   <- readRDS(paste0('./Data/CopyKat/MDS5q_',Sample,'_copykat_CNA',K_neighors,'_results_chr5_MDS.rds'))
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
}

# add cluster info for elder
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
Annotation_columns$Sample <- stringr::str_extract(rownames(Annotation_columns), '(?<=_)([A-Za-z0-9]+$)')
annotation_cols[['Sample']] <- c('#495867', '#577399', '#bdd5ea','#f7f7ff', '#272729')
names(annotation_cols[['Sample']]) <-  c("SMD34459", "SMD35109", "SMD35303", "SMD37209", 'Elder')


# sample_distribution <- data.frame(Cell = NULL, Cluster= NULL, Sample = NULL)
# tmp <-  setNames(as.data.frame(MDSclust_results[['clustered_cells']]), 'Cluster')
# tmp$Cell <- rownames(tmp)
# tmp$Sample <- 'MDS'
# sample_distribution <- rbind(sample_distribution, tmp)
# tmp <-  setNames(as.data.frame(ELDERclust_results[['clustered_cells']]), 'Cluster')
# tmp$Cell <- rownames(tmp)
# tmp$Sample <- 'Elder'
# sample_distribution <- rbind(sample_distribution, tmp)

# sample_distribution$Cluster_name <- paste0('Cluster_', sample_distribution$Sample, '_', sample_distribution$Cluster)


plotter <- setNames(as.data.frame(table(all_clusters$Cluster_ID)), c('Cluster','Number_of_cells'))
plotter$Phenotype <- ifelse(stringr::str_detect(plotter$Cluster, '_MDS_'), annotation_cols[['Phenotype']][['MDS']], annotation_cols[['Phenotype']][['Elder']])
rownames(plotter) <- paste(plotter$Cluster)




Pearson_Rank  <- readRDS('./Data/Pearson_Rank_signature.rds')
# Pearson_Rank  <- readRDS('./Data/Pearson_Rank_Relu.rds')
# Zscore_LogFC  <- readRDS('./Data/Zscore_LogFC_signature.rds')
# Signature     <- readRDS('./Data/Signature_signature.rds')
# SignatureRelu <- readRDS('./Data/SignatureRelu_signature.rds')
# Signature_cor <- readRDS('./Data/Signature_cor_signature.rds')
Pearson_Rank$Cell_id <- gsub('-', '.', rownames(Pearson_Rank))
Pearson_Rank$Sample <- stringr::str_extract(Pearson_Rank$Cell_id, '^[A-Z0-9]+')
# remove the 'normal' MDS samples
elder_samples <- c('GSM5460411', 'GSM5460412', 'GSM5460413')
Pearson_Rank <- Pearson_Rank[Pearson_Rank$Sample %in% c(all_samples, elder_samples),]

# Pearson_Rank_mds <- Pearson_Rank[Pearson_Rank$Sample == Sample,]
Pearson_Rank <- merge(Pearson_Rank, all_clusters[, c('Cell_id', 'Cluster_ID')], by='Cell_id')

# Pearson_Rank_mds$Cluster <- paste0('Cluster_MDS_', Pearson_Rank_mds$Cluster)

# Pearson_Rank_elder <- Pearson_Rank[Pearson_Rank$Cell_id %in% names(ELDERclust_results[['clustered_cells']]),]
# Pearson_Rank_elder$Cluster <- apply(Pearson_Rank_elder, 1, FUN=function(x){ELDERclust_results[['clustered_cells']][x['Cell_id']]})
# Pearson_Rank_elder$Cluster <- paste0('Cluster_Elder_', Pearson_Rank_elder$Cluster)

# Pearson_Rank_samples <- rbind(Pearson_Rank_mds, Pearson_Rank_elder)

Pearson_Rank_samples2 <- reshape2::dcast(Pearson_Rank[,c('Pearson_Rank', 'Cluster_ID', 'Cell_id') ], Cell_id~Cluster_ID, value.var = 'Pearson_Rank')
Pearson_Rank_samples2 <- Pearson_Rank_samples2[, !colnames(Pearson_Rank_samples2) %in% c('Cluster_MDS_NA', 'Cell_id')]

# row_ann1 <- ComplexHeatmap::rowAnnotation(chr_arm = Annotation_rows$chr_arm, 
# 										locus=Annotation_rows$locus, 
# 										col = annotation_cols, 
# 										show_legend = c(TRUE, FALSE))



plotter <- plotter[match(colnames(cluster_results), plotter$Cluster),]

Pearson_Rank_samples2 <- Pearson_Rank_samples2[, colnames(cluster_results)]
cols <- ifelse(stringr::str_detect(colnames(Pearson_Rank_samples2), 'MDS'), '#7A0403FF', '#30123BFF')

Pearson_Rank_samples2_lines <- setNames(as.data.frame(apply(Pearson_Rank_samples2, 2, FUN=function(x){mean(x, na.rm=TRUE)})), 'Mean')
Pearson_Rank_samples2_lines$SD <- apply(Pearson_Rank_samples2, 2, FUN=function(x){stats::sd(x, na.rm=TRUE)})
Pearson_Rank_samples2_lines$max <- Pearson_Rank_samples2_lines$Mean + Pearson_Rank_samples2_lines$SD 
Pearson_Rank_samples2_lines$min <- Pearson_Rank_samples2_lines$Mean - Pearson_Rank_samples2_lines$SD  
Pearson_Rank_samples2_lines$SD <- NULL

# Pearson_Rank_samples2_lines <- t(Pearson_Rank_samples2_lines)

top_ann2 <-ComplexHeatmap::HeatmapAnnotation(  Phenotype = Annotation_columns$Phenotype, 
												Sample= Annotation_columns$Sample,
												N_cells = ComplexHeatmap::anno_barplot(plotter$Number_of_cells,
													which='column', height = unit(2, "cm"), gp = grid::gpar(fill = plotter$Phenotype)),
												Signature = ComplexHeatmap::anno_boxplot(Pearson_Rank_samples2, 
													which='column', height = unit(1, "cm"), gp = grid::gpar(fill = cols), outline = FALSE),
												signature_line= ComplexHeatmap::anno_lines(Pearson_Rank_samples2_lines, gp = grid::gpar(col = c('#9d0208','#370617', '#370617'),lty=c(1, 1, 1), alpha=c(1, 0.5, 0.5))),
												col = annotation_cols)

pdf(paste0('./Plots/CopyKat/HM_Clustering_All.pdf'), width=8.25, height=11.75)
print(ComplexHeatmap::Heatmap(as.matrix(cluster_results), 
				name ='CNVA', 
				cluster_rows = FALSE, 
				column_km=8, 
				col = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
				show_column_names = FALSE,
				show_row_names = FALSE,
				show_column_dend = FALSE,
				top_annotation = top_ann2))
dev.off()
