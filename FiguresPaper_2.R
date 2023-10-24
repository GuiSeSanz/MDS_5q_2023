
library(ggplot2)
library(Seurat)


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

tmp <- as.data.frame(sc_data[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]])
table(tmp[,1], tmp[,2])

Idents(sc_data) <- 'integrated_snn_res.0.8_sub'

clusters <- setNames(as.data.frame(sc_data$integrated_snn_res.0.8_sub), c('Cluster'))


Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog', '????', 'Platelets')


MDS_sample_colors <- setNames(c('#495867', '#577399', '#bdd5ea','#f7f7ff'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))

Control_Elder_MDS_colors <- setNames(c('#3e115c','#840032', '#e59500', '#85a6b2', '#495867', '#577399', '#bdd5ea'), 
									c('Control', 'Elder', 'MDS', "SMD34459", "SMD35109", "SMD35303", "SMD37209"))



coords   <- as.data.frame(sc_data@reductions$umap@cell.embeddings)
coords_ann <- merge(setNames(sc_data[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]], c('Cluster','Cluster_names')), coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')




coords <- as.data.frame(sc_data@reductions$umap@cell.embeddings)
# annotation resolution is 0.4
coords <- merge(coords, setNames(as.data.frame(sc_data$integrated_snn_res.0.8), c('Cluster')), by=0)


Sample <- 'SMD37209'
chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMatSMD37209VsSMD29402.rds'))
coords_smp <- coords[stringr::str_detect(coords$Row.names, Sample), ]
casper_results <- reshape2::melt(chrMat)
casper_results$value2 <- "neutral"
casper_results$value2[casper_results$value > 0] <- "amplification"
casper_results$value2[casper_results$value < 0] <- "deletion"
casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
	"deletion", "neutral"))
casper_results$Var2 <- factor(casper_results$Var2, levels = colnames(chrMat))
casper_results$phenotype <- stringr::str_extract(casper_results$Var1, '^[A-Z0-9^]+')
casper_results$phenotype <- ifelse(casper_results$phenotype == 'SMD37209', 'SMD37209', 'Elder')


# pdf('./Plots/Test_casper.pdf')
casper_results_2 <- data.frame(Var1=NULL, Var2 = NULL, value = NULL, value2 = NULL, phenotype = NULL)
for (SMP in c('SMD37209', 'SMD34459', 'SMD35109', 'SMD35303')){
	message(SMP)
	chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat',SMP,'VsSMD29402.rds'))
	coords_smp <- coords[stringr::str_detect(coords$Row.names, SMP), ]
	tmp <- reshape2::melt(chrMat)
	tmp$value2 <- "neutral"
	tmp$value2[tmp$value > 0] <- "amplification"
	tmp$value2[tmp$value < 0] <- "deletion"
	tmp$value2 <- factor(tmp$value2, levels = c("amplification", 
		"deletion", "neutral"))
	tmp$Var2 <- factor(tmp$Var2, levels = colnames(chrMat))
	tmp$phenotype <- stringr::str_extract(tmp$Var1, '^[A-Z0-9^]+')
	# tmp$phenotype <- ifelse(tmp$phenotype == SMP, SMP, 'Elder')
	# Get_only the chr5
	tmp <- tmp[tmp$Var2 %in% c('5q', '5p'), ]
	if (SMP == 'SMD37209'){
		casper_results_2 <- tmp
	} else {
		casper_results_2 <- rbind(casper_results_2, tmp[tmp$phenotype== SMP, ])
	}
	# print(get_deletions_chromosomes(tmp) + facet_wrap(~phenotype, ncol=1))
}
# dev.off()

col_fun = circlize::colorRamp2(seq(-0.3, 0.3, by=(0.3-(-0.3))/(100-1)),
 , colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100))

# col_fun = circlize::colorRamp2(seq(-0.3, 0.3, by=(0.3-(-0.3))/(100-1)),
#  , colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100))

SCALED = FALSE
if (SCALED == TRUE){
	K_neighors <- 80
	CopyKat_MDS <- get_copykat_HM('SMD37209', K_neighors, anno_title=FALSE)
	dev.off()
	CopyKat_control <- get_copykat_HM('GSM5460411', K_neighors)
	dev.off()
	test_mds <- CopyKat_MDS$raw
	test_ctr <- CopyKat_control$raw
	colnames(test_mds) <- paste0('SAMPLE_', colnames(test_mds))
	colnames(test_ctr) <- paste0('CONTROL_', colnames(test_ctr))
	test <- cbind(test_mds, test_ctr)
	test_one_lin <- scale(matrix(test, nrow=ncol(test)*nrow(test), ncol=1))
	test_scaled <- matrix(test_one_lin,nrow=nrow(test), ncol=ncol(test))
	rownames(test_scaled) <- rownames(test)
	colnames(test_scaled) <- colnames(test)
	min_val <- round(min(test_scaled), 2)
	max_val <- round(max(test_scaled), 2)
	test_mds <- test_scaled[, colnames(test_scaled) %in% colnames(test_mds)]
	test_ctr <- test_scaled[, colnames(test_scaled) %in% colnames(test_ctr)]
	# test <- t(scale(t(cbind(test_mds, test_ctr))))
	# min_val <- round(min(test), 2)
	# max_val <- round(max(test), 2)
	# test_mds <- test[, colnames(test) %in% colnames(test_mds)]
	# test_ctr <- test[, colnames(test) %in% colnames(test_ctr)]
	ht_list = ComplexHeatmap::Heatmap(test_mds, name = 'SMD37209',cluster_rows = FALSE, 
						column_km=4, heatmap_legend_param = list(at = c(min_val, 0, max_val), direction = "horizontal"),
						col =col_fun,
						show_column_names = FALSE,,
						show_row_names = FALSE,
						top_annotation = CopyKat_MDS$top_anno,
						show_heatmap_legend = TRUE) + 
						ComplexHeatmap::Heatmap(test_ctr, name = 'CNA',cluster_rows = FALSE, 
						column_km=4, heatmap_legend_param = list(at = c(min_val, 0, max_val), direction = "horizontal"),
						col = col_fun,
						show_column_names = FALSE,
						show_row_names = FALSE,
						top_annotation = CopyKat_control$top_anno)
}else{
# Not scaled
	K_neighors <- 80
	CopyKat_MDS <- get_copykat_HM('SMD37209', K_neighors, anno_title=FALSE)
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
		all_bands <- stringr::str_extract(rownames(CopyKat_MDS$raw), '(?<=_)[\\w\\d\\.]+')
		position <- grep(band , all_bands)[1]
		positions <- c(positions, position)
	}

	ha = ComplexHeatmap::rowAnnotation(Locus = ComplexHeatmap::anno_mark(at = positions, labels = bands_2_show, side = "left",
	labels_gp = grid::gpar(col = "black", fontsize = 7.5), padding = 0.7))
	ht_list = ComplexHeatmap::Heatmap(CopyKat_MDS$raw,cluster_rows = FALSE, name = 'CNA Score',
						column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
						col =col_fun,
						show_column_names = FALSE,
						show_row_names = FALSE, 
						top_annotation = CopyKat_MDS$top_anno,
						left_annotation = ha,
						show_heatmap_legend = TRUE) + 
						ComplexHeatmap::Heatmap(CopyKat_control$raw,cluster_rows = FALSE, 
						column_km=4, heatmap_legend_param = list(direction = "horizontal", at = c(-0.3, 0, 0.3)),
						col = col_fun,
						show_column_names = FALSE,
						show_row_names = FALSE,
						top_annotation = CopyKat_control$top_anno,
						show_heatmap_legend = FALSE)
}



cluster_2_keep <- 1
selected_cells <- CopyKat_MDS$column_clusters[[cluster_2_keep]]
selected_cells <- CopyKat_MDS$column_names[selected_cells]
selected_cells <- selected_cells[stringr::str_detect(selected_cells, 'MDS')]
all_cells <- CopyKat_MDS$sample_distribution
selected_cells <- all_cells[all_cells$Cluster_name %in% selected_cells, 'Cell']
selected_cells <- sub('\\.', '-',selected_cells)

selected_casper <- paste(casper_results[casper_results$phenotype == 'SMD37209' &casper_results$Var2 == '5q' & casper_results$value2 == 'deletion', 'Var1'])
# data <- list(CopyKat = selected_cells, CASPER = selected_casper)
# venn_sample <- ggvenn::ggvenn(data,
# 	fill_color = destiny::cube_helix(length(data)),
# 	stroke_size = 0.4,
# 	show_percentage = TRUE,
# 	fill_alpha = 0.4,
# 	stroke_color = 'white',
# 	stroke_alpha = 1,
# 	stroke_linetype = 'solid',
# 	text_color = 'black',
# 	set_name_size = 4,
# 	text_size = 2,
# 	label_sep = ','
# 	)+ theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 7))
	
# pdf(paste0('./Plots/PaperFigures/Fig2_venn_Sample.pdf'), width=8.25, height=11.75)
# 	cowplot::plot_grid(
# 		grid::grid.grabExpr(ComplexHeatmap::draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", 
# 							annotation_legend_side = "bottom")),
# 		cowplot::plot_grid(
# 			cowplot::plot_grid(
# 				get_deletions_chromosomes(casper_results[casper_results$phenotype == 'SMD37209', ]) + theme(legend.position = 'none'),
# 				get_deletions_chromosomes(casper_results[casper_results$phenotype == 'Elder', ])+ theme(legend.position = 'none'),
# 			nrow=2),
# 			venn_sample,
# 		ncol=2),
# 	nrow=3, rel_heights=c(0.4,0.2,0.4))

# dev.off()




results_CK <- list()
plotter_list <- list()
for (Sample in c("SMD34459", "SMD35109", "SMD35303", "SMD37209")){
	tmp <- get_copykat_HM(Sample, 80)
	results_CK[[Sample]] <- tmp
	plotter_list[[Sample]] <- tmp$HM
	dev.off()
}

pdf(paste0('./Plots/PaperFigures/Fig2_SUPP.pdf'), height=10)
	cowplot::plot_grid(plotlist=plotter_list, 
	labels=names(plotter_list))
dev.off()



all_5q_selected_cells <- data.frame(Cell_id=NULL, Sample=NULL)
for (Sample in c("SMD34459", "SMD35109", "SMD35303", "SMD37209")){
	if(Sample %in% c("SMD35303", 'SMD35109')){cluster_2_keep <-  c(1,2)}else{cluster_2_keep <-  c(1)}
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


all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds')
results_CK <- readRDS('./Data/CopyKat/results_CK.rds')
plotter_list <- readRDS('./Data/CopyKat/plotter_list.rds')
cells_selected_CASPER <- readRDS('./Data/cells_selected_CASPER.rds')

# saveRDS(cells_selected_CASPER , './Data/cells_selected_CASPER.rds')
cells_selected_CASPER <- readRDS('./Data/cells_selected_CASPER.rds')

get_gvenn <- function(data, scale_Y=TRUE){
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
		axis.ticks.y=element_blank())
	if(scale_Y){
	p <- p + scale_y_continuous(limits = c(-1, 1.5))
	}
	return(p)
}

plot_list <- list()
for(Sample in c("SMD34459", "SMD35109", "SMD35303", "SMD37209")){
plot_list[[Sample]] <- get_gvenn(list(
						CopyKat = all_5q_selected_cells[all_5q_selected_cells$Sample == Sample, 'Cell_id'],
						CASPER = cells_selected_CASPER[cells_selected_CASPER$Sample == Sample, 'Cell_id']))
}

# Boxplots of genes in deleted zone 

sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
DefaultAssay(sc_data_MDS5q_subseted) <- 'SCT'
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
genes_2_check <- c('CD74', 'RPS14', 'BTF3', 'COX7C', 'HINT1', 'RPS23')

pdf('./Plots/CD74_RPS14.pdf')
VlnPlot(object = sc_data_MDS5q_subseted, features = genes_2_check, pt.size = 0,
split.by = 'cell_5q', group.by='Cluster_names')
dev.off()

plotter <- as.data.frame(sc_data_MDS5q_subseted@assays$SCT@data)
plotter <- t(plotter[genes_2_check,])
plotter <- merge(plotter, sc_data_MDS5q_subseted[[c('cell_5q', 'Cluster_names')]], by=0)

bplot <- ggplot(reshape2::melt(plotter), aes(x=Cluster_names, y=value , fill=cell_5q)) + 
geom_boxplot() + scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme_classic() + facet_wrap(~variable, ncol=1) +
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top', strip.background = element_blank(), strip.placement = "outside")


#Pseudobulk comparison
normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

data_5q    <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_elder <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

genes_2_check <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/genes_in_the_deleted_region/13-33-Table_1.tsv', sep='\t', header=T)
genes_2_check <- genes_2_check$Gene.name
genes_2_check <- genes_2_check[genes_2_check!= '']

data_5q_list <- SplitObject(data_5q, split.by='Sample')
data_5q_list <- lapply(data_5q_list, FUN=function(x){
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudo <- setNames(as.data.frame(rowSums(counts)), unique(x$Sample))
		return(pseudo)
})

data_elder_list <- SplitObject(data_elder, split.by='Sample')
data_elder_list <- lapply(data_elder_list, FUN=function(x){
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudo <- setNames(as.data.frame(rowSums(counts)), unique(x$Sample))
		return(pseudo)
})
data_5q_list <-  dplyr::bind_cols(data_5q_list)
data_elder_list <- dplyr::bind_cols(data_elder_list)
# TMP normalization
data_5q_list_tpm <- apply(data_5q_list, 2, normalize_TPM)
data_elder_list_tpm <- apply(data_elder_list, 2, normalize_TPM)

data_5q_list_tpm <- data_5q_list_tpm[rownames(data_5q_list_tpm) %in% genes_2_check, ]
data_elder_list_tpm <- data_elder_list_tpm[rownames(data_elder_list_tpm) %in% genes_2_check, ]

data_5q_list_tpm <- colSums(data_5q_list_tpm)
data_elder_list_tpm <- colSums(data_elder_list_tpm)


plotter2 <- rbind(reshape2::melt(data_5q_list_tpm), reshape2::melt(data_elder_list_tpm))
plotter2$genotype <- ifelse(stringr::str_detect(rownames(plotter2), '^SMD+'), 'MDS\n(n=4)', 'Healthy\n(n=3)')
plotter2$Sample <- rownames(plotter2)
psuedobulk <- ggplot(plotter2, aes(x=genotype, y=value, fill=genotype)) +
	labs(y ='Sum of counts in region')+
	geom_boxplot() + geom_point() + 
	# ggsignif::geom_signif(   comparisons = list(c("MDS", "Control")),
    # map_signif_level = TRUE)+
	ggrepel::geom_text_repel(aes(label=Sample), size=3) +
	scale_fill_manual(values=c('MDS\n(n=4)'='#e59500',  'Healthy\n(n=3)'='#840032')) +
	theme_classic() +
	theme(axis.title.x=element_blank(), legend.position='none', 
	axis.text.x=element_text(angle=45, vjust=1, hjust=1),
	axis.text.y=element_text(angle=90, vjust=0, hjust=0.5, size=8))



all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') # <- 
results_CK <- readRDS('./Data/CopyKat/results_CK.rds')
plotter_list <- readRDS('./Data/CopyKat/plotter_list.rds')
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')


real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
real_5q_cells <- all_5q_selected_cells[all_5q_selected_cells$Cell_id %in% real_5q_cells,]



plotter_selection <- coords_ann
plotter_selection$is_5q <- ifelse(plotter_selection$Row.names %in% real_5q_cells$Cell_id, '5q', 'Other')
plt_selection <- as.data.frame.matrix(table(plotter_selection$Sample, plotter_selection$is_5q))
pct_selection <- round(plt_selection/rowSums(plt_selection), 3)*100

prop_5q <- ggplot(plotter_selection, aes(x=Cluster_names, fill=is_5q)) + geom_bar() + facet_wrap(~Sample, nrow=1) + theme_classic()  +
		scale_fill_manual(name = 'Cell type', values=c('#586994', '#A2ABAB')) + labs(y='Number of Cells') +
		theme(legend.position='none',
		axis.text = element_text(family = "Helvetica", size = 7), 
		axis.title = element_text(family = "Helvetica", size = 9), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 9, family = "Helvetica"),
		legend.text = element_text(size = 7, family = "Helvetica"))



plotter_casper <- coords_ann
plotter_casper$is5q <- ifelse(plotter_casper$Row.names %in% cells_selected_CASPER$Cell_id, '5q', 'Other')
plt_casper <-  as.data.frame.matrix(table(plotter_casper$Sample, plotter_casper$is5q))
pct_casper <- round(plt_casper/rowSums(plt_casper), 3)*100


plotter_CopyKat <- coords_ann
plotter_CopyKat$is5q <- ifelse(plotter_CopyKat$Row.names %in% all_5q_selected_cells$Cell_id, '5q', 'Other')
plt_CopyKat <-  as.data.frame.matrix(table(plotter_CopyKat$Sample, plotter_CopyKat$is5q))
pct_CopyKat <- round(plt_CopyKat/rowSums(plt_CopyKat), 3)*100


deletion_percentages <- data.frame(Samples = c('SMD34459',	'SMD35109',	'SMD35303',	'SMD37209'), 
								Karyotype = c(43,75,90,30))
deletion_percentages <- merge(deletion_percentages, setNames(pct_selection[, '5q', drop=FALSE], c('Selected cells')), by.x= 'Samples', by.y=0)
deletion_percentages <- merge(deletion_percentages, setNames(pct_casper[, '5q', drop=FALSE], c('CASPER')), by.x= 'Samples', by.y=0)
deletion_percentages <- merge(deletion_percentages, setNames(pct_CopyKat[, '5q', drop=FALSE], c('CopyKat')), by.x= 'Samples', by.y=0)

MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))




selected_cells_Vs_Karyotype <- ggplot(reshape2::melt(deletion_percentages), aes(x = variable, y = value, group = Samples)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
			geom_hline(yintercept=deletion_percentages[, 'Karyotype'], linetype="dashed", color = "grey", size=0.5, alpha =0.5) +
			geom_line() + 
			geom_point(size = 2, aes(color = Samples)) + 
			ylim(0, 100) +
			scale_color_manual(values=MDS_sample_colors) +
			labs(y='Percentage of 5q cells') +
			theme_classic() + 
			guides(color = guide_legend(nrow=4, size=3, byrow=TRUE, title.position = 'top'))+
			theme(legend.position='none', 
			axis.text = element_text(family = "Helvetica", size = 7), 
			axis.title = element_text(family = "Helvetica", size = 9), 
			axis.text.x= element_text(angle = 45, vjust = 1, hjust=1),
			axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
			legend.title = element_text(size = 7, family = "Helvetica"),
			legend.text = element_text(size = 6, family = "Helvetica"),
			legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm'))



casper_results_2_tmp <- casper_results_2
casper_results_2_tmp$phenotype <- ifelse(casper_results_2_tmp$phenotype == 'SMD29402', 'Control', casper_results_2_tmp$phenotype)
casper_results_all_SMP <- ggplot(casper_results_2_tmp, aes(x=Var2, fill=value2)) +
	geom_bar(position='fill') + theme_classic()  +
	scale_fill_manual(values=c('neutral'="#e5e5e5", 'amplification'="#6a040f", 'deletion'="#1d3557"))+ 
	theme(legend.position='bottom', axis.title = element_text(family = "Helvetica"), 
		axis.text = element_text(family = "Helvetica", size = 8), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 9, family = "Helvetica", face = "bold"),
		legend.text = element_text(size = 8, family = "Helvetica"), 
		strip.background = element_blank(),
		strip.text= element_text(angle=45)) +
	labs(x = "Locus", y='Percent of cells') + 
	guides(fill = guide_legend(title='CNA')) + 
	facet_wrap(.~phenotype, nrow=1)



data_5q    <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

data_5q_list <- SplitObject(data_5q, split.by='Sample')
pseudoBulk <- NULL
for (sample in names(data_5q_list)){
	message(sample)
	data_5q <- data_5q_list[[sample]]
	tmp_counts <- as.data.frame(data_5q@assays$RNA@counts)
	is_5q <- data_5q[['cell_5q']]
	is_5q <- rownames(is_5q[is_5q$cell_5q=='del5q', , drop=FALSE])
	tmp_counts_del5q  <- tmp_counts[, colnames(tmp_counts) %in% is_5q]
	tmp_counts_normal <- tmp_counts[, !colnames(tmp_counts) %in% is_5q]
	tmp_counts_del5q  <- setNames(as.data.frame(rowSums(tmp_counts_del5q)), 
								paste0(unique(data_5q$Sample), '_del5q'))
	tmp_counts_normal <- setNames(as.data.frame(rowSums(tmp_counts_normal)), 
								paste0(unique(data_5q$Sample), '_normal'))
	tmp_counts_del5q  <- normalize_TPM(tmp_counts_del5q)
	tmp_counts_normal <- normalize_TPM(tmp_counts_normal)
	if(is.null(pseudoBulk)){
		pseudoBulk <- cbind(tmp_counts_del5q, tmp_counts_normal)
	}
	else{
		pseudoBulk <- cbind(pseudoBulk, tmp_counts_del5q, tmp_counts_normal)
	}
}

genes_2_check2 <- c('CD74','BTF3','COX7C','HINT1','RPS23','RPS14')
pseudoBulk <- pseudoBulk[rownames(pseudoBulk) %in% genes_2_check2, ]
pseudoBulk <- colSums(pseudoBulk)
pseudoBulk <- reshape2::melt(pseudoBulk)
pseudoBulk$Sample <- stringr::str_extract(rownames(pseudoBulk), '^[\\w\\d]+(?=_)')
pseudoBulk$genotype <- ifelse(stringr::str_detect(rownames(pseudoBulk), '_del5q'), 'del(5q)', 'Non\ndel(5q)')
pseudoBulk$genotype  <- factor(pseudoBulk$genotype, levels=c('del(5q)', 'Non\ndel(5q)'))


MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))

pseudoBulk_6_genes <- ggplot(pseudoBulk, aes(x=genotype, y=value, group=Sample)) +
	labs(y ='Sum of counts')+
	geom_line() + geom_point(aes(color=Sample)) +
	# ggrepel::geom_text_repel(aes(label=Sample), size=3) +
	scale_color_manual(values=MDS_sample_colors) +
	theme_classic() +
	theme(axis.title.x=element_blank(), 
	axis.text.x=element_text(angle=45, vjust=1, hjust=1),
	axis.text.y=element_text(angle=90, vjust=0, hjust=0.5, size=8),
	legend.position='none')


cowplot::get_legend(selected_cells_Vs_Karyotype + theme(legend.position='bottom'))
pdf('./Plots/PaperFigures/Fig2.pdf')
cowplot::plot_grid(
	cowplot::plot_grid(
			cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", annotation_legend_side = "bottom"))),
	cowplot::plot_grid(casper_results_all_SMP, 
					  get_gvenn(list( CopyKat = all_5q_selected_cells$Cell_id, CASPER = cells_selected_CASPER$Cell_id)),
					  nrow=2, rel_heights=c(1, 0.5)),
	ncol=2), 
	cowplot::plot_grid(psuedobulk, pseudoBulk_6_genes, selected_cells_Vs_Karyotype, ncol=3), 
nrow=2, rel_heights=c(1, 0.5))
dev.off()






# pdf(paste0('./Plots/PaperFigures/Fig2.pdf'), width=8.25, height=11.75/2)
# 	cowplot::plot_grid(
# 		grid::grid.grabExpr(ComplexHeatmap::draw(ht_list, merge_legend = TRUE, heatmap_legend_side = "bottom", 
# 							annotation_legend_side = "bottom")),
# 		cowplot::plot_grid(
# 			cowplot::plot_grid(
# 				get_deletions_chromosomes(casper_results[casper_results$phenotype == 'SMD37209', ]) + theme(legend.position = 'none'),
# 				get_deletions_chromosomes(casper_results[casper_results$phenotype == 'Elder', ])+ theme(legend.position = 'none'),
# 			nrow=2),
# 			# get_gvenn(list(CopyKat = all_5q_selected_cells$Cell_id, 
# 			# 		CASPER = cells_selected_CASPER$Cell_id)),
# 			get_gvenn(list( CopyKat = all_5q_selected_cells$Cell_id, 
# 				CASPER = cells_selected_CASPER$Cell_id),
# 				# Signature = rownames(Signatures_SC[Signatures_SC$Case == '5q' & Signatures_SC$Pearson_Rank > mean(Signatures_SC$Pearson_Rank),])),
# 				FALSE),
# 		ncol=2, labels= c('B', 'C')),
# 		NULL, 
# 	nrow=2, ncol=1, rel_heights=c(0.4,0.2), labels=c('A', NULL))

# dev.off()


