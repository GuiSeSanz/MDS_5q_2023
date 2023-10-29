

library(ggplot2)
library(Seurat)

# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4

get_umap_continous <- function(data, color_val){
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.8, size = 0.6) + theme_classic() + 
		viridis::scale_color_viridis( name=color_val) + 
		theme(legend.position='bottom', text = element_text(family = "Helvetica", size = 7), 
				line = element_blank(), axis.title = element_blank(), axis.text = element_blank(), 
				plot.title = element_text(hjust = 0.5)) 
		return(p)
}

get_umap <- function(data, color_val, mapped_colors, legend_row=3){
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.7, size = 0.1) + theme_classic() + 
		scale_color_manual(values=mapped_colors) + 
		theme(legend.position='bottom', 
		text = element_text(family = "Helvetica", size = 7), 
		line = element_blank(),
		title = element_blank(), 
		axis.text.x =element_blank(), 
		axis.text.y=element_blank(),  
		legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm')) + 
		guides(color = guide_legend(nrow=legend_row, byrow=TRUE, override.aes = list(size=3, alpha=1)))
	return(p)
}

get_dotplot <- function(data, features){
	data@active.ident <- factor(data@active.ident, levels=names(Cluster_colors))
	p <- DotPlot(data, features = features, col.min=0) + coord_flip() +
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	# viridis::scale_colour_viridis(option="plasma") +
	# scale_colour_gradient(low = 'blue', high='red') +
    scale_size_continuous(range = c(0.1,4))+
	scale_colour_gradient2(low = scales::muted("blue"), high = '#db3e25', mid ='white') +
	theme(legend.position='right', axis.title = element_text(family = "Helvetica", size = 7), 
	 	axis.text = element_text(family = "Helvetica", size = 7), 
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
																which='column', height = unit(2, "cm"), gp = grid::gpar(fill = cols), ylim =c(0,1500)),
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

# ===================================================
# FIGURE 1 <- SINGLE CELL DESCRIPTION
# ===================================================


# sc_data <- readRDS('./Data/all_seurat_integrated_sct.rds')
# sc_data <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_integrated_sct_subcluster.rds')
SAMPLE_NAME <- '5qSamples'
sc_data <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', SAMPLE_NAME, '_Annotated_final.rds'))

tmp <- as.data.frame(sc_data[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]])
table(tmp[,1], tmp[,2])

Idents(sc_data) <- 'integrated_snn_res.0.8_sub'

clusters <- setNames(as.data.frame(sc_data$integrated_snn_res.0.8_sub), c('Cluster'))
# clusters_12 <- setNames(as.data.frame(sc_data$integrated_snn_res.0.8_sub12), c('Cluster'))
# clusters_17 <- setNames(as.data.frame(sc_data$integrated_snn_res.0.8_sub17), c('Cluster'))
# clusters_12 <- clusters_12[clusters_12$Cluster %in% c('12_0', '12_2'),, drop=FALSE]
# clusters_17 <- clusters_17[clusters_17$Cluster %in% c('17_1'),, drop=FALSE]
# clusters_sub <- rbind(clusters_12, clusters_17)
# # remove the reassigned cells
# clusters <- clusters[!rownames(clusters) %in% rownames(clusters_sub),, drop=FALSE]
# clusters <- clusters[!clusters$Cluster %in% c(12, 17),, drop=FALSE]
# # add the reassignde cells 
# clusters <- rbind(clusters, clusters_sub)

# tmp <- RunUMAP(sc_data,  features = VariableFeatures(object = sc_data))


# Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722')
# names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog')


Cluster_colors <- c(
	'#A4CEE2',
	'#1F77B5',
	'#2B6E38', 
	'#BBD58D', 
	'#842758', 
	'#E1BB6D', 
	'#C3B1CC', 
	'#673F95', 
	'#010000', 
	'#E99C9A', 
	'#E0A630', 
	'#D33C28', 
	'#992223', 
	'#A76031')


names(Cluster_colors) <- c(
	'HSC', 
	'LMPP', 
	'GMP', 
	'Granulocyte',
	'Monocytes', 
	'DendriticCell', 
	'CLP', 
	'pro-B', 
	'T', 
	'MEP', 
	'MK_Prog', 
	'EarlyErythroid',
	'LateErythroid', 
	'Basophil')




cellType_legend_order <- c()

MDS_sample_colors <- setNames(c('#495867', '#577399', '#bdd5ea','#f7f7ff'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("Patient_1", "Patient_2", "Patient_3", "Patient_4"))
# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4

# MDS_sample_colors <- c('#355070', '#6d597a', '#b56576', '#e56b6f')
# MDS_sample_colors <- c('#f79f79', '#f7d08a', '#e3f09b', '#87b6a7')
# MDS_sample_colors <- c('#ed6a5a', '#f4f1bb', '#9bc1bc', '#5d576b')
# MDS_sample_colors <- c('#e8d6cb', '#d0ada7', '#ad6a6c', '#5d2e46')

Control_Elder_MDS_colors <- setNames(c('#002642','#840032', '#e59500'), c('Control', 'Elder', 'MDS'))

coords   <- as.data.frame(sc_data@reductions$umap@cell.embeddings)
coords_ann <- merge(setNames(sc_data[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]], c('Cluster','Cluster_names')), coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')


patient_dist_cell_count <- ggplot(coords_ann, aes(x=Cluster_names, fill=Sample)) + 
		geom_bar( position='stack', colour="black") + 
		theme_classic() + scale_fill_manual(values=MDS_sample_colors, name='Samples') + 
		theme(axis.text=element_text(family = "Helvetica", size=8),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=8),
			legend.title=element_text(size=7),
			panel.grid.major.y = element_line( size=.1, color="grey" ) ) + 
		labs(y = "Number of Cells")




Idents(sc_data) <- 'Cluster_names'
DefaultAssay(sc_data) <- "SCT"
markers <- read.table('./Data/Annotation/encyclopedia_table.tsv', sep='\t', header=TRUE)

markers_2_plot <- c('Erythrocytes', 'Monocytes', 'HSC/MPP', 'GMP', "Pro-B cells", "Early erythroid", "Late erythroid" ,  "Granulocytes", "DC")
features_2_plot <- strsplit(as.character(markers[ markers$Low.hierarchy.cell.types %in% markers_2_plot, 'Curated.markers']), split=',')
features_2_plot <- unique(trimws(unlist(features_2_plot)))
features_2_plot <- c('AVP', 'HOPX', 'CRHBP', 'SATB1', 'SPINK2', 'LSP1', 'FLT3','MPO', 'CTSG', 'PRTN3','AZU1', 'ELANE','FCER1G', 'CSTA', 'LYZ', 'IL3RA', 'IRF7', 'IRF8','JCHAIN', 'IKZF1','DNTT', 'ADA', 'LTB','CD79A', 'EBF1', 'VPREB1','PBX1', 'PLEK', 'DAD1','CNRIP1', 'SLC40A1', 'PKIG','KLF1', 'AHSP', 'HBB','RUNX1', 'HDC', 'MS4A3')





# proliferation <- cbind(setNames(as.data.frame(sc_data$S.Score), c('S Score')), setNames(as.data.frame(sc_data$G2M.Score), c('G2M Score')))
# proliferation <- merge(proliferation, coords_ann, by.x=0, by.y='Row.names')


# proliferation_s <- proliferation[sample(nrow(proliferation), 12000),]
# proliferation_s <- proliferation_s[order(proliferation_s$`S Score`),]
# proliferation_s   <- get_umap_continous(proliferation_s,   'S Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))
# proliferation_g2m <- proliferation[sample(nrow(proliferation), 12000),]
# proliferation_g2m <- proliferation[order(proliferation_g2m$`G2M Score`),]
# proliferation_g2m <- get_umap_continous(proliferation_g2m, 'G2M Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))

# pdf('./Plots/Integrated/subclusters_cicle.pdf')
# get_umap_continous(proliferation[proliferation$Cluster %in% c('17_1'),],   'S Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE)) + 
# geom_point(data = proliferation[!proliferation$Cluster %in% c('17_1'),], aes(x=UMAP_1, y=UMAP_2),size=0.2, alpha=0.3, color='grey') + ggtitle('Cluster 17_1 S Score')
# get_umap_continous(proliferation[proliferation$Cluster %in% c('17_1'),], 'G2M Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))+ 
# geom_point(data = proliferation[!proliferation$Cluster %in% c('17_1'),], aes(x=UMAP_1, y=UMAP_2),size=0.2, alpha=0.3, color='grey')+ ggtitle('Cluster 17_1 G2M Score')
# dev.off()


# pdf('./Plots/PaperFigures/Fig1_2.pdf', width=8.25, height=11.75)
# cowplot::plot_grid(
# 	cowplot::plot_grid(
# 		get_umap(coords_ann, 'Cluster_names', Cluster_colors), get_dotplot(sc_data_subseted , features_2_plot),
# 	ncol = 2, rel_widths = c(0.5, 0.5)
# 	), 
# 	cowplot::plot_grid(
# 		cowplot::plot_grid(
# 			get_umap(coords_ann, 'Sample', MDS_sample_colors, legend_row=1), patient_dist_cell_count + theme(legend.position = 'none'),
# 			nrow=2),
# 		cowplot::plot_grid(proliferation_s, proliferation_g2m, ncol=1)
# 	),
# nrow=3, rel_heights = c(0.3, 0.3, 0.3)
# )

# dev.off()



# sample_umaps <- list()
# for (Sample in  c("SMD34459", "SMD35109", "SMD35303", "SMD37209")){
# 	tmp <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', Sample, '/', Sample, '_seurat_obj_norm.rds'))
# 	tmp_coords <- as.data.frame(tmp@reductions$umap@cell.embeddings)
# 	tmp_coords <- tmp_coords[rownames(tmp_coords) %in% coords_ann$Row.names,]
# 	tmp_coords <- merge(tmp_coords, coords_ann[, c('Row.names', 'Cluster_names')], by.x=0, by.y='Row.names')
# 	sample_umaps[[Sample]] <- get_umap(tmp_coords, 'Cluster_names', Cluster_colors) + 
# 							  theme(legend.position= 'none', plot.title = element_text(hjust = 0.5)) + 
# 							  ggtitle(Sample)
# }


# pdf('./Plots/PaperFigures/Fig1_1.pdf', width=8.25, height=11.75)
# cowplot::plot_grid(
# 	cowplot::plot_grid(
# 		get_umap(coords_ann, 'Cluster_names', Cluster_colors), get_dotplot(sc_data_subseted , features_2_plot),
# 	ncol = 2, rel_widths = c(0.5, 0.5)), 
# 	cowplot::plot_grid(plotlist=sample_umaps, nrow=1),
# 	cowplot::plot_grid(
# 			get_umap(coords_ann, 'Sample', MDS_sample_colors, legend_row=1), patient_dist_cell_count + theme(legend.position = 'none'),
# 	nrow=1),
# 	cowplot::plot_grid(proliferation_s, proliferation_g2m, nrow=1),
# nrow=4, rel_heights = c(1, 0.75,1,1))

# dev.off()



elder_integrated <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
tmp <- as.data.frame(elder_integrated[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]])
table(tmp[,1], tmp[,2])


populations <- rbind(setNames(as.data.frame(elder_integrated$Cluster_names), 'Cluster_names'), 
					 setNames(as.data.frame(sc_data$Cluster_names), 'Cluster_names'))

populations$Sample <- stringr::str_extract(rownames(populations),'(^[A-Z0-9]+)(?=_|-)')
pop_prop <- table(populations$Cluster_names, populations$Sample)
pop_prop <- apply(pop_prop, 2, FUN=function(x) (x/sum(x))*100)
pop_prop <- as.data.frame(reshape2::melt(pop_prop))
pop_prop$Case <- ifelse(stringr::str_detect(pop_prop$Var2, '^SMD'), 'MDS_5q', 'Elder')


for (pop in unique(pop_prop$Var1)){
	result <- wilcox.test(pop_prop[pop_prop$Var1 == pop & pop_prop$Case == 'MDS_5q', 'value'], pop_prop[pop_prop$Var1 == pop & pop_prop$Case == 'Elder', 'value'])
	if (result$p.value < 0.05){
		print(pop)
		print('Significant')
		print(result)
	}
}

data_summary <- function(data, varname, groupnames){
	#+++++++++++++++++++++++++
	# Function to calculate the mean and the standard deviation
	# for each group
	#+++++++++++++++++++++++++
	# data : a data frame
	# varname : the name of a column containing the variable
	#to be summariezed
	# groupnames : vector of column names to be used as
	# grouping variables
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

pop_prop$Case <- ifelse(pop_prop$Case == 'MDS_5q', 'del(5q) MDS', 'Healthy Sample')
pop_prop <- data_summary(pop_prop, 'value', c('Var1', 'Case'))

# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4

coords_ann$Sample <- ifelse(coords_ann$Sample == 'SMD34459', 'Patient_1', 
					 ifelse(coords_ann$Sample == 'SMD35109', 'Patient_2',
					 ifelse(coords_ann$Sample == 'SMD35303', 'Patient_3',
					 ifelse(coords_ann$Sample == 'SMD37209', 'Patient_4','None'))))

tmp_data <- coords_ann
tmp_data$Cluster_names <- factor(tmp_data$Cluster_names, levels=names(Cluster_colors))
patient_dist_cell_count <- ggplot(tmp_data, aes(x=Cluster_names, fill=Sample)) + 
		geom_bar( position='stack', colour="black") + 
		theme_classic() + scale_fill_manual(values=MDS_sample_colors, name='Samples') + 
		theme(axis.text=element_text(family = "Helvetica", size=7, color='black'),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=7),
			legend.title=element_text(size=9),
			legend.key.size=unit(0.5, "cm"),
			panel.grid.major.y = element_line( size=.1, color="grey" ),
			legend.position = 'none' ) + 
		labs(y = "Number of Cells")


# get the UMAPS per sample

sample_umaps <- list()
for (sample_name in c('SMD34459','SMD35109','SMD35303','SMD37209')){
	print(sample_name)
	tmp <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', sample_name, '/', sample_name, '_seurat_obj_norm.rds'))
	tmp <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2'))
	tmp$Sample <- sample_name
	sample_name <- ifelse(sample_name == 'SMD34459', 'Patient_1', 
				   ifelse(sample_name == 'SMD35109', 'Patient_2',
				   ifelse(sample_name == 'SMD35303', 'Patient_3',
				   ifelse(sample_name == 'SMD37209', 'Patient_4','None'))))
	tmp <- merge(tmp, coords_ann[, c('Row.names', 'Cluster_names')], by.x=0, by.y='Row.names')
	p <- get_umap(tmp, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 3) 
	p <- p + theme(legend.position= 'none', plot.title = element_text(hjust = 0.5)) + ggtitle(sample_name)
	sample_umaps[[sample_name]] <- p
}



# pct_cluster_per_sample <- ggplot(reshape2::melt(table(coords_ann$Sample, coords_ann$Cluster_names)),aes(x= Var2, y= value, fill= Var1)) + 
# geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=MDS_sample_colors, name='Samples') +
# theme_classic() + 
# theme(axis.text=element_text(family = "Helvetica", size=7),
# 	axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
# 	axis.title.x = element_blank(),
# 	axis.title = element_text(size=7, family = "Helvetica"),
# 	axis.line = element_line(colour = 'black', size = 0.5),
# 	legend.text=element_text(size=7),
# 	legend.title=element_text(size=7),
# 	legend.key.size=unit(0.5, "cm"),
# 	panel.grid.major.y = element_line( size=.1, color="grey" ),
# 	legend.position = 'none' ) + 
# labs(y = "Percentage of cells")

tmp_data <- reshape2::melt(table(coords_ann$Sample, coords_ann$Cluster_names))
tmp_data$Var2 <- factor(tmp_data$Var2, levels=c(rev(names(Cluster_colors))))
pct_cluster_per_sample <- ggplot(tmp_data,aes(x= value, y= Var1, fill= Var2, ,order=Var2)) + 
geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=Cluster_colors, name='Samples') +
theme_classic() + 
theme(axis.text=element_text(family = "Helvetica", size=7, color='black'),
	# axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
	axis.title.y = element_blank(),
	axis.title = element_text(size=7, family = "Helvetica"),
	axis.line = element_line(colour = 'black', size = 0.5),
	legend.text=element_text(size=7, color='black'),
	legend.title=element_text(size=9),
	legend.key.size=unit(0.5, "cm"),
	panel.grid.major.y = element_line( size=.1, color="grey" ),
	legend.position = 'none' ) + 
labs(x = "Proportion of cells")


tmp_data <- pop_prop
tmp_data$Var1 <- factor(tmp_data$Var1, levels=c(names(Cluster_colors)))
barplot_populations <- ggplot(tmp_data, aes(x=Var1, y=value, fill=Case)) + theme_classic() +
geom_bar(stat = "identity", position = "dodge", colour = "black", alpha = 0.8) + 
geom_errorbar(aes(ymin = value - sd, ymax = value + sd), width=.2, position=position_dodge(.9)) +
scale_fill_manual(values=c('#840032', '#0077B6'), name='Condition') +
		theme(axis.text=element_text(family = "Helvetica", size=7, color='black'),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1, color='black'), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica", color='black' ),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=7),
			legend.title=element_text(size=9),
			panel.grid.major.y = element_line( size=.1, color="grey" ),
			legend.position = 'none') + 
		labs(y = "Percentage of  of Cells")


# # mid size figure width=5.82, height=8.25
# # full size figure width=8.25, height=11.6
# pdf('./Plots/PaperFigures/Fig1.pdf', width=5.82, height=8.25)
# cowplot::plot_grid(
# 	cowplot::plot_grid(
# 		get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 4), 
# 		patient_dist_cell_count, 
# 	ncol=1, rel_heights = c(1, 0.8), labels = c('A', 'C')
# 	),
# 	cowplot::plot_grid(get_dotplot(sc_data , features_2_plot), 
# 						barplot_populations, 
# 	ncol=1, rel_heights = c(1, 0.8), labels = c('B', 'D')), 
# ncol=2)
# dev.off()




# full size figure width=8.3, height=11.7

sample_legend <- cowplot::get_legend(patient_dist_cell_count + theme(legend.position='bottom', legend.direction = "horizontal", legend.background=element_blank()) + 
guides(fill = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=0.5, alpha=0.9))))
elder_mds_legend <- cowplot::get_legend(barplot_populations + theme(legend.position='bottom', legend.direction = "horizontal", legend.background=element_blank()) + 
guides(fill = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=0.3, alpha=0.9))))

# pdf('./Plots/PaperFigures/Fig1.pdf', width=8.3, height=10)
# cowplot::plot_grid(
# 	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 3), 
# 					get_dotplot(sc_data , features_2_plot), 
# 	ncol=2, labels=c('A', 'B')),
# 	cowplot::plot_grid(
# 		cowplot::plot_grid(pct_cluster_per_sample, patient_dist_cell_count, barplot_populations, ncol=3, labels=c('C', 'D', 'E'),align='h' ), 
# 		cowplot::plot_grid(sample_legend, elder_mds_legend , ncol=2, rel_widths=c(2,1) ), 
# 		nrow=2, rel_heights=c(1,0.2)),
# nrow=2)
# dev.off()


# pdf('./Plots/PaperFigures/Fig1.pdf', width=8.3, height=10)
# cowplot::plot_grid(
# 	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 3), 
# 					get_dotplot(sc_data , features_2_plot), 
# 	ncol=2, labels=c('A', 'B')),
# 	cowplot::plot_grid(plotlist = sample_umaps, nrow=1, labels='C'),
# 	cowplot::plot_grid(
# 		cowplot::plot_grid(pct_cluster_per_sample, patient_dist_cell_count, barplot_populations, ncol=3, labels=c('D', 'E', 'F'),align='h' ), 
# 		cowplot::plot_grid(sample_legend, elder_mds_legend , ncol=2, rel_widths=c(2,1) ), 
# 		nrow=2, rel_heights=c(1,0.2)),
# nrow=3, rel_heights=c(1,0.5,0.75))
# dev.off()

pdf('./Plots/PaperFigures/Fig1.pdf', width=8.3, height=10)
cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_image("/home/tereshkova/data/gserranos/MDS/Plots/PaperFigures/Untitled-27.png", scale = 0.9)
,
cowplot::plot_grid(
	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 3), 
					cowplot::plot_grid(plotlist = sample_umaps, nrow=2),
	ncol=2, labels=c('A', 'B')),
	cowplot::plot_grid(
		cowplot::plot_grid(
			get_dotplot(sc_data , features_2_plot),
			cowplot::plot_grid(pct_cluster_per_sample ,
				cowplot::plot_grid(
					cowplot::plot_grid( patient_dist_cell_count, barplot_populations, ncol=2, labels = c('E', 'F')),
					cowplot::plot_grid(sample_legend, elder_mds_legend ,NULL,  ncol=3, rel_widths=c(2,1, 0.1) ),
				nrow=2, rel_heights=c(1,0.1), align='v'),
			nrow=2, rel_heights=c(0.3, 0.6)),
			ncol=2, labels=c('C', 'D'))
		),
nrow=2, rel_heights=c(0.7,1))
, nrow=2, rel_heights = c(0.1,0.9))

dev.off()



# newdotplotgenes
# hsc<-c("CRHBP", "HOPX", "KYT", "CD34")
# lmpp<-c("PTPRC", "FLT3", "PROM1", "SATB1")
# cc<-c("CDC20", "TOP2A")
# gmp<-c("CSF3R", "CTSG", "PRTN3", "MPO")
# granul<-c("ELANE", "AZU1", "CEBPA", "CST7")
# mono<-c("LYZ", "CSTA")
# dc<-c("IRF8", "IRF7", "IL3RA", "CLEC4", "IGKC", "SPIB", "BLINK")
# t<-c("JCHAIN", "PRSS2", "TXNP", "LTB")
# nk<-c("CXXC5", "HOXA9", "HOXA10")
# clp<-c("IL7R", "DNTT")
# prob<-c("VPREB1", "EBF1", "CD79A", "CD79B", "IGHM")
# mep<-c("NFE2", "HFS1", "TAL1")
# mk<-c("PBX1", "VWF", "FLI1", "ITGA22B", "GP1BA")
# ery<-c("GATA1", "HBD", "HBB", "CA1", "AHSP",  "KLF1")
# baso<-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")