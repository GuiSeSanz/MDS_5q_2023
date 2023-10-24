

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
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.7, size = 0.5) + theme_classic() + 
		scale_color_manual(values=mapped_colors) + 
		theme(legend.position='bottom', text = element_text(family = "Helvetica", size = 7), 
		line = element_blank(),title = element_blank(), axis.text.x =element_blank(), axis.text.y=element_blank()) + 
		guides(color = guide_legend(nrow=legend_row, byrow=TRUE, override.aes = list(size=3, alpha=0.9)))
	return(p)
}

get_dotplot <- function(data, features){
	p <- DotPlot(data, features = features, col.min=0) + coord_flip() +
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	# viridis::scale_colour_viridis(option="plasma") +
	# scale_colour_gradient(low = 'blue', high='red') +
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

	Control_Elder_MDS_colors <- setNames(c('#002642','#840032', '#e59500'), c('Control', 'Elder', 'MDS'))
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
													# which='column', height = unit(1.5, "cm"), gp = grid::gpar(fill = cols), outline = FALSE),
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
# FIGURE S1 <- SINGLE CELL DESCRIPTION ELDER INTEGRATED
# ===================================================


annotation_Sofia_08 <- list()
annotation_Sofia_08[["0"   ]] <- 'HSC'
annotation_Sofia_08[["1"   ]] <- 'HSC'
annotation_Sofia_08[["2"   ]] <- 'HSC'
annotation_Sofia_08[["3"   ]] <- 'Granulocite'
annotation_Sofia_08[["4"   ]] <- 'EarlyErythroid'
# annotation_Sofia_08[["5_0" ]] <- 
# annotation_Sofia_08[["6_0" ]] <- 
# annotation_Sofia_08[["6_1" ]] <- 
annotation_Sofia_08[["7"   ]] <- 'EarlyErythroid'
annotation_Sofia_08[["8_0" ]] <- 'CLP'
annotation_Sofia_08[["8_1" ]] <- 'LMPP'
annotation_Sofia_08[["8_2" ]] <- 'LMPP'
annotation_Sofia_08[["9"   ]] <- 'LateErythroid'
annotation_Sofia_08[["10"  ]] <- 'GMP'
annotation_Sofia_08[["11"  ]] <- 'MEP'
# annotation_Sofia_08[["12_0"]] <- 
annotation_Sofia_08[["12_1"]] <- 'Granulocite'
annotation_Sofia_08[["12_2"]] <- 'MEP'
# annotation_Sofia_08[["13_0"]] <- 
# annotation_Sofia_08[["13_1"]] <- 
annotation_Sofia_08[["14_0"]] <- 'Basophil'
# annotation_Sofia_08[["14_1"]] <- 
annotation_Sofia_08[["14_2"]] <- 'Basophil'
annotation_Sofia_08[["15"  ]] <- 'MEP'
annotation_Sofia_08[["16"  ]] <- 'EarlyErythroid'
annotation_Sofia_08[["17"  ]] <- 'LateErythroid'
annotation_Sofia_08[["18_0"]] <- 'DendriticCell'
annotation_Sofia_08[["18_1"]] <- 'DendriticCell'
annotation_Sofia_08[["19"  ]] <- 'Monocytes'
# annotation_Sofia_08[["20_0"]] <- 
annotation_Sofia_08[["21_0"]] <- 'pro-B'
# annotation_Sofia_08[["22"  ]] <- 
annotation_Sofia_08[["23"  ]] <- 'LateErythroid'




elder_integrated <-  readRDS(paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_subcluster.rds'))
# saveRDS(elder_integrated,'./Data/all_seurat_Elder_integrated_sct_NO_MITO_FILTER_1JUL.rds')

# this metric with give us an idea of the complexity of our dataset
# elder_integrated$mitoRatio <- PercentageFeatureSet(object = elder_integrated, features =c(base::grep('^MT', rownames(elder_integrated), value=T)))
# elder_integrated$mitoRatio <- elder_integrated$mitoRatio/100
# elder_integrated <- subset(elder_integrated, subset=mitoRatio < 0.15)
# saveRDS(elder_integrated, './Data/all_seurat_Elder_integrated_sct.rds')

clusters <- setNames(as.data.frame(elder_integrated$integrated_snn_res.0.8_sub), c('Cluster'))

DefaultAssay(elder_integrated)  <- "RNA"
Idents(elder_integrated) <- 'integrated_snn_res.0.8_sub'


pdf('./Plots/Integrated/Test.pdf')
plotter <- merge(setNames(as.data.frame(elder_integrated$integrated_snn_res.0.8_sub), 'Cluster'),setNames(as.data.frame(elder_integrated$S.Score), 'S.score'), by=0)
plotter <- merge(plotter, setNames(as.data.frame(elder_integrated$G2M.Score)    , 'G2M.Score'), by.x= 'Row.names', by.y=0)
plotter <- merge(plotter, setNames(as.data.frame(filtered_seurat$mitoRatio)    , 'mitoRatio'), by.x= 'Row.names', by.y=0)
plotter <- merge(plotter, setNames(as.data.frame(elder_integrated$RPSRatio)     , 'RPSRatio'), by.x= 'Row.names',  by.y=0)
ggplot(reshape2::melt(plotter), aes(x=Cluster, y=value, fill = Cluster)) + geom_violin()+ theme_classic() + theme(legend.position='none') + facet_wrap(~variable, scales='free')
dev.off()




coords   <- as.data.frame(elder_integrated@reductions$umap@cell.embeddings)

Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78')
# Cluster_colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(unique(unlist(annotation_Sofia_08))))
# names(Cluster_colors) <- unique(unlist(annotation_Sofia_08))
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocite","CLP","MEP","Basophil","T","DendriticCell")
clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_08),, drop=FALSE]
coords_ann <- merge(clusters, coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_08[as.character(x[['Cluster']])][[1]] })
ElderSample_colors <- c( '#264653', '#2a9d8f', '#e9c46a')
names(ElderSample_colors) <- unique(coords_ann$Sample)





clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_08),, drop=FALSE]
coords_ann <- merge(clusters, coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_08[as.character(x[['Cluster']])][[1]] })




patient_dist_cell_count <- ggplot(coords_ann, aes(x=Cluster_names, fill=Sample)) + 
		geom_bar( position='stack', colour="black") + 
		theme_classic() + scale_fill_manual(values=ElderSample_colors, name='Samples') + 
		theme(axis.text=element_text(family = "Helvetica", size=6),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=7),
			legend.title=element_text(size=7),
			panel.grid.major.y = element_line( size=.1, color="grey" ) ) + 
		labs(y = "# of Cells")

sc_data_subseted <- subset(elder_integrated, idents = c('5_0', '6_0', '6_1', '12_0', '13_0', '13_1', '14_1', '20_0', '22'), invert = TRUE)




Clusters <- setNames(as.data.frame(as.character(sc_data_subseted$integrated_snn_res.0.8_sub)), c('Cluster'))
Clusters$Annotated_Clusters <- apply(Clusters, 1, function(x){ annotation_Sofia_08[as.character(x[['Cluster']])][[1]] })
rownames(Clusters) <- names(sc_data_subseted$integrated_snn_res.0.8_sub)

sc_data_subseted$Annotated_Clusters <- Clusters$Annotated_Clusters

Idents(sc_data_subseted) <- 'Annotated_Clusters'

markers <- read.table('./Data/Annotation/encyclopedia_table.tsv', sep='\t', header=TRUE)

markers_2_plot <- c('Erythrocytes', 'Monocytes', 'HSC/MPP', 'GMP', "Pro-B cells", "Early erythroid", "Late erythroid" ,  "Granulocytes", "DC")
features_2_plot <- strsplit(as.character(markers[ markers$Low.hierarchy.cell.types %in% markers_2_plot, 'Curated.markers']), split=',')
features_2_plot <- unique(trimws(unlist(features_2_plot)))

features_2_plot <- c('AVP', 'HOPX', 'CRHBP', 'SATB1', 'LSP1', 'FLT3','MPO', 'CTSG', 'PRTN3','AZU1', 'ELANE','FCER1G', 'CSTA', 'LYZ', 'IL3RA', 'IRF7', 'IRF8','JCHAIN', 'IKZF1','DNTT', 'ADA', 'LTB','CD79A', 'EBF1', 'VPREB1','PBX1', 'PLEK', 'DAD1','CNRIP1', 'SLC40A1', 'PKIG','KLF1', 'AHSP', 'HBB','RUNX1', 'HDC', 'MS4A3')



proliferation <- cbind(setNames(as.data.frame(sc_data_subseted$S.Score), c('S Score')), setNames(as.data.frame(sc_data_subseted$G2M.Score), c('G2M Score')))
proliferation <- merge(proliferation, coords_ann, by.x=0, by.y='Row.names')

proliferation_s <- proliferation[sample(nrow(proliferation), 12000),]
proliferation_s <- proliferation_s[order(proliferation_s$`S Score`),]
proliferation_s   <- get_umap_continous(proliferation_s,   'S Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))
proliferation_g2m <- proliferation[sample(nrow(proliferation), 12000),]
proliferation_g2m <- proliferation[order(proliferation_g2m$`G2M Score`),]
proliferation_g2m <- get_umap_continous(proliferation_g2m, 'G2M Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))


sample_umaps <- list()
for (Sample in  unique(sc_data_subseted$Sample)){
	tmp <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/Normal_Data/', Sample, '_seurat_obj_norm.rds'))
	tmp_coords <- as.data.frame(tmp@reductions$umap@cell.embeddings)
	tmp_coords <- tmp_coords[rownames(tmp_coords) %in% coords_ann$Row.names,]
	tmp_coords <- merge(tmp_coords, coords_ann[, c('Row.names', 'Cluster_names')], by.x=0, by.y='Row.names')
	sample_umaps[[Sample]] <- get_umap(tmp_coords, 'Cluster_names', Cluster_colors) + 
							  theme(legend.position= 'none', plot.title = element_text(hjust = 0.5)) + 
							  ggtitle(Sample)
}


pdf('./Plots/PaperFigures/FigS1.pdf', width=8.25, height=11.75)
cowplot::plot_grid(

	cowplot::plot_grid(
		get_umap(coords_ann, 'Cluster_names', Cluster_colors), 
		cowplot::plot_grid(plotlist=sample_umaps, nrow=2), 
		proliferation_s,
	ncol=1, rel_heights = c(1, 0.8, 0.5)
	),
	cowplot::plot_grid(get_dotplot(sc_data_subseted , features_2_plot), 
						patient_dist_cell_count + theme(legend.position = 'top'), 
						proliferation_g2m, 
	ncol=1, rel_heights = c(1, 0.8, 0.5)),
ncol=2)

dev.off()




# ===================================================
# FIGURE S2 <- SINGLE CELL DESCRIPTION PrePost Treatment Integrated
# ===================================================


# HSC: 10
# LMPP: 11,1
# GMP: 13
# GRANULOCITES: 8
# MONOCITES: 3
# DC: NO HAY
# T: 9
# CLP: 9
# PROB: 9
# MEP: 5,6,17_2
# ERY-EARLY: 0,12,14
# ERY-LATE: 2,4,7,16,17_1
# BASO:18
# PLAQUETAS:21
# UNKNOWN: 19



annotation_Nerea_08 <- list()
annotation_Nerea_08[["0"   ]] <- 'EarlyErythroid'
annotation_Nerea_08[["1"   ]] <- 'LMPP'
annotation_Nerea_08[["2"   ]] <- 'LateErythroid'
annotation_Nerea_08[["3"   ]] <- 'Monocytes'
annotation_Nerea_08[["4"   ]] <- 'LateErythroid'
annotation_Nerea_08[["5"   ]] <- 'MEP'
annotation_Nerea_08[["6"   ]] <- 'MEP'
annotation_Nerea_08[["7"   ]] <- 'EarlyErythroid'
annotation_Nerea_08[["8"   ]] <- 'Granulocite'
annotation_Nerea_08[["9_0" ]] <- 'CLP'
annotation_Nerea_08[["9_1" ]] <- 'CLP'
annotation_Nerea_08[["10"  ]] <- 'HSC'
annotation_Nerea_08[["11"  ]] <- 'LMPP'
annotation_Nerea_08[["12"  ]] <- 'EarlyErythroid'
annotation_Nerea_08[["13"  ]] <- 'GMP'
annotation_Nerea_08[["14"  ]] <- 'EarlyErythroid'
annotation_Nerea_08[["16"  ]] <- 'EarlyErythroid'
annotation_Nerea_08[["17_1"]] <- 'EarlyErythroid'
annotation_Nerea_08[["17_2"]] <- 'MEP'
annotation_Nerea_08[["18"  ]] <- 'Basophil'
annotation_Nerea_08[["19_0"]] <- '????'
annotation_Nerea_08[["19_1"]] <- '????'
annotation_Nerea_08[["21"  ]] <- 'Platelets'






pre_post_5q <- readRDS(paste0(getwd(), '/Data/','PrePost_integratedsct_subcluster.rds'))

DefaultAssay(pre_post_5q)  <- "RNA"
Idents(pre_post_5q) <- 'integrated_snn_res.0.8_sub'



coords   <- as.data.frame(pre_post_5q@reductions$umap@cell.embeddings)

Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#03071e', '#006d77')
# Cluster_colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(unique(unlist(annotation_Sofia_08))))
# names(Cluster_colors) <- unique(unlist(annotation_Sofia_08))
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocite","CLP","MEP","Basophil","T","DendriticCell", '????', 'Platelets')

clusters <- setNames(as.data.frame(pre_post_5q$integrated_snn_res.0.8_sub), c('Cluster'))

clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Nerea_08),, drop=FALSE]
coords_ann <- merge(clusters, coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Nerea_08[as.character(x[['Cluster']])][[1]] })
PrePostSample_colors <- c( '#606c38', '#dda15e')
names(PrePostSample_colors) <- unique(coords_ann$Sample)



patient_dist_cell_count <- ggplot(coords_ann, aes(x=Cluster_names, fill=Sample)) + 
		geom_bar( position='stack', colour="black") + 
		theme_classic() + scale_fill_manual(values=PrePostSample_colors, name='Samples') + 
		theme(axis.text=element_text(family = "Helvetica", size=6),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=7),
			legend.title=element_text(size=7),
			panel.grid.major.y = element_line( size=.1, color="grey" ) ) + 
		labs(y = "# of Cells")

#remove the clusters
sc_data_subseted <- subset(pre_post_5q, idents = c('15_0', '15_1', '17_0', '20'), invert = TRUE)

Clusters <- setNames(as.data.frame(as.character(sc_data_subseted$integrated_snn_res.0.8_sub)), c('Cluster'))
Clusters$Annotated_Clusters <- apply(Clusters, 1, function(x){ annotation_Nerea_08[as.character(x[['Cluster']])][[1]] })
rownames(Clusters) <- names(sc_data_subseted$integrated_snn_res.0.8_sub)

sc_data_subseted$Annotated_Clusters <- Clusters$Annotated_Clusters

Idents(sc_data_subseted) <- 'Annotated_Clusters'

markers <- read.table('./Data/Annotation/encyclopedia_table.tsv', sep='\t', header=TRUE)

markers_2_plot <- c('Erythrocytes', 'Monocytes', 'HSC/MPP', 'GMP', "Pro-B cells", "Early erythroid", "Late erythroid" ,  "Granulocytes", "DC")
features_2_plot <- strsplit(as.character(markers[ markers$Low.hierarchy.cell.types %in% markers_2_plot, 'Curated.markers']), split=',')
features_2_plot <- unique(trimws(unlist(features_2_plot)))

features_2_plot <- c('AVP', 'HOPX', 'CRHBP', 'SATB1', 'LSP1', 'FLT3','MPO', 'CTSG', 'PRTN3','AZU1', 'ELANE','FCER1G', 'CSTA', 'LYZ', 'IL3RA', 'IRF7', 'IRF8','JCHAIN', 'IKZF1','DNTT', 'ADA', 'LTB','CD79A', 'EBF1', 'VPREB1','PBX1', 'PLEK', 'DAD1','CNRIP1', 'SLC40A1', 'PKIG','KLF1', 'AHSP', 'HBB','RUNX1', 'HDC', 'MS4A3')



proliferation <- cbind(setNames(as.data.frame(sc_data_subseted$S.Score), c('S Score')), setNames(as.data.frame(sc_data_subseted$G2M.Score), c('G2M Score')))
proliferation <- merge(proliferation, coords_ann, by.x=0, by.y='Row.names')



proliferation_s <- proliferation[sample(nrow(proliferation), 12000),]
proliferation_s <- proliferation_s[order(proliferation_s$`S Score`),]
proliferation_s   <- get_umap_continous(proliferation_s,   'S Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))
proliferation_g2m <- proliferation[sample(nrow(proliferation), 12000),]
proliferation_g2m <- proliferation[order(proliferation_g2m$`G2M Score`),]
proliferation_g2m <- get_umap_continous(proliferation_g2m, 'G2M Score') + theme(legend.position='bottom')+ guides(color = guide_colourbar(barwidth = 5, barheight = 0.2, ticks=FALSE))



sample_umaps <- list()
for (Sample in  unique(sc_data_subseted$Sample)){
	tmp <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', Sample, '/', Sample, '_seurat_obj_norm.rds'))
	tmp_coords <- as.data.frame(tmp@reductions$umap@cell.embeddings)
	tmp_coords <- tmp_coords[rownames(tmp_coords) %in% coords_ann$Row.names,]
	tmp_coords <- merge(tmp_coords, coords_ann[, c('Row.names', 'Cluster_names')], by.x=0, by.y='Row.names')
	sample_umaps[[Sample]] <- get_umap(tmp_coords, 'Cluster_names', Cluster_colors) + 
							  theme(legend.position= 'none', plot.title = element_text(hjust = 0.5)) + 
							  ggtitle(Sample)
}


pdf('./Plots/PaperFigures/FigS2.pdf', width=8.25, height=11.75)
cowplot::plot_grid(

	cowplot::plot_grid(
		get_umap(coords_ann, 'Cluster_names', Cluster_colors), 
		cowplot::plot_grid(plotlist=sample_umaps, nrow=2), 
		proliferation_s,
	ncol=1, rel_heights = c(1, 0.8, 0.5)
	),
	cowplot::plot_grid(get_dotplot(sc_data_subseted , features_2_plot), 
						patient_dist_cell_count + theme(legend.position = 'top'), 
						proliferation_g2m, 
	ncol=1, rel_heights = c(1, 0.8, 0.5)),
ncol=2)

dev.off()



# ===================================================
#  S3 <- PrePost 5q densities
# ===================================================

# annotation_Nerea_08 <- list()
# annotation_Nerea_08[["0"   ]] <- 'EarlyErythroid'
# annotation_Nerea_08[["1"   ]] <- 'LMPP'
# annotation_Nerea_08[["2"   ]] <- 'LateErythroid'
# annotation_Nerea_08[["3"   ]] <- 'Monocytes'
# annotation_Nerea_08[["4"   ]] <- 'LateErythroid'
# annotation_Nerea_08[["5"   ]] <- 'MEP'
# annotation_Nerea_08[["6"   ]] <- 'MEP'
# annotation_Nerea_08[["7"   ]] <- 'EarlyErythroid'
# annotation_Nerea_08[["8"   ]] <- 'Granulocite'
# annotation_Nerea_08[["9_0" ]] <- 'CLP'
# annotation_Nerea_08[["9_1" ]] <- 'CLP'
# annotation_Nerea_08[["10"  ]] <- 'HSC'
# annotation_Nerea_08[["11"  ]] <- 'LMPP'
# annotation_Nerea_08[["12"  ]] <- 'EarlyErythroid'
# annotation_Nerea_08[["13"  ]] <- 'GMP'
# annotation_Nerea_08[["14"  ]] <- 'EarlyErythroid'
# annotation_Nerea_08[["16"  ]] <- 'EarlyErythroid'
# annotation_Nerea_08[["17_1"]] <- 'EarlyErythroid'
# annotation_Nerea_08[["17_2"]] <- 'MEP'
# annotation_Nerea_08[["18"  ]] <- 'Basophil'
# annotation_Nerea_08[["19_0"]] <- '????'
# annotation_Nerea_08[["19_1"]] <- '????'
# annotation_Nerea_08[["21"  ]] <- 'Platelets'







pre_post_5q <- readRDS(paste0(getwd(), '/Data/','pre_post_5q_Annotated_final.rds'))
DefaultAssay(pre_post_5q)  <- "RNA"
Idents(pre_post_5q) <- 'integrated_snn_res.0.8_sub'


cells_selected_CASPER <- readRDS('./Data/CASPER/SelectedCells5q_PrePost_CASPER.rds')
all_5q_selected_CopyKat <- readRDS('./Data/CopyKat/all_5q_selected_cellsPrePost.rds')
real_5q_cells <- readRDS('./Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')

coords   <- as.data.frame(pre_post_5q@reductions$umap@cell.embeddings)

Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#03071e', '#006d77')
# Cluster_colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(unique(unlist(annotation_Sofia_08))))
# names(Cluster_colors) <- unique(unlist(annotation_Sofia_08))
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocite","CLP","MEP","Basophil","T","DendriticCell", '????', 'Platelets')

# clusters <- setNames(as.data.frame(pre_post_5q$integrated_snn_res.0.8_sub), c('Cluster'))
# clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Nerea_08),, drop=FALSE]
clusters <- pre_post_5q[[c('integrated_snn_res.0.8_sub', 'Cluster_names', 'Sample')]]
coords_ann <- merge(clusters, coords, by=0)
# coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
# coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Nerea_08[as.character(x[['Cluster']])][[1]] })
PrePostSample_colors <- c( '#606c38', '#dda15e')
names(PrePostSample_colors) <- unique(coords_ann$Sample)



# pre_post_real5q_cells <- pre_post_5q[pre_post_5q$Cell_id %in% real_5q_cells,]


point_legend <- cowplot::get_legend(ggplot(coords_ann[!coords_ann$Row.names %in% real_5q_cells, ], aes(x=UMAP_1, y=UMAP_2, color=Cluster_names)) + geom_point() + scale_color_manual(values=Cluster_colors) + guides(color = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=3, alpha=0.9))) + theme_classic() + 
theme(legend.position='bottom', axis.title = element_text(family = "Helvetica", size = 7), 
axis.text = element_text(family = "Helvetica", size = 5), 
axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
legend.title = element_text(size = 7, family = "Helvetica"),
legend.text = element_text(size = 5, family = "Helvetica")))

all_results <- ggplot() + 
				geom_point(coords_ann, mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 1, alpha = 1)+ 
				scale_color_manual(values=Cluster_colors)+
				geom_density_2d_filled(data= coords_ann[coords_ann$Row.names %in% real_5q_cells, ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)) , contour_var = "ndensity") +
				scale_discrete_manual("alpha",guide = "none", values =c(0, seq(0,0.8, by=0.1)))  + scale_fill_viridis_d(option = "inferno", guide='none') + theme_classic()+
				geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells, ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.2) + 
				theme(legend.position='none',  axis.title = element_text(family = "Helvetica", size = 7),  axis.text=element_blank(), axis.ticks=element_blank())


coords_ann$Sample <- factor(coords_ann$Sample , levels=c('SMD211420', 'SMD132114579'))
perSample_results <- ggplot() + 
			geom_point(coords_ann, mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 1, alpha = 1)+ 
			scale_color_manual(values=Cluster_colors) + 
			geom_density_2d_filled(data= coords_ann[coords_ann$Row.names %in% real_5q_cells, ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)),bins=20, contour_var = "ndensity") + 
			geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells, ], mapping=aes(x=UMAP_1, y=UMAP_2),bins=20, colour = "black", alpha = 0.2) + 
			scale_discrete_manual("alpha", guide = "none", values = c(0, seq(0.1, 1, by=(1-(0.1))/(20-2))))  + 
			scale_fill_viridis_d(option = "inferno") + theme_classic() + geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells, ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.3) + 
			theme(legend.position = 'none',  axis.title = element_text(family = "Helvetica", size = 7),  axis.text=element_blank(), axis.ticks=element_blank()) + 
			facet_wrap(~Sample, nrow=1)


plotter_selection <- coords_ann
plotter_selection$is_5q <- ifelse(plotter_selection$Row.names %in% real_5q_cells, '5q', 'Other')
plt_selection <- as.data.frame.matrix(table(plotter_selection$Sample, plotter_selection$is_5q))
pct_selection <- round(plt_selection/rowSums(plt_selection), 3)*100

plotter_selection$Sample <- factor(plotter_selection$Sample , levels=c('SMD211420', 'SMD132114579'))
prop_5q <- ggplot(plotter_selection, aes(x=Cluster_names, fill=is_5q)) + geom_bar(position="fill", stat="count") + facet_wrap(~Sample, nrow=1) + theme_classic()  +
		scale_fill_manual(name = 'Cell type', values=c('#586994', '#A2ABAB')) +
		theme(legend.position='bottom', axis.title = element_text(family = "Helvetica", size = 7), 
	 	axis.text = element_text(family = "Helvetica", size = 5), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 7, family = "Helvetica"),
		legend.text = element_text(size = 5, family = "Helvetica"))



plotter_casper <- coords_ann
plotter_casper$is5q <- ifelse(plotter_casper$Row.names %in% paste(cells_selected_CASPER$X1), '5q', 'Other')
plt_casper <-  as.data.frame.matrix(table(plotter_casper$Sample, plotter_casper$is5q))
pct_casper <- round(plt_casper/rowSums(plt_casper), 3)*100


plotter_CopyKat <- coords_ann
plotter_CopyKat$is5q <- ifelse(plotter_CopyKat$Row.names %in% all_5q_selected_CopyKat$Cell_id, '5q', 'Other')
plt_CopyKat <-  as.data.frame.matrix(table(plotter_CopyKat$Sample, plotter_CopyKat$is5q))
pct_CopyKat <- round(plt_CopyKat/rowSums(plt_CopyKat), 3)*100


deletion_percentages <- data.frame(Samples = c("SMD132114579","SMD211420"), 
								Karyotype = c(100,100))
deletion_percentages <- merge(deletion_percentages, setNames(pct_selection[, '5q', drop=FALSE], c('Selected cells')), by.x= 'Samples', by.y=0)
deletion_percentages <- merge(deletion_percentages, setNames(pct_casper[, '5q', drop=FALSE], c('CASPER')), by.x= 'Samples', by.y=0)
deletion_percentages <- merge(deletion_percentages, setNames(pct_CopyKat[, '5q', drop=FALSE], c('CopyKat')), by.x= 'Samples', by.y=0)

MDS_sample_colors <- setNames(c( '#606c38', '#dda15e'), unique(coords_ann$Sample))


selected_cells_Vs_Karyotype <- ggplot(reshape2::melt(deletion_percentages), aes(x = variable, y = value, group = Samples)) + 
								geom_line() + geom_point(size = 2, aes(color = Samples)) + ylim(0, 100) + geom_hline(yintercept=deletion_percentages[, 'Karyotype'], linetype="dashed", color = "grey", size=0.5, alpha =0.5) +
								scale_color_manual(values=MDS_sample_colors) +
								theme_classic() + 
								guides(color = guide_legend(nrow=2, byrow=TRUE, title.position = 'top'))+
								theme(legend.position='top', axis.title = element_text(family = "Helvetica", size = 7), 
								axis.text = element_text(family = "Helvetica", size = 5), 
								axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
								axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
								legend.title = element_text(size = 7, family = "Helvetica"),
								legend.text = element_text(size = 5, family = "Helvetica"))




pdf(paste0('./Plots/PaperFigures/FigS3.pdf'), width=8.25, height=11.75)
cowplot::plot_grid(
	cowplot::plot_grid(
		cowplot::plot_grid(
			all_results, 
			point_legend,
			nrow=2, rel_heights=c(0.9,0.1)), 
			cowplot::plot_grid(selected_cells_Vs_Karyotype,NULL,  nrow=2, rel_heights=c(0.9,0.1)),
		ncol=2, rel_widths = c(0.7, 0.3), align = 'h'),

	cowplot::plot_grid(
		perSample_results,
		prop_5q, nrow=2, rel_heights=c(0.6,0.4), align = 'v'),
	nrow=2, rel_widths=c(0.8,0.2),align = 'v')
dev.off()

