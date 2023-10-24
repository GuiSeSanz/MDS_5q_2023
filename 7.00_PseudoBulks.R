
library(ggplot2)
library(Seurat)


get_umap <- function(data, color_val, mapped_colors, legend_row=3){
	p <- ggplot(data, aes(x=UMAP_1, y=UMAP_2, color=get(color_val))) + geom_point(alpha=0.7, size = 0.5) + theme_classic() + 
		scale_color_manual(values=mapped_colors) + 
		theme(legend.position='bottom', text = element_text(family = "Helvetica", size = 7), 
		line = element_blank(),title = element_blank(), axis.text.x =element_blank(), axis.text.y=element_blank()) + 
		guides(color = guide_legend(nrow=legend_row, byrow=TRUE, override.aes = list(size=3, alpha=0.9)))
	return(p)
}


annotation_Sofia_04 <- list()
annotation_Sofia_04[['0' ]] <- 'HSC'
annotation_Sofia_04[['1' ]] <- 'HSC'
annotation_Sofia_04[['2' ]] <- 'EarlyErythroid'
annotation_Sofia_04[['3' ]] <- 'MEP'
annotation_Sofia_04[['4' ]] <- 'GMP'
annotation_Sofia_04[['5' ]] <- 'HSC'
annotation_Sofia_04[['6' ]] <- 'CLP'
annotation_Sofia_04[['7' ]] <- 'HSC'
annotation_Sofia_04[['8' ]] <- 'EarlyErythroid'
annotation_Sofia_04[['9' ]] <- 'HSC'
annotation_Sofia_04[['10']] <- 'Basophil'
annotation_Sofia_04[['11']] <- 'LMPP'
annotation_Sofia_04[['12']] <- 'DendriticCell'
annotation_Sofia_04[['13']] <- 'EarlyErythroid'
# annotation_Sofia_04[['14']] <- ''
# annotation_Sofia_04[['15']] <- ''


elder_integrated <- readRDS('./Data/all_seurat_Elder_integrated_sct.rds')
clusters <- setNames(as.data.frame(elder_integrated$integrated_snn_res.0.4), c('Cluster'))



coords   <- as.data.frame(elder_integrated@reductions$umap@cell.embeddings)

Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78')
# Cluster_colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(unique(unlist(annotation_Sofia_08))))
# names(Cluster_colors) <- unique(unlist(annotation_Sofia_08))
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocite","CLP","MEP","Basophil","T","DendriticCell")
clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_04),, drop=FALSE]
coords_ann <- merge(clusters, coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_04[as.character(x[['Cluster']])][[1]] })
ElderSample_colors <- c( '#264653', '#2a9d8f', '#e9c46a')
names(ElderSample_colors) <- unique(coords_ann$Sample)


pdf('./Plots/Integrated/Annotated_elders.pdf')
cowplot::plot_grid(
	cowplot::plot_grid(
		get_umap(coords_ann, 'Cluster_names', Cluster_colors),
		get_umap(coords_ann, 'Sample', ElderSample_colors),
		ncol=2),
	cowplot::plot_grid(
		ggplot(coords_ann, aes(x=Cluster_names, fill =Cluster_names)) + geom_bar() + 
		scale_fill_manual(values=Cluster_colors) + theme_classic() + theme( axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.position= 'none'),
		ggplot(coords_ann, aes(x=Sample, fill =Sample)) + geom_bar() + 
		scale_fill_manual(values=ElderSample_colors) + theme_classic() + theme( axis.text.x=element_text(angle=90, vjust = 0.5, hjust=1), legend.position= 'none'),
		ncol=2),
nrow=2, rel_heights=c(0.6,0.4))

dev.off()






# sc_data <- readRDS('./Data/all_seurat_integrated_sct.rds')
sc_data_MDS5q <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_integrated_sct_subcluster.rds')

annotation_Sofia_08_MDS5q <- list()
annotation_Sofia_08_MDS5q[['0' ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['1' ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['2' ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['3' ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['4' ]] <- 'pro-B'
annotation_Sofia_08_MDS5q[['5' ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['6' ]] <- 'Monocytes'
annotation_Sofia_08_MDS5q[['7' ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['8' ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['9' ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['10']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['11']] <- 'Granulocite'
annotation_Sofia_08_MDS5q[['12_0']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['12_1']] <-  'LateErythroid'
annotation_Sofia_08_MDS5q[['12_2']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['13']] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['14']] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['15']] <- 'CLP'
annotation_Sofia_08_MDS5q[['16']] <- 'MEP'
annotation_Sofia_08_MDS5q[['17_1']] <- 'GMP'
annotation_Sofia_08_MDS5q[['18']] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['19']] <- 'pro-B'
annotation_Sofia_08_MDS5q[['20']] <- 'Basophil'
annotation_Sofia_08_MDS5q[['21']] <- 'LMPP'
annotation_Sofia_08_MDS5q[['22']] <- 'T'
annotation_Sofia_08_MDS5q[['24']] <- 'Monocytes'
annotation_Sofia_08_MDS5q[['26']] <- 'DendriticCell'
annotation_Sofia_08_MDS5q[['27']] <- 'LateErythroid'


Idents(sc_data_MDS5q) <- 'integrated_snn_res.0.8_sub'
clusters <- setNames(as.data.frame(sc_data_MDS5q$integrated_snn_res.0.8_sub), c('Cluster'))
coords   <- as.data.frame(sc_data_MDS5q@reductions$umap@cell.embeddings)
Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78')
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocite","CLP","MEP","Basophil","T","DendriticCell")
MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_08_MDS5q),, drop=FALSE]
coords_ann <- merge(clusters, coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_08_MDS5q[as.character(x[['Cluster']])][[1]] })

all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') # <- 
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')
real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)

# Creating pseudo-bulk samples
# The most obvious differential analysis is to look for changes in expression between conditions.
#We perform the DE analysis separately for each label to identify cell type-specific transcriptional effects

# we need to create the aggregate of counts per celltype and per Sample

# remove the cells from the non-selected clusters 
sc_data_MDS5q$Cell_id <-  colnames(sc_data_MDS5q)
Idents(sc_data_MDS5q) <- 'Cell_id'
sc_data_MDS5q_subseted <- subset(sc_data_MDS5q, idents=coords_ann$Row.names)
tmp <- coords_ann[, 'Cluster_names', drop=FALSE]
rownames(tmp) <- coords_ann$Row.names
sc_data_MDS5q_subseted$cell_type <- tmp


norm_data_5qMDS <- as.data.frame(sc_data_MDS5q@assays$SCT@data)

security_check <- as.data.frame.matrix(table(coords_ann$Cluster_names, coords_ann$Sample))
coords_ann$is_5q <- ifelse(coords_ann$Row.names %in% real_5q_cells , 'is5q', 'Other')

	# pdf('./Plots/PaperFigures/Test.pdf')
	# ggplot(coords_ann, aes(x=Cluster_names, fill=is_5q)) + geom_bar(position='dodge', alpha = 0.8) +
	# theme_classic() + theme(axis.text.x = element_text( angle=90, vjust=0.5, hjust=1),) + facet_wrap(~Sample, scales='free') + 
	# scale_fill_manual(values=c('#0a9396', '#bb3e03')) + geom_hline(yintercept=10, linetype="dashed", color = "grey", size=0.5, alpha=0.8)
	# dev.off()



	# for (Sample in unique(sc_data_MDS5q$Sample)){
	# 	for (Cell_type in sort(unique(coords_ann$Cluster_names))){
	# 		cells_2_keep_healthy <- coords_ann[coords_ann$Cluster_names==Cell_type & coords_ann$Sample == Sample & !coords_ann$Row.names %in% real_5q_cells, ]
	# 		cells_2_keep_mds5q   <- coords_ann[coords_ann$Cluster_names==Cell_type & coords_ann$Sample == Sample &  coords_ann$Row.names %in% real_5q_cells, ]
	# 		if( any(nrow(cells_2_keep_healthy) < 10, nrow(cells_2_keep_mds5q) < 10)){
	# 			message(paste0('WARNING: Skipping ', Cell_type, ' for ', Sample, ' : not enough cells in one condition'))
	# 			next
	# 		}
	# 		if(nrow(cells_2_keep_healthy) + nrow(cells_2_keep_mds5q) != security_check[Cell_type, Sample]){
	# 			message(paste0('WARNING: dimmensions for ', Sample, 'and', Cell_type, 'dont match!'))
	# 		}
	# 		tmp_healthy <- setNames(as.data.frame(
	# 			rowSums(norm_data_5qMDS[, colnames(norm_data_5qMDS) %in% cells_2_keep_healthy$Row.names]))
	# 			, paste0(Sample, '_',Cell_type))
	# 		tmp_mds5q <- setNames(as.data.frame(
	# 			rowSums(norm_data_5qMDS[, colnames(norm_data_5qMDS) %in% cells_2_keep_mds5q$Row.names]))
	# 			, paste0(Sample, '_',Cell_type))
			
	# 		if(Sample ==  unique(sc_data_MDS5q$Sample)[1]){

	# 			tmp_healthy <- setNames(tmp_healthy, c(paste0(names(tmp_healthy), '_healthy')))
	# 			tmp_mds5q <- setNames(tmp_mds5q, c(paste0(names(tmp_mds5q), '_5q')))
	# 			all <- cbind(tmp_mds5q, tmp_healthy)
	# 		}else{

	# 			tmp_healthy <- setNames(tmp_healthy, c(paste0(names(tmp_healthy), '_healthy')))
	# 			tmp_mds5q <- setNames(tmp_mds5q, c(paste0(names(tmp_mds5q), '_5q')))
	# 			all <- cbind(all, tmp_healthy)
	# 			all <- cbind(all, tmp_mds5q)
	# 		}
	# 	}
	# }


	# design <- data.frame(Sample   = factor(stringr::str_extract(colnames(all), '^[A-Z0-9]+')),
	# 					Condition = factor(stringr::str_extract(colnames(all), '[a-z0-9]+$')), 
	# 					Cluster   = factor(stringr::str_extract(colnames(all), '(?<=_)([A-z\\-]+)(?=_)')))



	# library(edgeR)
	# library(scran)



	# dea <- pseudoBulkDGE(all, 
	# 	col.data=design, 
	# 	sample = design$Sample,
	# 	label = design$Cluster, 
	# 	condition =design$Condition, 
	# 	design = ~ 0 + Condition, 
	# 	coef=NULL,
	# 	contrast = c(1, -1))

	# # for (name in names(dea)){
	# # 	test <- dea[[name]]
	# # 	test <- test[order(test$PValue),]
	# # 	print(head(test))
	# # }

	# is.de <- decideTestsPerLabel(dea, threshold=0.05)
	# summarizeTestsPerLabel(is.de)


	# test_de_ann <- coords_ann[coords_ann$Cluster_names == name, ]
	# test_de <- norm_data_5qMDS[rownames(test)[1],test_de_ann$Row.names]
	# test_de <- reshape2::melt(test_de)
	# test_de <- merge(test_de, test_de_ann , by.x='variable', by.y='Row.names')

	# pdf('./Plots/Integrated/Test_pseudo.pdf')
	# ggplot(test_de ,aes(x=is_5q, y=value)) + geom_boxplot()
	# edgeR::plotBCV(metadata(test)$y)
	# dev.off()




	# de.results <- pseudoBulkDGE(all, 
	# 	col.data=design, 
	# 	sample = design$Sample,
	# 	label = design$Cluster, 
	# 	condition =design$Condition, 
	# 	design = ~ 0 + Condition, 
	# 	contrast="Condition5q - Conditionhealthy" 
	# )


	# test_de_ann <- coords_ann[coords_ann$Cluster_names == 'MEP', ]
	# test_de <- norm_data_5qMDS[rownames(de.results$MEP[order(de.results$MEP$PValue),])[1],test_de_ann$Row.names]
	# test_de <- reshape2::melt(test_de)
	# test_de <- merge(test_de, test_de_ann , by.x='variable', by.y='Row.names')

	# pdf('./Plots/Integrated/Test_pseudo.pdf')
	# ggplot(test_de ,aes(x=is_5q, y=value)) + geom_boxplot()
	# edgeR::plotBCV(metadata(test)$y)
	# dev.off()


	# is.de <- decideTestsPerLabel(de.results, threshold=0.05)
	# summarizeTestsPerLabel(is.de)


	# pbmc.sce <- as.SingleCellExperiment(sc_data_MDS5q)

	# pbmc.sce <- pbmc.sce[, colnames(pbmc.sce) %in% coords_ann$Row.names]
	# pbmc.sce$Cell_types <- coords_ann$Cluster_names





###########
# LIBRA 
###########

norm_data_5qMDS <- as.data.frame(sc_data_MDS5q_subseted@assays$SCT@data)

sc_data_MDS5q_subseted$replicate <- sc_data_MDS5q_subseted$Sample
sc_data_MDS5q_subseted$label <- ifelse( colnames(sc_data_MDS5q_subseted) %in% real_5q_cells, '5q', 'other')
coords_ann$is_5q  <- ifelse( coords_ann$Row.names %in% real_5q_cells, '5q', 'other') 

#
	# DefaultAssay(object = sc_data_MDS5q_subseted) <- "SCT"

	# DE = Libra::run_de(sc_data_MDS5q_subseted)
	# DE_df <- as.data.frame(DE)

	# DE_df <- DE_df[DE_df$p_val_adj < 0.05,]


	# model.matrix(~ group, data = targets)


	# plot_genes_from_DE <- function(de_results, annotation = coords_ann, norm_data =norm_data_5qMDS){
	# 	plot_list <- list()
	# 	for (CT in unique(de_results$cell_type)){
	# 		genes_2_check <- de_results[de_results$cell_type == CT, 'gene']
	# 		tmp_CT <- norm_data_5qMDS[rownames(norm_data_5qMDS) %in% genes_2_check, colnames(norm_data_5qMDS) %in% annotation[annotation$Cluster_names == CT, 'Row.names']]
	# 		for (gene in genes_2_check){
	# 			tmp <- reshape2::melt( tmp_CT[gene, ] )
	# 			tmp$is_5q <- ifelse(tmp$variable %in% annotation[annotation$is_5q=='5q', 'Row.names'], '5q', 'No 5q')
	# 			plot <- ggplot(tmp, aes(x=is_5q, y=value, fill=is_5q)) + geom_violin(alpha=0.8)+ theme_classic() + 
	# 					ggtitle(paste0('Gene: ', gene, '  adj.pVal: ', round(de_results[de_results$gene==gene,'p_val_adj'], 3), 
	# 					'\nLogFC: ', round(de_results[de_results$gene==gene,'avg_logFC'], 2)))
	# 			plot_list[[paste0(CT, '_', gene)]] <- plot
	# 		}
	# 	}
	# 	return(plot_list)
	# }

	# plots <- plot_genes_from_DE(DE_df)

	# pdf('./Plots/Integrated/Test_pseudo_res.pdf')
	# for (i in seq(1, 17, by=4)){
	# 	print(cowplot::plot_grid(plots[[i]],
	# 	plots[[i+1]],
	# 	plots[[i+2]],
	# 	plots[[i+3]], ncol=2))
	# }
	# dev.off()

	# Idents(sc_data_MDS5q_subseted) <- 'cell_type'
	# VlnPlot(sc_data_MDS5q_subseted, features='TREML1', split.by='label')



	# Idents(sc_data_MDS5q_subseted) <- 'cell_type'

	# MEP_5q    <- coords_ann[coords_ann$Cluster_names == 'MEP' & coords_ann$is_5q == 'other', 'Row.names']
	# MEP_non5q <- coords_ann[coords_ann$Cluster_names == 'MEP' & coords_ann$is_5q == '5q', 'Row.names']
	# aa <- Seurat::FindMarkers(sc_data_MDS5q_subseted, slot='counts',test.use= 'DESeq2', cells.1 = MEP_5q, cells.2 = MEP_non5q )


	# Idents(sc_data_MDS5q_subseted) <- 'label'
	# aa <- Seurat::FindMarkers(sc_data_MDS5q_subseted, slot='counts',test.use= 'DESeq2', ident.1 = '5q', ident.2 = 'other', group.by='label',)





#

sc_data_MDS5q_subseted <- readRDS('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')

DefaultAssay(sc_data_MDS5q_subseted)  <- "SCT"

all_results <- list()
for( Sample in unique(sc_data_MDS5q_subseted$Sample)){
	Idents(sc_data_MDS5q_subseted) <- 'Sample'
	print(Sample)
	tmp_smp <- subset(sc_data_MDS5q_subseted, idents =Sample)
	for ( CT in unique(sc_data_MDS5q_subseted$Cluster_names) ){
		print(paste0(Sample, ' -- ', CT ))
		Idents(tmp_smp) <- 'Cluster_names'
		tmp <- subset(tmp_smp, idents = CT)
		Idents(tmp) <- "cell_5q"
		if(all(table(tmp$cell_5q)[['del5q']] > 10 & table(tmp$cell_5q)[['normal']] >10)){
			results <- FindMarkers(tmp, slot='counts', test.use = 'DESeq2', ident.1 = "del5q", ident.2 = "normal", verbose = FALSE)
			all_results[[Sample]][[CT]] <- results
		}
	}
}


saveRDS(all_results, './Data/Results_DE_per_CET_5q_healthy_per_patient.rds')

library(fgsea)
all_results <- readRDS('./Data/Results_DE_per_CET_5q_healthy_per_patient.rds')

pathways_c2 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c2.cp.v7.5.1.symbols.gmt')
pathways_c1 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c1.all.v7.5.1.symbols.gmt')

get_tornado <- function(results){
	p <- ggplot(results, aes(x = NES, y = reorder(pathway, NES), fill = is_pos)) + 
			geom_bar(stat = "identity", position = "identity") + theme_classic() +
			ylab('Pathways in C2')+
			theme(legend.position="none", text = element_text(family = "Helvetica", size = 7)) + 
			scale_fill_manual(values=c('#006e90', '#cc2936'))
	return(p)
}

get_volcano <- function(results, SampleName, CType){
	p <-  EnhancedVolcano::EnhancedVolcano(tmp,
				lab = rownames(results),
				x = 'avg_log2FC', captionLabSize=7,
				y = 'p_val_adj',col=c('black', 'black', 'black', 'red3'),
				subtitle = "5q Vs healthy", axisLabSize = 5, title=paste0(SampleName, ' -- ', CType),
				pCutoff = 5e-3, labSize=2,
				FCcutoff = 0.2) +
				theme(legend.position='none', text = element_text(family = "Helvetica", size = 2),
				axis.title= element_text(family = "Helvetica", size = 7))
	return(p)
}

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
	p <- p+ scale_y_continuous(limits = c(-1, 1.5))
	}
	return(p)
}

get_filtered_values <- function(data, pval=0.05, logfc=0.2){
	data <- data[rowSums(is.na(data))==0,]
	data <- data[data$p_val_adj <pval & abs(data$avg_log2FC) >logfc,]
	return(data)
}

plots <- list()
for( Sample in names(all_results)){
	tmp_smp <- all_results[[Sample]]
	for( CT in names(tmp_smp)){
		tmp <- tmp_smp[[CT]]
		# remove rows with NAs
		tmp <- tmp[rowSums(is.na(tmp))==0,]
		# print(summary(tmp$p_val_adj))
		tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_log2FC) > 0.2,]
		if(nrow(tmp)!=0){
			tmp <- tmp[order(tmp$avg_log2FC, decreasing=TRUE),]
			ranks <- setNames(tmp$avg_log2FC, rownames(tmp))
			results_c1 <- fgsea::fgsea(pathways_c1, ranks, minSize=15, maxSize = 500)
			results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 500)
			volcano <- get_volcano(tmp, Sample, CT)
			results_c1 <- results_c1[results_c1$padj < 0.05,]
			results_c1 <- results_c1[order(results_c1$NES, decreasing=TRUE),]
			results_c1$is_pos <- ifelse(results_c1$NES<=0, 'Neg', 'Pos')
			results_c1$pathway <- factor(results_c1$pathway)
			results_c2 <- results_c2[results_c2$padj < 0.05,]
			results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
			results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
			results_c2$pathway <- factor(results_c2$pathway)
			if(nrow(results_c1) != 0){
				print(paste0('C1:  ',Sample, ' -- ', CT))
			}
			if(nrow(results_c2) != 0){
				print(paste0('C2:  ',Sample, ' -- ', CT))
				plots[[paste0(Sample, ' -- ', CT)]]<- cowplot::plot_grid(
																	volcano,
																	get_tornado(results_c2),
																	ncol=1)
				}else{
					plots[[paste0(Sample, ' -- ', CT)]] <- volcano
				}
		}
	}
}

ven_list <- list()
for( CT in names(all_results[[1]])){
	tmp <- rownames(get_filtered_values(all_results[['SMD34459']][[CT]], 0.05, 0.2))
	results_venn <- list(
		SMD34459 = rownames(get_filtered_values(all_results[['SMD34459']][[CT]], 0.05, 0.2)),
		SMD35109 = rownames(get_filtered_values(all_results[['SMD35109']][[CT]], 0.05, 0.2)),
		SMD37209 = rownames(get_filtered_values(all_results[['SMD37209']][[CT]], 0.05, 0.2)),
		SMD35303 = rownames(get_filtered_values(all_results[['SMD35303']][[CT]], 0.05, 0.2))
	)
	ven_list[[CT]] <- get_gvenn(results_venn)

}




pdf('./Plots/Integrated/Test_tornado.pdf')
for( pl in plots){
	print(pl)
}
dev.off()


for(Sample in names(all_results)){
	tmp <- all_results[[Sample]]
	tmp <- lapply(tmp, get_filtered_values)
	WriteXLS::WriteXLS(tmp, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/',Sample, 'SingleCell_DE.xlsx'), SheetNames = names(tmp),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

}





################################################################################
# Compare results between pseudoBulk and single cell DE 
################################################################################
# data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') # <- 
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')

real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
cells_to_discard <- c(setdiff(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id), 
					  setdiff(cells_selected_CASPER$Cell_id, all_5q_selected_cells$Cell_id))

# NONE OF THE CELLS TO DISCARD ARE IN THE SEURAT OBJ.

expression <- as.data.frame(data_5q@assays$RNA@counts)
# rownames(expression) <- gsub('-', '\\.', rownames(expression))



# merged_seurat <- merge(elder_integrated_subseted, sc_data_MDS5q_subseted)

metadata <- setNames(data_5q[[c('cell_5q', 'Sample', 'Cluster_names')]], 
			c('label', 'replicate', 'cell_type'))

library(Libra)

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_flt_01 <- DE[DE$p_val_adj < 0.1,]
DE_flt_005 <- DE[DE$p_val_adj < 0.05,]

pdf('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_PseudoBulk_qvalHist.pdf')
hist(DE$p_val_adj, breaks=100)
dev.off()

Idents(data_5q) <- 'Sample'
sample_objs <- SplitObject(data_5q, split.by = "Sample")


de_resulst_Patient_DESeq2 <- list()
de_resulst_Patient_MAST <- list()
for (smp in names(sample_objs)){
	print(smp)
	tmp <- sample_objs[[smp]]
	DefaultAssay(tmp) <- 'SCT'
	Idents(tmp) <- 'cell_5q'
	de_resulst_Patient_DESeq2[[smp]] <- Seurat::FindMarkers(tmp, slot='counts',test.use= 'DESeq2', ident.1 = 'del5q', ident.2 = 'normal' )
	# de_resulst_Patient_MAST[[smp]]   <- Seurat::FindMarkers(tmp, assay='SCT', test.use= 'MAST', ident.1 = 'del5q', ident.2 = 'normal' )
}

saveRDS(de_resulst_Patient_DESeq2, '/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientNoCellType.rds')
# saveRDS(de_resulst_Patient_MAST, '/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_MAST_PatientByPatientNoCellType.rds')

results_flt_DESeq2 <- lapply(de_resulst_Patient_DESeq2,  function(x) na.omit(x[x$p_val_adj < 0.05,]))
WriteXLS::WriteXLS(results_flt_DESeq2, 
ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientNoCellType_Filtered.xlsx'),
SheetNames = names(results_flt_DESeq2),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)


# results_flt_MAST <- lapply(de_resulst_Patient_MAST,  function(x) na.omit(x[x$p_val_adj < 0.05,]))
# WriteXLS::WriteXLS(results_flt_MAST, 
# ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_MAST_PatientByPatientNoCellType_Filtered.xlsx'),
# SheetNames = names(results_flt_MAST),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)



# Also per cell type
Idents(data_5q) <- 'Sample'
sample_objs <- SplitObject(data_5q, split.by = "Sample")

de_resulst_Patient_DESeq2 <- list()
de_resulst_Patient_MAST <- list()
for (smp in names(sample_objs))
{
	print(smp)
	for (cluster in sort(unique(data_5q$Cluster_names)))
	{	
		print(cluster)
		DefaultAssay(tmp_smp) <- 'SCT'
		Idents(tmp_smp) <- 'Cluster_names'
		tmp <- subset(tmp_smp, idents=cluster)
		Idents(tmp) <- 'cell_5q'
		if(nrow(as.data.frame(table(tmp$cell_5q)))>1 & all(as.data.frame(table(tmp$cell_5q))$Freq>3)){
			de_resulst_Patient_DESeq2[[paste0(smp, '_', cluster)]] <- Seurat::FindMarkers(tmp, slot='counts',test.use= 'DESeq2', ident.1 = 'del5q', ident.2 = 'normal' )
		}else{
			print('Note enough cells in one condition')
		}
		# de_resulst_Patient_MAST[[paste0(smp, '_', cluster)]]   <- Seurat::FindMarkers(tmp, assay='SCT',  test.use= 'MAST',   ident.1 = 'del5q', ident.2 = 'normal' )
	}
}

saveRDS(de_resulst_Patient_DESeq2, '/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType.rds')
# saveRDS(de_resulst_Patient_MAST, '/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_MAST_PatientByPatientByCellType.rds')

results_flt_DESeq2 <- lapply(de_resulst_Patient_DESeq2,  function(x) na.omit(x[x$p_val_adj < 0.05,]))
WriteXLS::WriteXLS(results_flt_DESeq2, 
ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType_Filtered.xlsx'),
SheetNames = names(results_flt_DESeq2),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)


# results_flt_MAST <- lapply(de_resulst_Patient_MAST,  function(x) na.omit(x[x$p_val_adj < 0.05,]))
# WriteXLS::WriteXLS(results_flt_MAST, 
# ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_MAST_PatientByPatientByCellType_Filtered.xlsx'),
# SheetNames = names(results_flt_MAST),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

de_resulst_Patient_DESeq2 <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType.rds')
results_flt_DESeq2 <- lapply(de_resulst_Patient_DESeq2,  function(x) na.omit(x[x$p_val_adj < 0.05,]))

pdf('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType_Filtered_VENN.pdf')
for (cluster in sort(unique(data_5q$Cluster_names)))
{	
	list_2_plot <- list()
	for (smp in names(sample_objs))
	{	
		experiment_name <- paste0(smp, '_', cluster)
		tmp <- rownames(results_flt_DESeq2[[experiment_name]])
		list_2_plot[[smp]] <- tmp
	}
	print(ggvenn::ggvenn(list_2_plot,
		fill_color = destiny::cube_helix(length(names(sample_objs))),
		stroke_size = 0.4,
		show_percentage = TRUE,
		fill_alpha = 0.4,
		stroke_color = 'white',
		stroke_alpha = 1,
		stroke_linetype = 'solid',
		text_color = 'black',
		set_name_size = 4,
		text_size = 2)+ggtitle(cluster))
	# print(UpSetR::upset(UpSetR::fromList(list_2_plot), order.by = "freq"))
}
dev.off()

# check overlapping pathways
data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
de_resulst_Patient_DESeq2 <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType.rds')
results_flt_DESeq2 <- lapply(de_resulst_Patient_DESeq2,  function(x) na.omit(x[x$p_val_adj < 0.05,]))
pathways_c2 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c2.cp.v7.5.1.symbols.gmt')

get_venn <- function(list_2_plot, title)
{
	venn <- ggvenn::ggvenn(list_2_plot,
		fill_color = destiny::cube_helix(length(unique(data_5q$Sample))),
		stroke_size = 0.4,
		show_percentage = TRUE,
		fill_alpha = 0.4,
		stroke_color = 'white',
		stroke_alpha = 1,
		stroke_linetype = 'solid',
		text_color = 'black',
		set_name_size = 4,
		text_size = 2)+ggtitle(title)
	return(venn)
}

pdf('./Plots/DESeq2_pathway_overlap_DEG.pdf')
for (cluster in sort(unique(data_5q$Cluster_names)))
{	
	plotter_BP <- list()
	plotter_MF <- list()
	plotter_CC <- list()
	plotter_KG <- list()
	for (smp in unique(stringr::str_extract(names(de_resulst_Patient_DESeq2), '^[A-Z0-9]+')))
	{	
		print(paste0(cluster, '____', smp))
		experiment_name <- paste0(smp, '_', cluster)
		tmp <- rownames(results_flt_DESeq2[[experiment_name]])
		if(length(tmp)>1)
		{
			names = clusterProfiler::bitr(sub('\\.','-', tmp ), fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
			results_BP = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_MF = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_CC = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_KG = as.data.frame(clusterProfiler::enrichKEGG(names[, 'ENTREZID'], organism = "hsa",pAdjustMethod = "BH",pvalueCutoff  = 0.05))
			plotter_BP[[smp]] <-  results_BP[results_BP$p.adjust < 0.05, 'Description']
			plotter_MF[[smp]] <-  results_MF[results_MF$p.adjust < 0.05, 'Description']
			plotter_CC[[smp]] <-  results_CC[results_CC$p.adjust < 0.05, 'Description']
			plotter_KG[[smp]] <-  results_KG[results_KG$p.adjust < 0.05, 'Description']
		}else{
			plotter_BP[[smp]] <- runif(1)
			plotter_MF[[smp]] <- runif(1)
			plotter_CC[[smp]] <- runif(1)
			plotter_KG[[smp]] <- runif(1)
		}
	}
	print(cowplot::plot_grid(
		get_venn(plotter_BP, paste0(cluster, ':: Biological Process')),
		get_venn(plotter_MF, 'Molecular Function'),
		get_venn(plotter_CC, 'Cellular Component'),
		get_venn(plotter_KG, 'KEGG'),
		ncol = 2, nrow = 2)
	)
}
dev.off()



get_tornado <- function(results, title){
	p <- ggplot(results, aes(x = NES, y = reorder(pathway, NES), fill = is_pos)) + 
			geom_bar(stat = "identity", position = "identity") + theme_classic() +
			ylab('Pathways in C2')+
			theme(legend.position="none", text = element_text(family = "Helvetica", size = 7), axis.text.y= element_text(size = 5),) + 
			scale_fill_manual(values=c(Neg='#006e90', Pos='#cc2936')) + ggtitle(title)
	return(p)
}

pdf('./Plots/DESeq2_Tornado.pdf')
for (experiment_name in names(results_flt_DESeq2))
{	
	tmp <- results_flt_DESeq2[[experiment_name]]
	tmp <- tmp[order(tmp$avg_log2FC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_log2FC, rownames(tmp))
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
		message('printing')
		print(get_tornado(results_c2, experiment_name))
	}
}
dev.off()


# check if patients are grouped the way I want :D
de_resulst_Patient_DESeq2 <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/DE_results/5qVsNon5q_SingleCell_DESeq2_PatientByPatientByCellType.rds')
Idents(data_5q) <- 'Cluster_names'
DefaultAssay(data_5q)<- 'SCT'

plotter_sum        <- data.frame(Sample= unique(data_5q$Sample))
plotter_mean       <- data.frame(Sample= unique(data_5q$Sample))
plotter_sum_split  <- data.frame(Sample= c("SMD34459_del5q", "SMD34459_normal", "SMD35109_del5q", "SMD35109_normal", "SMD35303_del5q" , "SMD35303_normal", "SMD37209_del5q", "SMD37209_normal"))
plotter_mean_split <- data.frame(Sample= c("SMD34459_del5q", "SMD34459_normal", "SMD35109_del5q", "SMD35109_normal", "SMD35303_del5q" , "SMD35303_normal", "SMD37209_del5q", "SMD37209_normal"))


for(cell_type in unique(data_5q$Cluster_names))
{
	print(cell_type)
	tmp <- subset(data_5q, ident=cell_type)
	results <- de_resulst_Patient_DESeq2[grepl(cell_type, names(de_resulst_Patient_DESeq2))]
	results <- lapply(results, function(x) na.omit(x[x$p_val_adj < 0.05,]))
	results <- list(de_score = c(unname(unlist(lapply(results, function(x) rownames(x))))))
	nbin = switch(  
		cell_type,  
		"LateErythroid" = 24,
		"EarlyErythroid" = 24,
		"LMPP" = 24,
		"MK_Prog" = 20,
		"MEP" = 24,
		"GMP" = 24,
		"CLP" = 24,
		"HSC" = 24,
		"Basophil" = 15,
		"T" = 15,
		"Granulocyte" = 24,
		"Monocytes" = 24,
		"DendriticCell" = 15,
		"pro-B" = 24
	)  
	tmp <- AddModuleScore(
		object = tmp,
		features = results,
		ctrl = 5,
		name = 'de_score', 
		nbin=nbin)
	data <-  tmp[[c('Sample', 'de_score1')]]
	data_split <-  tmp[[c('Sample', 'cell_5q', 'de_score1')]]
	data_split$pheno <- paste0(data_split$Sample, '_', data_split$cell_5q)

	data_split_sum  <- setNames(aggregate(data_split$de_score1, by=list(Sample=data_split$pheno), FUN=sum),  c('Sample', cell_type))
	data_split_mean <- setNames(aggregate(data_split$de_score1, by=list(Sample=data_split$pheno), FUN=mean), c('Sample', cell_type))
	data_sum        <- setNames(aggregate(data$de_score1, by=list(Sample=data$Sample), FUN=sum),  c('Sample', cell_type))
	data_mean       <- setNames(aggregate(data$de_score1, by=list(Sample=data$Sample), FUN=mean), c('Sample', cell_type))

	plotter_mean_split <- merge(plotter_mean_split, data_split_mean, by='Sample', all.x=TRUE)
	plotter_sum_split  <- merge(plotter_sum_split,  data_split_sum,  by='Sample', all.x=TRUE)
	plotter_sum        <- merge(plotter_sum, data_sum,   by='Sample', all.x=TRUE)
	plotter_mean       <- merge(plotter_mean, data_mean, by='Sample', all.x=TRUE)
}

rownames(plotter_sum) <- plotter_sum$Sample
rownames(plotter_mean) <- plotter_mean$Sample
rownames(plotter_mean_split) <- plotter_mean_split$Sample
rownames(plotter_sum_split) <- plotter_sum_split$Sample

plotter_sum <- plotter_sum[,-1]
plotter_mean <- plotter_mean[,-1]
plotter_mean_split <- plotter_mean_split[,-1]
plotter_sum_split <- plotter_sum_split[,-1]


pdf('./Plots/Clustering_Samples_DESeq2.pdf')
ann_colors <- c('#2a9d8f', '#f4a261')
names(ann_colors) <- c('RollerCoaster', 'Flat')
ann_colors <- list(Pattern=ann_colors)

annotation_row <- data.frame(Sample = rownames(plotter_sum))
annotation_row$Pattern <- ifelse(grepl('SMD34459|SMD37209', annotation_row$Sample), 'RollerCoaster', 'Flat')
annotation_row_split <- data.frame(Sample = rownames(plotter_sum_split))
annotation_row_split$Pattern <- ifelse(grepl('SMD34459|SMD37209', annotation_row_split$Sample), 'RollerCoaster', 'Flat')
rownames(annotation_row) <- annotation_row$Sample
rownames(annotation_row_split) <- annotation_row_split$Sample

annotation_row <- annotation_row[,-1, drop=FALSE]
annotation_row_split <- annotation_row_split[,-1, drop=FALSE]
colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100)
pheatmap::pheatmap(plotter_sum,  color=colors, annotation_row=annotation_row, annotation_colors=ann_colors, main='Sum')
pheatmap::pheatmap(plotter_mean, color=colors, annotation_row=annotation_row, annotation_colors=ann_colors, main='Mean')
pheatmap::pheatmap(plotter_sum_split,  color=colors, annotation_row=annotation_row_split, annotation_colors=ann_colors, main='Sum')
pheatmap::pheatmap(plotter_mean_split, color=colors, annotation_row=annotation_row_split, annotation_colors=ann_colors, main='Mean')
dev.off()


get_overlap_results <- function(result_list)
{	
	if(length(result_list)==0){
		return(NULL)
	}
	rst <- data.frame(SAMPLE = NULL, PATH = NULL)
	for( smp in names(result_list)){
		tmp <- result_list[[smp]]
		if(length(tmp)==0){
			tmp <- runif(1)
		}
		rst <- rbind(rst, data.frame(SAMPLE = smp, PATH = tmp))
	}
	results <- data.frame(SAMPLES=NULL, PATH=NULL)
	for (path in unique(rst$PATH)){
		tmp <- rst[rst$PATH==path,]
		aa <- paste(tmp$SAMPLE, collapse='|')
		results <- rbind(results, data.frame(SAMPLES=aa, PATH=path))
	}
	return(results)
}

get_venn <- function(list_2_plot, title)
{
	venn <- ggvenn::ggvenn(list_2_plot,
		fill_color = destiny::cube_helix(length(unique(data_5q$Sample))),
		stroke_size = 0.4,
		show_percentage = TRUE,
		fill_alpha = 0.4,
		stroke_color = 'white',
		stroke_alpha = 1,
		stroke_linetype = 'solid',
		text_color = 'black',
		set_name_size = 4,
		text_size = 2)+ggtitle(title)
	return(venn)
}

pdf('./Plots/DESeq2_pathway_overlap_DEG_POS.pdf')
common_plotter_BP <- list()
common_plotter_MF <- list()
common_plotter_CC <- list()
common_plotter_KG <- list()
for (cluster in sort(unique(data_5q$Cluster_names)))
{	
	plotter_BP <- list()
	plotter_MF <- list()
	plotter_CC <- list()
	plotter_KG <- list()
	for (smp in unique(stringr::str_extract(names(de_resulst_Patient_DESeq2), '^[A-Z0-9]+')))
	{	
		print(paste0(cluster, '____', smp))
		experiment_name <- paste0(smp, '_', cluster)
		tmp <- results_flt_DESeq2[[experiment_name]]
		tmp <- rownames(tmp[tmp$avg_log2FC>0,])
		if(length(tmp)>1)
		{
			names = clusterProfiler::bitr(sub('\\.','-', tmp ), fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
			results_BP = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_MF = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_CC = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_KG = as.data.frame(clusterProfiler::enrichKEGG(names[, 'ENTREZID'], organism = "hsa",pAdjustMethod = "BH",pvalueCutoff  = 0.05))
			plotter_BP[[smp]] <-  results_BP[results_BP$p.adjust < 0.05, 'Description']
			plotter_MF[[smp]] <-  results_MF[results_MF$p.adjust < 0.05, 'Description']
			plotter_CC[[smp]] <-  results_CC[results_CC$p.adjust < 0.05, 'Description']
			plotter_KG[[smp]] <-  results_KG[results_KG$p.adjust < 0.05, 'Description']
		}else{
			plotter_BP[[smp]] <- runif(1)
			plotter_MF[[smp]] <- runif(1)
			plotter_CC[[smp]] <- runif(1)
			plotter_KG[[smp]] <- runif(1)
		}
		common_plotter_MF[[cluster]] <- get_overlap_results(plotter_MF)
		common_plotter_BP[[cluster]] <- get_overlap_results(plotter_BP)
		common_plotter_KG[[cluster]] <- get_overlap_results(plotter_KG)
		common_plotter_CC[[cluster]] <- get_overlap_results(plotter_CC)
	}
	print(cowplot::plot_grid(
		get_venn(plotter_BP, paste0(cluster, ':: Biological Process')),
		get_venn(plotter_MF, 'Molecular Function'),
		get_venn(plotter_CC, 'Cellular Component'),
		get_venn(plotter_KG, 'KEGG'),
		ncol = 2, nrow = 2)
	)
}
dev.off()

WriteXLS::WriteXLS(common_plotter_MF, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_POS_pathways_MF.xlsx'), SheetNames = names(common_plotter_MF),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_BP, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_POS_pathways_BP.xlsx'), SheetNames = names(common_plotter_BP),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_KG, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_POS_pathways_KG.xlsx'), SheetNames = names(common_plotter_KG),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_CC, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_POS_pathways_CC.xlsx'), SheetNames = names(common_plotter_CC),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)






pdf('./Plots/DESeq2_pathway_overlap_DEG_NEG.pdf')
common_plotter_BP <- list()
common_plotter_MF <- list()
common_plotter_CC <- list()
common_plotter_KG <- list()
for (cluster in sort(unique(data_5q$Cluster_names)))
{	
	plotter_BP <- list()
	plotter_MF <- list()
	plotter_CC <- list()
	plotter_KG <- list()
	for (smp in unique(stringr::str_extract(names(de_resulst_Patient_DESeq2), '^[A-Z0-9]+')))
	{	
		print(paste0(cluster, '____', smp))
		experiment_name <- paste0(smp, '_', cluster)
		tmp <- results_flt_DESeq2[[experiment_name]]
		tmp <- rownames(tmp[tmp$avg_log2FC<0,])
		if(length(tmp)>1)
		{
			names = clusterProfiler::bitr(sub('\\.','-', tmp ), fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
			results_BP = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_MF = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_CC = as.data.frame(clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE))
			results_KG = as.data.frame(clusterProfiler::enrichKEGG(names[, 'ENTREZID'], organism = "hsa",pAdjustMethod = "BH",pvalueCutoff  = 0.05))
			plotter_BP[[smp]] <-  results_BP[results_BP$p.adjust < 0.05, 'Description']
			plotter_MF[[smp]] <-  results_MF[results_MF$p.adjust < 0.05, 'Description']
			plotter_CC[[smp]] <-  results_CC[results_CC$p.adjust < 0.05, 'Description']
			plotter_KG[[smp]] <-  results_KG[results_KG$p.adjust < 0.05, 'Description']
		}else{
			plotter_BP[[smp]] <- runif(1)
			plotter_MF[[smp]] <- runif(1)
			plotter_CC[[smp]] <- runif(1)
			plotter_KG[[smp]] <- runif(1)
		}
	common_plotter_MF[[cluster]] <- get_overlap_results(plotter_MF)
	common_plotter_BP[[cluster]] <- get_overlap_results(plotter_BP)
	common_plotter_KG[[cluster]] <- get_overlap_results(plotter_KG)
	common_plotter_CC[[cluster]] <- get_overlap_results(plotter_CC)
	}
	print(cowplot::plot_grid(
		get_venn(plotter_BP, paste0(cluster, ':: Biological Process')),
		get_venn(plotter_MF, 'Molecular Function'),
		get_venn(plotter_CC, 'Cellular Component'),
		get_venn(plotter_KG, 'KEGG'),
		ncol = 2, nrow = 2)
	)
}

dev.off()

WriteXLS::WriteXLS(common_plotter_MF, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_NEG_pathways_MF.xlsx'), SheetNames = names(common_plotter_MF),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_BP, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_NEG_pathways_BP.xlsx'), SheetNames = names(common_plotter_BP),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_KG, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_NEG_pathways_KG.xlsx'), SheetNames = names(common_plotter_KG),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(common_plotter_CC, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/common_NEG_pathways_CC.xlsx'), SheetNames = names(common_plotter_CC),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

