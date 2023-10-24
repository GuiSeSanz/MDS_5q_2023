
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

get_prop_tables_per_sample <- function(Sobj){
	print(
		prop.table(table(Sobj$Sample, Sobj$cell_5q), margin=1)*100
	)
}


post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
all_5q_depleted_cells_COPYKAT <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')
#  There are less cells as some clusters have disapeared based on the annotation

# elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')



PrePost_CASPER <- readRDS('./Data/CASPER/SelectedCells5q_PrePost_CASPER.rds')
PrePost_Copykat <- readRDS( './Data/CopyKat/all_5q_selected_cellsPrePost.rds')

## 
	# get_prop_tables_per_sample(post_data)
	# get_prop_tables_per_sample(MDS_data)
	# get_prop_tables_per_sample(pre_post_data)


	# coords_ann <- rbind(FetchData(post_data,     vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')), 
	# 					FetchData(pre_post_data, vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')))

	# coords_ann <- coords_ann[coords_ann$Sample != 'SMD211420',]

	# Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722')
	# names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog')


	# Enrichment_results <- data.frame(Sample=NULL, Cluster=NULL, pval=NULL)
	# for (Sample in sort(unique(coords_ann$Sample))){
	# 	tmp_Smp <- coords_ann[coords_ann$Sample == Sample,]
	# 	all_cells <- nrow(tmp_Smp)
	# 	all_cells_5q <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'del5q',])
	# 	all_cells_normal <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'normal',])

	# 	for (CT in sort(unique(coords_ann$Cluster_names))){
	# 		tmp_CT <- tmp_Smp[tmp_Smp$Cluster_names ==CT,]
	# 		# phyper(q, m, n, k)
	# 		pval_phyper <- phyper(sum(tmp_CT$cell_5q == 'del5q') -1 , all_cells_5q, all_cells_normal, nrow(tmp_CT), lower.tail=FALSE)
	# 		Enrichment_results <- rbind(Enrichment_results, data.frame(Sample=Sample, Cluster=CT, pval=pval_phyper))
	# 	}
	# }
	# Enrichment_results_df <- reshape2::dcast(Enrichment_results, Sample ~ Cluster, value.var = 'pval')
	# rownames(Enrichment_results_df) <-  Enrichment_results_df$Sample
	# Enrichment_results_df$Sample <- NULL







	# formatter <- function(...){
	#   function(x) format(round(x, 1), ...)
	# }

	# plot_list_ratios <- list()
	# for (Sample in sort(unique(coords_ann$Sample))){
	# 	tmp <- coords_ann[coords_ann$Sample == Sample,]
	# 	tmp <- table(tmp$Cluster_names, tmp$cell_5q)
	# 	ratio_all <- sum(tmp[,1]) / sum(tmp[, 2])
	# 	tmp <- setNames(as.data.frame(tmp[,1] / tmp[, 2]), c('Ratio'))
	# 	tmp$Cluster_names <- factor(rownames(tmp))
	# 	tmp$Sample <- Sample
	# 	labels_bold <- Enrichment_results[Enrichment_results$Sample == Sample & Enrichment_results$pval < 0.05, 'Cluster']
	# 	labels_bold <- ifelse(levels(tmp$Cluster_names) %in% labels_bold, 'bold', 'plain')
	# 	p <- ggplot(tmp, aes(x=Cluster_names, y=Ratio, group = Sample)) + 
	# 	geom_hline(yintercept = ratio_all, color='grey', linetype='dotted') +
	# 	geom_point(color='dodgerblue4') + geom_line(color='dodgerblue4') +
	# 	theme_classic() + labs(x = "Cell types", y = "Ratio 5qCells/normalcells") +
	# 	theme(legend.position='top', 
	# 	axis.text = element_text(family = "Helvetica", size = 7), 
	# 	axis.title = element_text(family = "Helvetica", size = 9), 
	# 	axis.text.x = element_text(angle = 90,family = "Helvetica",  vjust = 0.5, hjust=1, face=labels_bold),
	# 	legend.title = element_text(size = 9, family = "Helvetica"),
	# 	legend.text = element_text(size = 7, family = "Helvetica"),
	# 	axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin()) +
	# 	facet_wrap(~Sample)
	# 	if (Sample != sort(unique(coords_ann$Sample))[1]){
	# 		p <- p + theme(axis.title.y=element_blank())
	# 	}
	# 	if(Sample != 'FS-0634-post'){
	# 		p + 	scale_y_continuous(labels = formatter(nsmall = 2))
	# 	}
	# 	plot_list_ratios[[Sample]] <- p
	# }

	# prop_5q <- ggplot(coords_ann[coords_ann$Sample!='SMD211420', ], aes(x=Cluster_names, fill=cell_5q)) + geom_bar() + facet_wrap(~Sample, nrow=1) + theme_classic()  +
	# 		scale_fill_manual(name = 'Cell type', values=c('#586994', '#A2ABAB')) + labs(y='Number of Cells') +
	# 		theme(legend.position='none',
	# 		axis.text = element_text(family = "Helvetica", size = 7), 
	# 		axis.title = element_text(family = "Helvetica", size = 9), 
	# 		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
	# 		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
	# 		legend.title = element_text(size = 7, family = "Helvetica"),
	# 		legend.text = element_text(size = 5, family = "Helvetica"))


	#### PLOT 1
		# pdf('./Plots/PaperFigures/FigPOST.pdf')
		# ggplot() + 
		# 		geom_point(coords_ann[coords_ann$cell_5q == 'normal', ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 1, alpha = 1)+ 
		# 		scale_color_manual(values=Cluster_colors) + #xlim(-13.5, 12) + ylim(-7.2, 8.5) +
		# 		geom_density_2d_filled(data= coords_ann[coords_ann$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)),bins=20, contour_var = "ndensity") + 
		# 		geom_density_2d(data= coords_ann[coords_ann$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2),bins=20, colour = "black", alpha = 0.2) + 
		# 		scale_discrete_manual("alpha", guide = "none", values = c(0, seq(0.1, 1, by=(1-(0.1))/(20-2))))  + 
		# 		scale_fill_viridis_d(option = "inferno") + theme_classic() + 
		# 		geom_density_2d(data= coords_ann[coords_ann$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.3) + 
		# 		theme(legend.position = 'none',  axis.title = element_text(family = "Helvetica", size = 7),  axis.text=element_blank(), axis.ticks=element_blank()) + 
		# 		facet_wrap(~Sample, nrow=2)

		# ggplot() + 
		# 		geom_point(coords_ann[coords_ann$cell_5q == 'normal', ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 1, alpha = 1)+ 
		# 		scale_color_manual(values=Cluster_colors) + 
		# 		geom_point(coords_ann[coords_ann$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2), size = 1, alpha = 0.7, shape=21, color='black', fill='#FCFFA4FF') +facet_wrap(~Sample, nrow=2)



		# print(prop_5q)

		# test <- as.matrix(-log(Enrichment_results_df))
		# myBreaks <- c(seq(min(test), 3, length.out=ceiling(9/2) + 1), 
		#               seq(max(test)/9, max(test), length.out=floor(9/2)))
		# pheatmap::pheatmap(-log(Enrichment_results_df), 
		# 				   cluster_rows = FALSE, 
		# 				   cluster_cols = FALSE, 
		# 				   show_rownames = TRUE, 
		# 				   show_colnames = TRUE, 
		# 				   cellheight = 10,
		# 				   cellwidth = 10,
		# 				   color = RColorBrewer::brewer.pal(name='Reds', n=9), 
		# 				   breaks=myBreaks)

		# cowplot::plot_grid(plotlist= plot_list_ratios)
		# dev.off()




coords_ann <- rbind(FetchData(post_data,     vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')), 
					FetchData(pre_post_data, vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')))


all_5q_depleted_cells_COPYKAT <- readRDS('./Data/CopyKat/all_5q_CopyKat_depleted_cells_POST.rds')
all_5q_depleted_cells_COPYKAT$Cell_id <- gsub('\\.', '-', all_5q_depleted_cells_COPYKAT$Cell_id)
all_5q_depleted_cells_CASPER <- readRDS('./Data/CASPER/all_5q_depleted_cells_POST.rds')

plotter_casper <- coords_ann
plotter_casper$is5q <- ifelse(rownames(plotter_casper) %in% c(all_5q_depleted_cells_CASPER,  paste(PrePost_CASPER$X1)), 'del5q', 'Other')
plt_casper <-  as.data.frame.matrix(table(plotter_casper$Sample, plotter_casper$is5q))
pct_casper <- round(plt_casper/rowSums(plt_casper), 3)*100

rownames(pct_casper) <- ifelse(rownames(pct_casper)== 'FS-0406-post', 'Partial Responder', 
						ifelse(rownames(pct_casper)== 'FS-0634-post', 'Complete Responder', 
						ifelse(rownames(pct_casper)== 'SMD132114579', 'Non Responder','UPS' )))

plotter_CopyKat <- coords_ann
plotter_CopyKat$is5q <- ifelse(rownames(plotter_casper) %in% c(all_5q_depleted_cells_COPYKAT$Cell_id, PrePost_Copykat$Cell_id), 'del5q', 'Other')
plt_CopyKat <-  as.data.frame.matrix(table(plotter_CopyKat$Sample, plotter_CopyKat$is5q))
pct_CopyKat <- round(plt_CopyKat/rowSums(plt_CopyKat), 3)*100


rownames(pct_CopyKat) <- ifelse(rownames(pct_CopyKat)== 'FS-0406-post', 'Partial Responder', 
						ifelse(rownames(pct_CopyKat)== 'FS-0634-post', 'Complete Responder', 
						ifelse(rownames(pct_CopyKat)== 'SMD132114579', 'Non Responder','UPS' )))


get_prop_tables_per_sample(post_data)
get_prop_tables_per_sample(MDS_data)
get_prop_tables_per_sample(pre_post_data)

deletion_percentages <- data.frame(Samples = c('Partial Responder',	'Complete Responder',	'Non Responder'), 
								Karyotype = c(50,0,100))
deletion_percentages <- merge(deletion_percentages, setNames(pct_casper[, 'del5q', drop=FALSE], c('CASPER')), by.x= 'Samples', by.y=0)
deletion_percentages <- merge(deletion_percentages, setNames(pct_CopyKat[, 'del5q', drop=FALSE], c('CopyKat')), by.x= 'Samples', by.y=0)

pct_selection <- reshape2::dcast(rbind(as.data.frame(get_prop_tables_per_sample(pre_post_data)), as.data.frame(get_prop_tables_per_sample(post_data))), Var1~Var2)
pct_selection <- setNames(pct_selection[pct_selection$Var1 != 'SMD211420', ], c('Samples', 'Selected\ncells', 'Other'))
pct_selection$Samples <- ifelse(pct_selection$Samples == 'FS-0406-post', 'Partial Responder', 
						 ifelse(pct_selection$Samples == 'FS-0634-post', 'Complete Responder', 
						 ifelse(pct_selection$Samples == 'SMD132114579', 'Non Responder','UPS')))

deletion_percentages <- merge(deletion_percentages, pct_selection[, c('Samples', 'Selected\ncells')], by= 'Samples')





deletion_percentages$Samples <- ifelse(deletion_percentages$Samples== 'Partial Responder', 'Partial Responder\nPatient_5', 
							   ifelse(deletion_percentages$Samples== 'Complete Responder', 'Complete Responder\nPatient_6', 
							   ifelse(deletion_percentages$Samples== 'Non Responder', 'Non Responder\nPatient_7', 'UPS' )))

MDS_sample_colors <- setNames(c('#d62828', '#f77f00', '#003049'), c('Partial Responder\nPatient_5',  'Complete Responder\nPatient_6', 'Non Responder\nPatient_7'))

plotter <- reshape2::melt(deletion_percentages)
plotter$variable <- factor(plotter$variable, levels = c('Karyotype', 'Selected\ncells', 'CASPER', 'CopyKat'))
selected_cells_Vs_Karyotype <- ggplot(plotter, aes(x = variable, y = value, group = Samples)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
			geom_hline(yintercept=deletion_percentages[, 'Karyotype'], linetype="dashed", color = "grey", size=0.5, alpha =0.5) +
			geom_line() + 
			geom_point(size = 2, aes(color = Samples)) + 
			ylim(0, 100) +
			scale_color_manual(values=MDS_sample_colors) +
			labs(y='Percentage of 5q cells') +
			theme_classic() + 
			guides(color = guide_legend(ncol=1, size=3, byrow=TRUE, title.position = 'top'))+
			theme(legend.position='right', 
				axis.text = element_text(family = "Helvetica", size = 7), 
				axis.title = element_text(family = "Helvetica", size = 9), 
				axis.text.x= element_text(angle = 45, hjust=1),
				axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
				legend.title = element_text(size = 7, family = "Helvetica"),
				legend.text = element_text(size = 6, family = "Helvetica"),
				legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm'))




coords_ann <- rbind(FetchData(post_data,     vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')), 
					FetchData(pre_post_data, vars=c('UMAP_1', 'UMAP_2', 'Sample', 'Cluster_names', 'cell_5q')))

coords_ann <- coords_ann[coords_ann$Sample != 'SMD211420',]

coords_ann$Sample <- ifelse(coords_ann$Sample == 'FS-0406-post', 'Partial Responder', 
					 ifelse(coords_ann$Sample == 'FS-0634-post', 'Complete Responder', 
					 ifelse(coords_ann$Sample == 'SMD132114579', 'Non Responder','UPS')))

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




get_geom_point_density <- function(plotter){
	xmin <- min(plotter$UMAP_1) - 3
	xmax <- max(plotter$UMAP_1) + 3
	ymin <- min(plotter$UMAP_2) - 3
	ymax <- max(plotter$UMAP_2) + 3 
	ggplot() + 
		geom_point(plotter[plotter$cell_5q == 'normal', ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 0.3, alpha = 1)+ 
		scale_color_manual(values=Cluster_colors) + xlim(xmin, xmax) + ylim(ymin, ymax) +
		geom_density_2d_filled(data= plotter[plotter$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)),bins=20, contour_var = "ndensity") + 
		geom_density_2d(data= plotter[plotter$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2),bins=20, colour = "black", alpha = 0.2) + 
		scale_discrete_manual("alpha", guide = "none", values = c(0, seq(0.1, 1, by=(1-(0.1))/(20-2))))  + 
		scale_fill_viridis_d(option = "inferno") + theme_void() + 
		geom_density_2d(data= plotter[plotter$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.3) + 
		theme(legend.position = 'none',  axis.title = element_blank(),  axis.text=element_blank(), axis.ticks=element_blank(), plot.title = element_text(hjust = 0.5, size = 8))
}

get_geom_point_Del_point <- function(plotter){
	xmin <- min(plotter$UMAP_1) - 3
	xmax <- max(plotter$UMAP_1) + 3
	ymin <- min(plotter$UMAP_2) - 2 
	ymax <- max(plotter$UMAP_2) + 2 
	ggplot() + 
			geom_point(plotter[plotter$cell_5q == 'normal', ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 0.3, alpha = 1)+ 
			scale_color_manual(values=Cluster_colors) + theme_void() + xlim(xmin, xmax) + ylim(ymin, ymax) +
			geom_point(plotter[plotter$cell_5q == 'del5q', ], mapping=aes(x=UMAP_1, y=UMAP_2), size = 1, alpha = 0.7, shape=21, color='black', fill='#FCFFA4FF') +
			theme(legend.position = 'none',  axis.title =element_blank(),  axis.text=element_blank(), axis.ticks=element_blank(),plot.title = element_text(hjust = 0.5, size = 8))
	}



plotter_selection <- coords_ann
plotter_selection$Sample <- ifelse(plotter_selection$Sample== 'Partial Responder', 'Partial Responder\nPatient_5', 
							ifelse(plotter_selection$Sample== 'Complete Responder', 'Complete Responder\nPatient_6', 
							ifelse(plotter_selection$Sample== 'Non Responder', 'Non Responder\nPatient_7', 'UPS' )))
plotter_selection$Sample <- factor(plotter_selection$Sample, levels=c('Partial Responder\nPatient_5', 'Complete Responder\nPatient_6', 'Non Responder\nPatient_7'))
plotter_selection$Cluster_names <- factor(plotter_selection$Cluster_names, levels=names(Cluster_colors))
plotter_selection$cell_5q <- ifelse(plotter_selection$cell_5q == 'del5q', 'del(5q)', 'non-del(5q)')
prop_5q <- ggplot(plotter_selection, aes(x=Cluster_names, fill=cell_5q)) + geom_bar(color='black') + facet_wrap(~Sample, nrow=1) + theme_classic()  +
		scale_fill_manual(name = 'Genotype', values=c('#586994', '#A2ABAB')) + labs(y='Number of Cells') +
		theme(
		legend.position='right',
		axis.text = element_text(family = "Helvetica", size = 7, color='black'), 
		axis.title = element_text(family = "Helvetica", size = 9), 
		axis.text.x= element_text(angle = 45, hjust=1),
		axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
		legend.title = element_text(size = 7, family = "Helvetica"),
		legend.text = element_text(size = 6, family = "Helvetica"),
		legend.key.size = unit(0.5, 'cm'), 
		strip.background = element_blank(), 
		strip.text = element_text(size =8, family = "Helvetica"), 
		panel.grid.major.y = element_line(color = "grey",size = 0.2,linetype = 2),
		panel.background = element_blank(),
		plot.background = element_blank())



Enrichment_results <- data.frame(Sample=NULL, Cluster=NULL, pval=NULL)
samples <- sort(unique(coords_ann$Sample))
samples <- samples[samples!='SMD211420']
for (Sample in samples){
	tmp_Smp <- coords_ann[coords_ann$Sample == Sample,]
	all_cells <- nrow(tmp_Smp)
	all_cells_5q <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'del5q',])
	all_cells_normal <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'normal',])

	for (CT in sort(unique(coords_ann$Cluster_names))){
		tmp_CT <- tmp_Smp[tmp_Smp$Cluster_names ==CT,]
		# phyper(q, m, n, k)
		pval_phyper <- phyper(sum(tmp_CT$cell_5q == 'del5q') -1 , all_cells_5q, all_cells_normal, nrow(tmp_CT), lower.tail=FALSE)
		Enrichment_results <- rbind(Enrichment_results, data.frame(Sample=Sample, Cluster=CT, pval=pval_phyper))
	}
}
Enrichment_results_df <- reshape2::dcast(Enrichment_results, Sample ~ Cluster, value.var = 'pval')
rownames(Enrichment_results_df) <-  Enrichment_results_df$Sample
Enrichment_results_df$Sample <- NULL


test <- as.matrix(-log(Enrichment_results_df))
myBreaks <- c(seq(min(test), 3, length.out=ceiling(9/2) + 1), 
              seq(max(test)/9, max(test), length.out=floor(9/2)))
mat <- as.matrix(-log(Enrichment_results_df))
enrichment_5q <- pheatmap::pheatmap(mat, 
				   cluster_rows = FALSE, 
				   cluster_cols = FALSE, 
				   show_rownames = TRUE, 
				   show_colnames = TRUE, 
				   cellheight = 10,
				   cellwidth = 10,
				   color = RColorBrewer::brewer.pal(name='Reds', n=9), 
				   breaks=myBreaks, silent=TRUE)


Enrichment_hm <-  t(as.matrix(-log(Enrichment_results_df)))
Enrichment_hm[Enrichment_hm > 10] <- 10

Enrichment_HM <- Heatmap(Enrichment_hm,
name = "-log \np.value",
col = circlize::colorRamp2(c(0, 1, 3, 5, 10), c("white", "white", "#FFF5F0", "#EF3B2C", "#67000D")),
heatmap_legend_param = list(
	title = "-log \np.value", at = c(0, 3, 5, 10), 
    labels = c('0', '3 (pval<0.05)' , '5', '10'),
	border = "grey"
),
rect_gp = gpar(col = "grey", lwd = 1),
cluster_rows = FALSE,
cluster_columns = FALSE, 
column_names_gp = grid::gpar(fontsize = 7, angle=45, fontfamily='Helvetica', color='black'),
row_names_gp = grid::gpar(fontsize = 7, fontfamily='Helvetica', color='black'),
width = ncol(Enrichment_hm)*unit(5, "mm"), 
height = nrow(Enrichment_hm)*unit(3, "mm"),
column_names_rot = 45) 




coords_ann$Cluster_names <-   factor(coords_ann$Cluster_names, levels=names(Cluster_colors))

legend <- cowplot::get_legend(get_geom_point_Del_point(coords_ann[coords_ann$Sample=='FS-0634-post', ]) + 
theme(legend.position = 'bottom', legend.title=element_blank(), legend.spacing.x = unit(1, 'mm'), legend.spacing.y = unit(1, 'mm')) + 
guides(colour = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=2.5))))
density_legend <-cowplot::get_legend( ggplot() + 
				geom_point(coords_ann[coords_ann$UMAP_1 <1 & coords_ann$UMAP_1>0, ], mapping=aes(x=UMAP_1, y=UMAP_2, color=UMAP_1), size = 1, alpha = 0.6)+ 
				viridis::scale_color_viridis(option = "inferno", limits=c(0,1), breaks =c(0,0.5,1), name= 'del(5q) cell\ndensity', alpha=0.95) + theme_classic()+
				# guides(fill=guide_legend(title="del(5q) cell\ndensity")) +
				theme(
					legend.position='right',  axis.title = element_text(family = "Helvetica", size = 9),
					legend.key.height= unit(0.6, 'cm'), legend.key.width= unit(0.2, 'cm'),
					axis.text=element_blank(), axis.ticks=element_blank(), legend.title=element_text(size=9, family='Helvetica')))

prop_5q_legend <- cowplot::get_legend(prop_5q)

pdf('./Plots/PaperFigures/Fig6.pdf')
cowplot::plot_grid(
	cowplot::plot_grid(
		cowplot::plot_grid(NULL, 
		cowplot::plot_grid( get_geom_point_density(coords_ann[coords_ann$Sample=='Partial Responder', ])   + ggtitle('Partial Responder\nPatient_5'), 
							get_geom_point_Del_point(coords_ann[coords_ann$Sample=='Complete Responder', ])+ ggtitle('Complete Responder\nPatient_6'),
							get_geom_point_density(coords_ann[coords_ann$Sample=='Non Responder',])   + ggtitle('Non Responder\nPatient_7'), 
							nrow=1),
		density_legend,
		nrow=1, rel_widths=c(0.1, 1, 0.2)),
		legend, 
	ncol=1, rel_heights = c(0.7, 0.1)), 
cowplot::plot_grid(prop_5q + theme(legend.position='none'), prop_5q_legend, ncol=2, rel_widths=c(1, 0.1)),
cowplot::plot_grid(selected_cells_Vs_Karyotype ,
	cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(Enrichment_HM,background = "transparent"))), 
	ncol=2, labels = c('C', 'D')),
# cowplot::plot_grid(enrichment_5q$gtable),
ncol=1, rel_heights=c(0.5, 0.4, 0.4), labels = c('A', 'B'))
dev.off()





# pdf('./Plots/PaperFigures/Figure7.pdf')
# legend <- cowplot::get_legend(get_geom_point_Del_point(coords_ann[coords_ann$Sample=='FS-0634-post', ]) + theme(legend.position = 'bottom')
# + guides(colour = guide_legend(override.aes = list(size=4))))

# cowplot::plot_grid(
# 	cowplot::plot_grid(
# 		cowplot::plot_grid( get_geom_point_density(coords_ann[coords_ann$Sample=='FS-0406-post', ])   + ggtitle('FS-0406-post'), 
# 							get_geom_point_Del_point(coords_ann[coords_ann$Sample=='FS-0634-post', ])+ ggtitle('FS-0634-post'),
# 							# get_geom_point_density(coords_ann[coords_ann$Sample=='SMD211420',])      + ggtitle('SMD211420'), 
# 							get_geom_point_density(coords_ann[coords_ann$Sample=='SMD132114579',])   + ggtitle('SMD132114579'), 
# 							nrow=1),
# 		legend, 
# 	ncol=1, rel_heights = c(0.5, 0.1)), 
# cowplot::plot_grid(
# 	prop_5q+ facet_wrap(~Sample, ncol=2, nrow=2),
# 	cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(Enrichment_HM))),
# 	ncol=2),
# # cowplot::plot_grid(enrichment_5q$gtable),
# ncol=1, rel_heights=c(0.5, 0.3, 0.2), labels = c('A', 'B', 'C'))
# selected_cells_Vs_Karyotype

# dev.off()


