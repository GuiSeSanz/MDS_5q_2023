
library(Seurat)
library(ggplot2)
library(ComplexHeatmap)

# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4

all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') # <- 
results_CK <- readRDS('./Data/CopyKat/results_CK.rds')
plotter_list <- readRDS('./Data/CopyKat/plotter_list.rds')
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')


real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
real_5q_cells <- all_5q_selected_cells[all_5q_selected_cells$Cell_id %in% real_5q_cells,]



SAMPLE_NAME <- '5qSamples'
sc_data <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', SAMPLE_NAME, '_Annotated_final.rds'))


Idents(sc_data) <- 'integrated_snn_res.0.8_sub'


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



MDS_sample_colors <- setNames(c('#495867', '#577399', '#bdd5ea','#f7f7ff'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
# MDS_sample_colors <- c('#355070', '#6d597a', '#b56576', '#e56b6f')
# MDS_sample_colors <- c('#f79f79', '#f7d08a', '#e3f09b', '#87b6a7')
# MDS_sample_colors <- c('#ed6a5a', '#f4f1bb', '#9bc1bc', '#5d576b')
# MDS_sample_colors <- c('#e8d6cb', '#d0ada7', '#ad6a6c', '#5d2e46')

Control_Elder_MDS_colors <- setNames(c('#002642','#840032', '#e59500'), c('Control', 'Elder', 'MDS'))

coords   <- as.data.frame(sc_data@reductions$umap@cell.embeddings)
coords_ann <- merge(setNames(sc_data[[c('integrated_snn_res.0.8_sub', 'Cluster_names')]], c('Cluster','Cluster_names')), coords, by=0)
coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')



point_legend <- cowplot::get_legend(ggplot(coords_ann[!coords_ann$Row.names %in% real_5q_cells$Cell_id, ], aes(x=UMAP_1, y=UMAP_2, color=Cluster_names)) + 
geom_point() + scale_color_manual(values=Cluster_colors) + 
guides(color = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=3, alpha=0.9))) + theme_classic() + 
theme(legend.position='bottom', axis.title = element_text(family = "Helvetica", size = 5), 
axis.text = element_text(family = "Helvetica", size = 5), 
axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
legend.title = element_blank(),
legend.text = element_text(size = 7, family = "Helvetica"), 
legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm')))



all_results <- ggplot() + 
				geom_point(coords_ann[!coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 0.6, alpha = 1)+ 
				scale_color_manual(values=Cluster_colors)+ xlim(-13.5, 12) + ylim(-7.2, 8.5) +
				geom_density_2d_filled(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)) , contour_var = "ndensity") +
				scale_discrete_manual("alpha",guide = "none", values =c(0, seq(0,0.8, by=0.1)))  + scale_fill_viridis_d(option = "inferno", guide='none') + theme_void()+
				geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.2) + 
				theme(legend.position='none',  #axis.title = element_text(family = "Helvetica", size = 7),  
				axis.title=element_blank(),
				axis.text=element_blank(), axis.ticks=element_blank(), legend.title=element_blank())



# perSample_results <- ggplot() + 
# 			geom_point(coords_ann[!coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 1, alpha = 1)+ 
# 			scale_color_manual(values=Cluster_colors) + xlim(-13.5, 12) + ylim(-9, 9) +
# 			geom_density_2d_filled(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)),bins=20, contour_var = "ndensity") + 
# 			geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2),bins=20, colour = "black", alpha = 0.2) + 
# 			scale_discrete_manual("alpha", guide = "none", values = c(0, seq(0.1, 1, by=(1-(0.1))/(20-2))))  + 
# 			scale_fill_viridis_d(option = "inferno") + theme_classic() + geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.3) + 
# 			theme(legend.position = 'none',  axis.title = element_text(family = "Helvetica", size = 7),  axis.text=element_blank(), axis.ticks=element_blank()) + 
# 			facet_wrap(~Sample, nrow=1)


plotter_selection <- coords_ann
plotter_selection$is_5q <- ifelse(plotter_selection$Row.names %in% real_5q_cells$Cell_id, '5q', 'Normal')
plt_selection <- as.data.frame.matrix(table(plotter_selection$Sample, plotter_selection$is_5q))
pct_selection <- round(plt_selection/rowSums(plt_selection), 3)*100


plotter_selection$Sample <- ifelse(plotter_selection$Sample == 'SMD34459', 'Patient_1', 
					 ifelse(plotter_selection$Sample == 'SMD35109', 'Patient_2',
					 ifelse(plotter_selection$Sample == 'SMD35303', 'Patient_3',
					 ifelse(plotter_selection$Sample == 'SMD37209', 'Patient_4','None'))))

plotter_selection$Cluster_names <- factor(plotter_selection$Cluster_names, levels=names(Cluster_colors))
plotter_selection$is_5q <- ifelse(plotter_selection$is_5q == '5q', 'del(5q)', 'non del(5q)')
prop_5q <- ggplot(plotter_selection, aes(x=Cluster_names, fill=is_5q)) + 
		geom_bar(, color='black') + facet_wrap(~Sample, nrow=1) + theme_classic()  +
		scale_fill_manual(name = 'Cell type', values=c('#586994', '#A2ABAB')) + labs(y='Number of Cells') +
		theme(legend.position='bottom',
		axis.text = element_text(family = "Helvetica", size = 7, color='black'), 
		axis.title = element_text(family = "Helvetica", size = 9, color='black'), 
		axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
		axis.title.x=element_blank(),
		legend.key.size = unit(5, 'mm'),
		legend.box="horizontal", legend.margin=margin(), 
		legend.title = element_text(size = 9, family = "Helvetica"),
		legend.text = element_text(size = 7, family = "Helvetica"),
		panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',colour = "grey"),
		panel.background = element_blank(),
		plot.background = element_blank(),
		strip.background = element_blank())



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
# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4


selected_cells_Vs_Karyotype <- ggplot(reshape2::melt(deletion_percentages), aes(x = variable, y = value, group = Samples)) + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10))+
			geom_hline(yintercept=deletion_percentages[, 'Karyotype'], linetype="dashed", color = "grey", size=0.5, alpha =0.5) +
			geom_line() + 
			geom_point(size = 2, aes(color = Samples)) + 
			ylim(0, 100) +
			scale_color_manual(values=MDS_sample_colors) +
			labs(y='Percentage of 5q cells') +
			theme_classic() + 
			guides(color = guide_legend(nrow=2, size=3, byrow=TRUE, title.position = 'top'))+
			theme(legend.position='bottom', 
			axis.text = element_text(family = "Helvetica", size = 7), 
			axis.title = element_text(family = "Helvetica", size = 9), 
			axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
			axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin(), 
			legend.title = element_text(size = 7, family = "Helvetica"),
			legend.text = element_text(size = 6, family = "Helvetica"),
			legend.spacing.x = unit(0, 'cm'), legend.spacing.y = unit(0, 'cm'))




data_5q_annotated <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_2_plot <- data_5q_annotated[[c('Sample', 'Cluster_names', 'cell_5q')]]
data_2_plot_5q     <- table(data_2_plot[data_2_plot$cell_5q == 'del5q', ])
data_2_plot_normal <- table(data_2_plot[data_2_plot$cell_5q == 'normal',])

Enrichment_results <- data.frame(Sample=NULL, Cluster=NULL, pval=NULL)
for (Sample in sort(unique(data_2_plot$Sample))){
	tmp_Smp <- data_2_plot[data_2_plot$Sample == Sample,]
	all_cells <- nrow(tmp_Smp)
	all_cells_5q <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'del5q',])
	all_cells_normal <- nrow(tmp_Smp[tmp_Smp$cell_5q == 'normal',])

	for (CT in sort(unique(data_2_plot$Cluster_names))){
		tmp_CT <- tmp_Smp[tmp_Smp$Cluster_names ==CT,]
		# phyper(q, m, n, k)
		pval_phyper <- phyper(sum(tmp_CT$cell_5q == 'del5q') -1 , all_cells_5q, all_cells_normal, nrow(tmp_CT), lower.tail=FALSE)
		Enrichment_results <- rbind(Enrichment_results, data.frame(Sample=Sample, Cluster=CT, pval=pval_phyper))
	}
}

formatter <- function(...){
  function(x) format(round(x, 1), ...)
}
plot_list_ratios <- list()
for (Sample in sort(unique(data_2_plot$Sample))){
	tmp <- data_2_plot[data_2_plot$Sample == Sample,]
	tmp <- table(tmp$Cluster_names, tmp$cell_5q)
	ratio_all <- sum(tmp[,1]) / sum(tmp[, 2])
	tmp <- setNames(as.data.frame(tmp[,1] / tmp[, 2]), c('Ratio'))
	tmp$Cluster_names <- factor(rownames(tmp))
	tmp$Sample <- Sample
	labels_bold <- Enrichment_results[Enrichment_results$Sample == Sample & Enrichment_results$pval < 0.05, 'Cluster']
	labels_bold <- ifelse(levels(tmp$Cluster_names) %in% labels_bold, 'bold', 'plain')
	p <- ggplot(tmp, aes(x=Cluster_names, y=Ratio, group = Sample)) + 
	geom_hline(yintercept = ratio_all, color='grey', linetype='dotted') +
	geom_point(color='dodgerblue4') + geom_line(color='dodgerblue4') +
	scale_y_continuous(labels = formatter(nsmall = 2))+
	theme_classic() + labs(x = "Cell types", y = "Ratio 5qCells/normalcells") +
	theme(legend.position='top', 
	axis.text = element_text(family = "Helvetica", size = 7), 
	axis.title = element_text(family = "Helvetica", size = 9), 
	axis.text.x = element_text(angle = 90,family = "Helvetica",  vjust = 0.5, hjust=1, face=labels_bold),
	legend.title = element_text(size = 9, family = "Helvetica"),
	legend.text = element_text(size = 7, family = "Helvetica"),
	axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin()) +
	facet_wrap(~Sample)
	if (Sample != sort(unique(data_2_plot$Sample))[1]){
		p <- p + theme(axis.title.y=element_blank())
	}
	plot_list_ratios[[Sample]] <- p
}

proliferations_ratios <- FetchData(data_5q_annotated, c('Cluster_names', 'cell_5q', 'S.Score', 'G2M.Score', 'Phase'))

proliferations_ratios$Sample <- stringr::str_extract(rownames(proliferations_ratios), '[\\w]+(?=_)')
data <- reshape2::melt(proliferations_ratios[proliferations_ratios$Cluster_names %in% c('EarlyErythroid', 'LateErythroid', 'GMP', 'MK_Prog', 'LMPP', 'MEP'),])
data$Cluster_names <-  stringr::str_replace(data$Cluster_names, '(?<=\\w)E', ' E')
proliferations_ratios_plot <-  ggplot(data, aes(x=Cluster_names, y= value, fill=cell_5q)) + 
 geom_boxplot(alpha=0.9) + 
labs(y='Proliferation Scores') +
ggpubr::stat_compare_means(aes(group=cell_5q), label = "p.signif", method="wilcox.test", hide.ns=FALSE, paired=F, label.y=1.3) + 
theme_classic() +  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype')  + 
scale_x_discrete(labels=scales::label_wrap(10)) + 
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.title.x = element_blank(),
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, family = "Helvetica", size = 7),
		axis.title = element_text(size=9, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text( family = "Helvetica",size=7),
		legend.title=element_text( family = "Helvetica",size=9),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'right', strip.placement = "outside")

# pdf('./Plots/PaperFigures/Test.pdf.pdf', width=7, height=7)
# proliferations_ratios_plot + facet_wrap(Sample~variable, ncol=2)
# dev.off()

results <- data.frame(SMP=NULL,CT=NULL,Score=NULL, pval=NULL)
for (SMP in unique(data$Sample)){
	for (CT in unique(data$Cluster_names)){
		for (Score in c('G1', 'G2M', 'S')){
			tmp <- data[data$Sample == SMP & data$Cluster_names == CT & data$Phase == Score,]
			if(length(unique(tmp$cell_5q)) < 2){
				next
			}
			wt <- wilcox.test(tmp[tmp$cell_5q == 'normal', 'value'], tmp[tmp$cell_5q == 'del5q', 'value'], paired = F)
			results <- rbind(results, data.frame(SMP, CT, Score, pval=wt$p.value))
		}

	}
}

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(data_5q_annotated) <- 'SCT'
data_5q_annotated_bis <- CellCycleScoring(data_5q_annotated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


tmp <- coords_ann[!coords_ann$Row.names %in% real_5q_cells$Cell_id & coords_ann$UMAP_1 <1 & coords_ann$UMAP_1>0, ]
density_legend <-cowplot::get_legend( ggplot() + 
				geom_point(tmp, mapping=aes(x=UMAP_1, y=UMAP_2, color=UMAP_1), size = 1, alpha = 0.6)+ 
				viridis::scale_color_viridis(option = "inferno", limits=c(0,1), breaks =c(0,0.5,1), name= 'del(5q) cell\ndensity', alpha=0.95) + theme_classic()+
				# guides(fill=guide_legend(title="del(5q) cell\ndensity")) +
				theme(
					legend.position='right',  axis.title = element_text(family = "Helvetica", size = 9),
					legend.key.height= unit(0.6, 'cm'), legend.key.width= unit(0.2, 'cm'),
					axis.text=element_blank(), axis.ticks=element_blank(), legend.title=element_text(size=9, family='Helvetica')))



Enrichment_results_df <- reshape2::dcast(Enrichment_results, Sample ~ Cluster, value.var = 'pval')
rownames(Enrichment_results_df) <-  Enrichment_results_df$Sample
Enrichment_results_df$Sample <- NULL


test <- as.matrix(-log(Enrichment_results_df))
myBreaks <- c(seq(min(test), 3, length.out=ceiling(9/2) + 1), 
              seq(max(test)/9, max(test), length.out=floor(9/2)))
enrichment_5q <- pheatmap::pheatmap(-log(Enrichment_results_df), 
				   cluster_rows = FALSE, 
				   cluster_cols = FALSE, 
				   show_rownames = TRUE, 
				   show_colnames = TRUE, 
				   cellheight = 10,
				   cellwidth = 10,
				   color = RColorBrewer::brewer.pal(name='Reds', n=9), 
				   breaks=myBreaks, silent=TRUE)


enrichment_5q_CH <- Heatmap(-log(Enrichment_results_df), name='-log(Pvalue)',
col = RColorBrewer::brewer.pal(9,"Reds"), 
rect_gp = gpar(col = "grey", lwd = 1),
cluster_rows = FALSE,
cluster_columns = FALSE
) 

# SMD34459 <- Patient_1
# SMD35109 <- Patient_2
# SMD35303 <- Patient_3
# SMD37209 <- Patient_4

Enrichment_hm <-  t(as.matrix(-log(Enrichment_results_df)))
Enrichment_hm[Enrichment_hm > 10] <- 10
colnames(Enrichment_hm) <- c('Patient_1', 'Patient_2', 'Patient_3', 'Patient_4')
Enrichment_hm <- Enrichment_hm[c(names(Cluster_colors)),]

Enrichment_HM <- Heatmap(Enrichment_hm,
col = circlize::colorRamp2(c(0, 1, 3, 5, 10), c("white", "white", "#FFF5F0", "#EF3B2C", "#67000D")),
heatmap_legend_param = list(
	title_gp = grid::gpar(fontsize=8, fontface="bold"),
	labels_gp=grid::gpar(fontsize=7),
	legend_height = unit(3, "cm"),
	title = "-log \np.value", at = c(0, 3, 5, 10), 
    labels = c('0', '3\n(pval<0.05)' , '5', '10'),
	border = "black"
),
# name = "-log \np.value",
rect_gp = gpar(col = "grey", lwd = 1),
cluster_rows = FALSE,
cluster_columns = FALSE, 
column_names_gp = grid::gpar(fontsize = 7, angle=45),
row_names_gp = grid::gpar(fontsize = 7),
width = ncol(Enrichment_hm)*unit(5, "mm"), 
column_names_rot = 45) 


coords_ann$Sample <- ifelse(coords_ann$Sample == 'SMD34459', 'Patient_1', 
					 ifelse(coords_ann$Sample == 'SMD35109', 'Patient_2',
					 ifelse(coords_ann$Sample == 'SMD35303', 'Patient_3',
					 ifelse(coords_ann$Sample == 'SMD37209', 'Patient_4','None'))))

perSample_results <- ggplot() + 
			geom_point(coords_ann[!coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, color=Cluster_names), size = 0.3, alpha = 1)+ 
			scale_color_manual(values=Cluster_colors) + xlim(-13.5, 12) + ylim(-9, 9) +
			geom_density_2d_filled(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2, alpha = after_stat(level)),bins=20, contour_var = "ndensity") + 
			geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2),bins=20, colour = "black", alpha = 0.2) + 
			scale_discrete_manual("alpha", guide = "none", values = c(0, seq(0.1, 1, by=(1-(0.1))/(20-2))))  + 
			scale_fill_viridis_d(option = "inferno") + theme_classic() + geom_density_2d(data= coords_ann[coords_ann$Row.names %in% real_5q_cells$Cell_id, ], mapping=aes(x=UMAP_1, y=UMAP_2), colour = "black", alpha = 0.3) + 
			theme(legend.position = 'none',  axis.text=element_blank(), axis.ticks=element_blank(), 
			axis.title=element_blank(), line = element_blank(), strip.background = element_blank()) + 
			facet_wrap(~Sample, ncol=2) 





pdf(paste0('./Plots/PaperFigures/Fig3_axisGood.pdf'), width=8.25, height=6.75)
cowplot::plot_grid(
	cowplot::plot_grid(
		cowplot::plot_grid(
				all_results, 
				perSample_results,
				density_legend, 
		ncol=3, rel_widths = c(0.6,0.4, 0.1), labels=c('A', 'B', '')),
	point_legend, 
	nrow=2, rel_heights=c(0.9,0.1)),
cowplot::plot_grid(prop_5q, NULL, 
					cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(Enrichment_HM))),
					ncol=3, rel_widths=c(1, 0.1,0.4), labels=c('C', 'D', '')), 
nrow=2, rel_heights=c(0.5,0.5),align = 'h', axis = "lr")
dev.off()




