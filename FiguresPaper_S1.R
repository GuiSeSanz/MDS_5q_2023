
library(Seurat)
library(ggplot2)


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
		guides(color = guide_legend(nrow=legend_row, byrow=TRUE, override.aes = list(size=3, alpha=0.9)))
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

elder_integrated <-  readRDS(paste0(getwd(), '/Data/','elder_Samples_Annotated_final.rds'))
coords_ann <- FetchData(elder_integrated, vars=c('UMAP_1', 'UMAP_2', 'Cluster_names', 'Sample'))

features_2_plot <- c('AVP', 'HOPX', 'CRHBP', 'SATB1', 'LSP1', 'FLT3','MPO', 'CTSG', 'PRTN3','AZU1', 'ELANE','FCER1G', 'CSTA', 'LYZ', 'IL3RA', 'IRF7', 'IRF8','JCHAIN', 'IKZF1','DNTT', 'ADA', 'LTB','CD79A', 'EBF1', 'VPREB1','PBX1', 'PLEK', 'DAD1','CNRIP1', 'SLC40A1', 'PKIG','KLF1', 'AHSP', 'HBB','RUNX1', 'HDC', 'MS4A3')

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

Idents(elder_integrated) <- 'Cluster_names'
DefaultAssay(elder_integrated) <- 'SCT'



sample_umaps <- list()
for (Sample in  unique(elder_integrated$Sample)){
	sample_name <- ifelse(Sample == 'GSM5460411', 'Healthy_1', 
 				   ifelse(Sample == 'GSM5460412', 'Healthy_2',
 				   ifelse(Sample == 'GSM5460413', 'Healthy_3', 'UPS')))
	tmp <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/Normal_Data/', Sample, '_seurat_obj_norm.rds'))
	tmp_coords <- as.data.frame(tmp@reductions$umap@cell.embeddings)
	tmp_coords <- tmp_coords[rownames(tmp_coords) %in% rownames(coords_ann),]
	tmp_coords <- merge(tmp_coords, coords_ann[,'Cluster_names', drop=FALSE], by=0)
	sample_umaps[[Sample]] <- get_umap(tmp_coords, 'Cluster_names', Cluster_colors) + 
							  theme(legend.position= 'none', plot.title = element_text(hjust = 0.5)) + 
							  ggtitle(sample_name)
}

coords_ann$Sample <- ifelse(coords_ann$Sample == 'GSM5460411', 'Healthy_1', 
 					 ifelse(coords_ann$Sample == 'GSM5460412', 'Healthy_2',
 					 ifelse(coords_ann$Sample == 'GSM5460413', 'Healthy_3', 'UPS')))


ElderSample_colors <- c( '#264653', '#2a9d8f', '#e9c46a')
names(ElderSample_colors) <- unique(coords_ann$Sample)


patient_dist_cell_count <- ggplot(coords_ann, aes(x=Cluster_names, fill=Sample)) + 
		geom_bar( position='stack', colour="black") + 
		theme_classic() + scale_fill_manual(values=ElderSample_colors, name='Samples') + 
		theme(axis.text=element_text(family = "Helvetica", size=7),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=7),
			legend.title=element_text(size=7),
			legend.key.size=unit(0.5, "cm"),
			panel.grid.major.y = element_line( size=.1, color="grey" ),
			legend.position = 'none' ) + 
		labs(y = "Number of Cells")


tmp_data <- reshape2::melt(table(coords_ann$Sample, coords_ann$Cluster_names))
tmp_data$Var2 <- factor(tmp_data$Var2, levels=c(rev(names(Cluster_colors))))
pct_cluster_per_sample <- ggplot(tmp_data,aes(x= value, y= Var1, fill= Var2, ,order=Var2)) + 
geom_bar(position="fill", stat="identity", colour="black") + scale_fill_manual(values=Cluster_colors, name='Samples') +
theme_classic() + 
theme(axis.text=element_text(family = "Helvetica", size=7),
	# axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
	axis.title.y = element_blank(),
	axis.title = element_text(size=7, family = "Helvetica"),
	axis.line = element_line(colour = 'black', size = 0.5),
	legend.text=element_text(size=7),
	legend.title=element_text(size=9),
	legend.key.size=unit(0.5, "cm"),
	panel.grid.major.y = element_line( size=.1, color="grey" ),
	legend.position = 'none' ) + 
labs(x = "Percentage of cells")



# sample_umaps_mod <- sample_umaps
# sample_umaps_mod[['dist']] <- pct_cluster_per_sample

# pdf('./Plots/PaperFigures/FigS1.pdf', width=8.25, height=11.75)
# cowplot::plot_grid(
# 	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(coords_ann$Cluster_names) ], 3) ,
# 					cowplot::plot_grid(plotlist=sample_umaps_mod, 
# 					nrow=2, labels=c('C')), 
# 	ncol=2, labels=c('A', 'B')),
# 	cowplot::plot_grid(
# 			get_dotplot(elder_integrated , features_2_plot), 
# 			patient_dist_cell_count + theme(legend.position='bottom'), 
# 		ncol=2, labels=c('', 'E')),
# nrow=2, rel_heights = c(1, 0.8))
# dev.off()



sample_umaps_mod <- sample_umaps
sample_umaps_mod[['legend']] <-cowplot::get_legend(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(coords_ann$Cluster_names) ], 8))

pdf('./Plots/PaperFigures/FigS1.pdf', width=8.25, height=11.75)
cowplot::plot_grid(
	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(coords_ann$Cluster_names) ], 3) + theme(legend.position='none'),
					cowplot::plot_grid(plotlist=sample_umaps_mod, 
					nrow=2), 
	ncol=2, labels=c('A', 'B')),
	cowplot::plot_grid(
			get_dotplot(elder_integrated , features_2_plot),
			cowplot::plot_grid(pct_cluster_per_sample,
			patient_dist_cell_count + theme(legend.position='bottom'),
			nrow=2, labels=c('D', 'E')), 
		ncol=2, labels=c('C', '')),
nrow=2, rel_heights = c(1, 0.8))
dev.off()








# # full size figure width=8.3, height=11.7
# sample_legend <- cowplot::get_legend(patient_dist_cell_count + theme(legend.position='bottom') + 
# guides(fill = guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=0.5, alpha=0.9))))
# elder_mds_legend <- cowplot::get_legend(barplot_populations + theme(legend.position='bottom') + 
# guides(fill = guide_legend(nrow=1, byrow=TRUE, override.aes = list(size=0.5, alpha=0.9))))
# pdf('./Plots/PaperFigures/Fig1.pdf', width=8.3, height=10)
# cowplot::plot_grid(
# 	cowplot::plot_grid(get_umap(coords_ann, 'Cluster_names', Cluster_colors[names(Cluster_colors) %in% unique(sc_data$Cluster_names) ], 4), 
# 					get_dotplot(sc_data , features_2_plot), 
# 	ncol=2, labels=c('A', 'B')),
# 	cowplot::plot_grid(
# 		cowplot::plot_grid(patient_dist_cell_count, pct_cluster_per_sample, barplot_populations, ncol=3, labels=c('C', 'D', 'E')), 
# 		cowplot::plot_grid(sample_legend, elder_mds_legend , ncol=2, rel_widths=c(2,1)), 
# 		nrow=2, rel_heights=c(1,0.2)),
# nrow=2)
