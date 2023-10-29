
library(Libra)
library(Seurat)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(colorRamp2)
library(enrichR)

get_enrichr_results <- function(data_set){
	all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
	for (CT in names(data_set)){
		message(CT)
		tmp <- data_set[[CT]]
		tmp_pos <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >0,'gene']
		tmp_neg <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC <0,'gene']
		enriched <- enrichr(tmp_pos, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
		if(nrow(enriched) > 1){
			all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
		}
		enriched <- enrichr(tmp_neg, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
		if(nrow(enriched) > 1){
			all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))
		}
	}
	return(all_results)
}


get_SimiC_densisites <- function(regulon, data_auc=df_auc){
	p <- ggplot(data_auc[data_auc$driver ==regulon,], aes(x=value, y=..scaled.., fill=TimePoint)) + #xlim(0,1) + 
	geom_density(alpha = 0.6, adjust = 1/2) + theme_classic() + 
	scale_fill_manual(values=color_map) + scale_y_continuous(breaks=seq(0,1,0.5)) +
	theme(legend.position = 'none') + 
	labs(y=regulon, x='Activity score') + 
	theme(axis.title.x=element_blank(), 
	axis.title.y=element_text(size=9, face='bold'))
	return(p)
}


cell_type_order <- c(
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
#



sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
# LIBRA uses the raw counts

expression <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
metadata <-	 setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type')
metadata$label <- sc_data_MDS5q_subseted$cell_5q
metadata$replicate <- sc_data_MDS5q_subseted$Sample

metadata$label <- factor(metadata$label, levels = c('del5q', 'normal'))
metadata <- metadata[rownames(metadata) %in% colnames(expression), ]

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_plotter <- DE[DE$gene %in% c('PRSS21', 'ZCCHC10', 'FIGN', 'SIL1',  'CCL5', 'MAP3K7CL',  'IGKC'), c('cell_type', 'gene', 'avg_logFC', 'p_val_adj')]



DE_plotter_avg_logFC <- reshape2::dcast(DE_plotter[, c('cell_type', 'gene', 'avg_logFC')], cell_type ~ gene)
DE_plotter_p_val_adj <- reshape2::dcast(DE_plotter[, c('cell_type', 'gene', 'p_val_adj')], cell_type ~ gene)
rownames(DE_plotter_p_val_adj) <- DE_plotter_p_val_adj$cell_type
rownames(DE_plotter_avg_logFC) <- DE_plotter_avg_logFC$cell_type
DE_plotter_p_val_adj <- DE_plotter_p_val_adj[, -1]
DE_plotter_avg_logFC <- DE_plotter_avg_logFC[, -1]

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm =TRUE)
  breaks[!duplicated(breaks, na.rm=TRUE)]
}
mat_breaks <- quantile_breaks(-log(DE_plotter_p_val_adj), n = 11)
 
# pdf('./Plots/DE_5qVsNon5q.pdf')
# log_Pval <- -log(DE_plotter_p_val_adj)
# log_Pval[log_Pval>10] <- 10
# print(pheatmap::pheatmap(t(log_Pval),cellwidth = 20, cellheight = 20,
# cluster_rows=FALSE, cluster_cols=FALSE, main='-log10(p_val_adj)', 
# color = colorRampPalette(c("white", "yellow", "red"))(50)))

# print(pheatmap::pheatmap(t(DE_plotter_avg_logFC),cellwidth = 20, cellheight = 20,
# cluster_rows=FALSE, cluster_cols=FALSE, main='avg_logFC', 
# color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="Spectral")))(100)))
# dev.off()

# HM_DE <- pheatmap::pheatmap(t(log_Pval),cellwidth = 20, cellheight = 20,
# cluster_rows=FALSE, cluster_cols=FALSE, main='-log10(p_val_adj)', 
# color = colorRampPalette(c("white", "yellow", "red"))(50))


terms_2_plot <- c("DNA replication (GO:0006260)", "DNA metabolic process (GO:0006259)", "regulation of transcription involved in G1/S transition of mitotic cell cycle (GO:0000083)", "DNA replication initiation (GO:0006270", "mitotic sister chromatid segregation (GO:0000070)", "mitotic DNA replication (GO:1902969)", "DNA replication checkpoint signaling (GO:0000076)", "mitotic chromosome condensation (GO:0007076)", "mitotic nuclear division (GO:0140014)", "DNA-dependent DNA replication (GO:0006261)", "DNA strand elongation involved in DNA replication (GO:0006271)","peptide biosynthetic process (GO:0043043)", "cytoplasmic translation (GO:0002181)", "translation (GO:0006412)", "cellular protein metabolic process (GO:0044267)", "ribosome biogenesis (GO:0042254)", "rRNA processing (GO:0006364)", "mitochondrial translational elongation (GO:0070125)", "mitochondrial translational termination (GO:0070126)", "translational termination (GO:0006415)", "mitochondrial translation (GO:0032543)", "positive regulation of apoptotic process (GO:0043065)","intrinsic apoptotic signaling pathway (GO:0097193)", 'p53 signaling pathway','DNA repair (GO:0006281)')

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human')

all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')
all_results_5qVsElder <- get_enrichr_results(all_results_psuedo_5qVsElder)


all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')
all_results_Non5qVsElder <- get_enrichr_results(all_results_psuedo_Non5qVsElder)


all_results_Non5qVsElder$Comparison <- 'non-del(5q) Vs Healthy'
all_results_5qVsElder$Comparison    <- 'del(5q) Vs Healthy'
all_results <- rbind(all_results_5qVsElder, all_results_Non5qVsElder)

tmp <- all_results
tmp[tmp$Combined.Score >200, 'Combined.Score'] <- 200
tmp$Combined.Score_size <- tmp$Combined.Score
tmp <- tmp[tmp$Term %in% terms_2_plot, ]
tmp$Term <- stringr::str_remove(tmp$Term,'\\([^)]*\\)')
tmp[tmp$is_pos == FALSE, 'Combined.Score'] <- tmp[tmp$is_pos == FALSE, 'Combined.Score']*-1
# get_order
for (term in unique(tmp$Term)){
	tmp[tmp$Term == term, 'Combined.Score_OVERALL'] <- sum(tmp[tmp$Term == term, 'Combined.Score'])
}
sort_order <-  unique(tmp[,c('Term', 'Combined.Score_OVERALL') ])
sort_order <- sort_order[order(sort_order$Combined.Score_OVERALL, decreasing=FALSE),'Term']
tmp$Term <- factor(tmp$Term, levels = sort_order)

tmp$cell_type <- factor(tmp$cell_type, levels=cell_type_order)

dotplot_DE <- ggplot() + 
geom_point(tmp, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) +
facet_wrap(~Comparison) +  
guides(fill=guide_colorbar(title='Combined\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
theme(
	legend.position='bottom',
	axis.text.y=element_text(size=8, color='black'),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size = 8, angle=90, vjust=1, hjust=1, color='black'),
	legend.key.size = unit(0.3, "cm"),
	legend.title=element_text(size=9, face='bold'), 
	legend.text=element_text(size=7))


log_Pval <- -log(DE_plotter_p_val_adj)
log_Pval[log_Pval>10] <- 10
log_Pval <- log_Pval[cell_type_order,]

logs_pval_res <- Heatmap(as.matrix(log_Pval),
col = RColorBrewer::brewer.pal(9,"Blues"),name = "-log \np.value",
rect_gp = gpar(col = "grey", lwd = 1),
cluster_rows = FALSE,
cluster_columns = FALSE, 
column_names_gp = grid::gpar(fontsize = 7, angle=45),
row_names_gp = grid::gpar(fontsize = 7),
width = ncol(log_Pval)*unit(10, "mm"), 
column_names_rot = 45) 


df_auc <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_df_auc_ThreeWay.rds')
MinMax <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_MinMax_clust_ThreeWay.rds')
# MinMax <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_MinMax_clust_only_5q_non5q.rds')

phenotypes <- c('del(5q)', 'non-del(5q)', 'Healthy')
color_map <- c('#003049', '#d62828', '#f77f00')
names(color_map) <- phenotypes

# Plot the heatmap of the regulatory dissimilarity score

df_auc$TimePoint <- ifelse(df_auc$TimePoint == '5q', 'del(5q)',
					ifelse(df_auc$TimePoint == 'elder', 'Healthy', 
					ifelse(df_auc$TimePoint == 'non5q', 'non-del(5q)','UPS')))


# p <- pheatmap::pheatmap(MinMax,color=viridis::viridis(50, direction = 1, option = "C"), fontsize=5, angle_col =45, cellwidth=10, silent=TRUE, treeheight_row=10, treeheight_col=0)

# TF_2_show <- c('ZNF451','YBX1','PSPC1','MIER1','ARID4A','KAT6B', 'RERE', 'KDM2A')
# ha = rowAnnotation(foo = anno_mark(at = match(TF_2_show,rownames(MinMax)), 
#     labels = TF_2_show))

col_fun = colorRamp2(seq(1:50), viridis::viridis(50, direction = 1, option = "C"))

lgd = Legend(at = c(0, 1), col_fun = col_fun, 
	title = "Activity\nDisimilarity\nScore", 
    labels = c("0", "1"))



hm <- ComplexHeatmap::Heatmap(as.matrix(MinMax), name='SimiC3W', col=viridis::viridis(50, direction = 1, option = "inferno")[1:48], 
width = ncol(MinMax)*unit(4.2, "mm"), 
height = nrow(MinMax)*unit(1.35, "mm"),
show_column_dend = FALSE, 
show_row_dend = FALSE, 
# row_dend_width = unit(0.4, "cm"),
column_names_gp = grid::gpar(fontsize = 6, angle=45),
row_names_gp = grid::gpar(fontsize = 4.5),
# row_names_rot = 45,
# column_title_rot = 45,
# right_annotation = ha, 
heatmap_legend_param=list(
	at = c(0, 1), col_fun = col_fun, 
	title = "Activity\nDisimilarity\nScore", 
	# direction = "vertical", 
	# title_position = "right",
	labels_gp = gpar(fontsize = 8),
	title_gp = gpar(fontsize = 8),
	labels = c("0", "1"), 
	legend_height = unit(3, "cm"),
	legend_width = unit(0, "cm")),
show_row_names=TRUE)


simic_legend <- cowplot::get_legend(get_SimiC_densisites('JARID2') + 
theme(legend.position='right')+ 
guides(fill=guide_legend(nrow=2,byrow=TRUE, title='Condition')))


log_Pval <- log_Pval[cell_type_order,]
logs_pval_res <- Heatmap(as.matrix(log_Pval),
col = circlize::colorRamp2(c(0, 1, 3, 5, 10), c("white", "white", "#f5d1bf", "#EF3B2C", "#67000D")),
heatmap_legend_param = list(
	direction = "horizontal",
	title = "-log Adj.Pvalue", at = c(0, 3, 5, 10), 
	labels = c('0', '3\n(Adj.Pval<0.05)' , '5', '10'),
	border = "black"
),
rect_gp = gpar(col = "grey", lwd = 1),
cluster_rows = FALSE,
cluster_columns = FALSE, 
column_names_gp = grid::gpar(fontsize = 8, angle=45, fontface='italic'),
row_names_gp = grid::gpar(fontsize = 8),
width  = ncol(log_Pval)*unit(4, "mm"), 
height = nrow(log_Pval)*unit(6, "mm"), 
column_names_rot = 90) 


pdf('./Plots/PaperFigures/Fig4.pdf', width=8.3, height=11.7)
cowplot::plot_grid(
	cowplot::plot_grid(
					cowplot::plot_grid(grid::grid.grabExpr(ComplexHeatmap::draw(logs_pval_res, 
					heatmap_legend_side = "bottom")))
					,
					dotplot_DE, 
			ncol=2, rel_widths = c(0.3, 0.7), labels=c('A', 'B')),
	cowplot::plot_grid(
		cowplot::plot_grid(
				get_SimiC_densisites('ZNF451') ,
				get_SimiC_densisites('YBX1')   ,
				get_SimiC_densisites('PSPC1')  ,
		ncol=3),
		cowplot::plot_grid(
				get_SimiC_densisites('JARID2') ,
				get_SimiC_densisites('IRF1')   ,
				get_SimiC_densisites('KAT6B')  ,
		ncol=3),
		cowplot::plot_grid(
				get_SimiC_densisites('RERE')  + theme(axis.title.x=element_text(size=9)),
				get_SimiC_densisites('KDM2A') + theme(axis.title.x=element_text(size=9)), 
				simic_legend, 
		ncol=3),
	nrow=3, labels=c('C', 'D', 'E')),
	nrow=2, rel_heights=c(0.55,0.45))
dev.off()