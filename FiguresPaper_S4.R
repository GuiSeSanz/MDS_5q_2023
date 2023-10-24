
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(colorRamp2)
library(enrichR)



df_auc <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_df_auc_ThreeWay.rds')
MinMax <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_MinMax_clust_ThreeWay.rds')
weigths <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_Weigths_ThreeWay.rds')
# MinMax <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_MinMax_clust_only_5q_non5q.rds')

phenotypes <- c('del(5q)', 'non-del(5q)', 'Healthy')
color_map <- c('#003049', '#d62828', '#f77f00')
names(color_map) <- phenotypes

# Plot the heatmap of the regulatory dissimilarity score

df_auc$TimePoint <- ifelse(df_auc$TimePoint == '5q', 'del(5q)',
					ifelse(df_auc$TimePoint == 'elder', 'Healthy', 
					ifelse(df_auc$TimePoint == 'non5q', 'non-del(5q)','UPS')))

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021")


get_SimiC_densisites <- function(regulon, data_auc=df_auc){
	p <- ggplot(data_auc[data_auc$driver ==regulon,], aes(x=value, y=..scaled.., fill=TimePoint)) + #xlim(0,1) + 
	geom_density(alpha = 0.6, adjust = 1/2) + theme_classic() + 
	scale_fill_manual(values=color_map) + scale_y_continuous(breaks=seq(0,1,0.5)) +
	theme(legend.position = 'none') + 
	labs(y=regulon, x='Activity score') + 
	theme(axis.title.x=element_blank(), 
	axis.title.y=element_text(size=8))
	return(p)
}




weigths_KDM2A <- as.character(weigths[weigths$driver == 'KDM2A' & weigths$value !=0,'target'])
enrichr_KDM2A <- enrichr(weigths_KDM2A, dbs)
weigths_RERE <-  as.character(weigths[weigths$driver == 'RERE'  & weigths$value !=0,'target'])
enrichr_RERE  <- enrichr(weigths_RERE, dbs)

enrichr_KDM2A <- do.call("rbind", enrichr_KDM2A)
enrichr_RERE  <- do.call("rbind", enrichr_RERE)

enrichr_KDM2A <- enrichr_KDM2A[enrichr_KDM2A$Adjusted.P.value < 0.05, ]
enrichr_RERE <- enrichr_RERE[enrichr_RERE$Adjusted.P.value < 0.05, ]
enrichr_KDM2A <- enrichr_KDM2A[order(enrichr_KDM2A$Adjusted.P.value), ]
enrichr_RERE <- enrichr_RERE[order(enrichr_RERE$Adjusted.P.value), ]

enrichr_KDM2A[enrichr_KDM2A$Combined.Score > 200, 'Combined.Score'] <- 200 
enrichr_RERE[enrichr_RERE$Combined.Score > 200, 'Combined.Score'] <- 200 

enrichr_RERE$Combined.Score_size <- enrichr_RERE$Combined.Score
enrichr_KDM2A$Combined.Score_size <- enrichr_KDM2A$Combined.Score

enrichr_RERE$Term <-  stringr::str_remove(enrichr_RERE$Term,'\\([^)]*\\)')
enrichr_KDM2A$Term <-  stringr::str_remove(enrichr_KDM2A$Term,'\\([^)]*\\)')

enrichr_RERE$TF <- 'RERE'
enrichr_KDM2A$TF <- 'KDM2A'

tmp <- rbind(enrichr_KDM2A[1:15,], enrichr_RERE[1:15,])
tmp$Term <- factor(tmp$Term, levels=unique(tmp$Term))
tmp$TF <- factor(tmp$TF, levels=c('RERE', 'KDM2A'))

dotplot_TFs <- ggplot() + 
geom_point(tmp, mapping= aes(y=Term, x=TF, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) + 
guides(fill=guide_colorbar(title='Combined\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
theme(axis.text.y=element_text(size=8),
	axis.title.x=element_blank(),
	axis.text.x=element_text(size = 8, angle=90, vjust=1, hjust=1),
	legend.position='right',
	legend.key.size = unit(0.3, "cm"),
	legend.title=element_text(size=8), 
	legend.text=element_text(size=6))





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


simic_legend <- cowplot::get_legend(get_SimiC_densisites('JARID2') + theme(legend.position='bottom'))

# pdf('./Plots/PaperFigures/Fig5.pdf', height=7.5)
# cowplot::plot_grid(
# 	cowplot::plot_grid(
# 		grid::grid.grabExpr(ComplexHeatmap::draw(hm, heatmap_legend_side = "left")), #
# 		# p$gtable,
# 		cowplot::plot_grid(
# 		get_SimiC_densisites('ZNF451') ,
# 		get_SimiC_densisites('YBX1')   ,
# 		get_SimiC_densisites('PSPC1')  ,
# 		get_SimiC_densisites('JARID2') ,
# 		get_SimiC_densisites('IRF1')   ,
# 		get_SimiC_densisites('KAT6B')  ,
# 		ncol=1, labels=c('B', '', '', 'C', '', ''),hjust = 0.02, vjust=0.9),
# 	ncol=2, labels=c('A', NULL), hjust = 0.02, vjust=0.9), 
# 	cowplot::plot_grid(
# 		get_SimiC_densisites('RERE')  + theme(axis.title.x=element_text(size=9)),
# 		get_SimiC_densisites('KDM2A') + theme(axis.title.x=element_text(size=9)), 
# 	ncol=2, labels=c('D', NULL), hjust = 0.02, vjust=0.5),
# 	simic_legend,
# nrow=3, rel_heights=c(0.85, 0.15, 0.05))
# dev.off()


pdf('./Plots/PaperFigures/FigS4.pdf', height=7.5)

cowplot::plot_grid(
	grid::grid.grabExpr(ComplexHeatmap::draw(hm, heatmap_legend_side = "left")), 
	dotplot_TFs, 
	ncol=2, labels = c('A', 'B'))
	
dev.off()