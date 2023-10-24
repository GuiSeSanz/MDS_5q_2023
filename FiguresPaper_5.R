
library(ggplot2)
library(ComplexHeatmap)
library(enrichR)




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
text_size = 2)
return(p)
}

get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)),
	order.by = "freq", empty.intersections = "on")
	# plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
	# 		nrow=2, align='hv', rel_heights = c(3,1),
	# 		rel_widths = c(2,3))
	return(plot)
}


get_HM <- function(data, title){
	pheatmap::pheatmap(as.data.frame.matrix(table(data$source, data$target)), 
	cluster_rows=FALSE, cluster_cols=FALSE, main = title, color=RColorBrewer::brewer.pal(9,"Blues"))
}


get_HM <- function(data, title){
	pheatmap::pheatmap(as.data.frame.matrix(table(data$source, data$target)), 
	cluster_rows=FALSE, cluster_cols=FALSE, main = title, color=RColorBrewer::brewer.pal(9,"Blues"))
}

data <- read.table('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/ResultsPerSample_CT_and_Condition.csv', sep='\t', header=TRUE, row.names=1)

data_flt <-  data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]

# cehck all common between them
data_elder <- data_flt[grepl('^GSM', data_flt$Sample),  ]
data_elder$Unique_Interaction <- paste(stringr::str_extract(data_elder$source, '[\\w]+(?=&)'), stringr::str_extract(data_elder$target, '[\\w]+(?=&)'), data_elder$ligand_complex, data_elder$receptor_complex, sep='_')

data_MDS <- data_flt[grepl('^SMD', data_flt$Sample),  ]
data_MDS$Unique_Interaction <- paste(stringr::str_extract(data_MDS$source, '[\\w]+(?=&)'), stringr::str_extract(data_MDS$target, '[\\w]+(?=&)'), data_MDS$ligand_complex, data_MDS$receptor_complex, sep='_')


data_MDS_flt   <- data_MDS[-log10(data_MDS$magnitude_rank) > 5, ]
data_elder_flt <- data_elder[-log10(data_elder$magnitude_rank) > 5, ]



ground_healthy <- data_elder_flt$Unique_Interaction
data_MDS_flt_NoT <- data_MDS_flt[!data_MDS_flt$target %in% c('T&normal', 'T&del5q') & !data_MDS_flt$source %in% c('T&normal', 'T&del5q'), ]
tmp <- data_MDS_flt_NoT[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,'Unique_Interaction'])
all_shared_comm <- Reduce(intersect,tmp)

comms_core_MDS <- data_MDS_flt_NoT[data_MDS_flt_NoT$Unique_Interaction %in% all_shared_comm , ]
comms_core_MDS_no_Healthy <- comms_core_MDS[!comms_core_MDS$Unique_Interaction %in% ground_healthy, ]

data_venn <- list('MDS unique interactions' = unique(comms_core_MDS$Unique_Interaction), 
				  'Healthy uniqueinteractions' = unique(data_elder$Unique_Interaction))

#Not DOne 
	# pdf('./Plots/PaperFigures/Figure6.pdf')
	# get_gvenn( data_venn)
	# dev.off()
	# all_sources <- sort(unique(stringr::str_remove(c(comms_core_MDS_no_Healthy$source, data_elder_flt$source), '&[\\w]+$')))
	# all_targets <- sort(unique(stringr::str_remove(c(comms_core_MDS_no_Healthy$target, data_elder_flt$target), '&[\\w]+$')))

	# Enrichment_results <- data.frame( Cluster=NULL, pval=NULL)
	# for (ct in unique(all_sources, all_targets)){
	# 	print(ct)
	# 	# q <- communications in MDS involving the ct
	# 	# m <- number of communications in MDS
	# 	# n <- number of communications in elder
	# 	# k <- communications in involving the ct
	# 	# phyper(q, m, n, k)
	# 	q <- unique(comms_core_MDS_no_Healthy[grepl(ct, comms_core_MDS_no_Healthy$source) | grepl(ct, comms_core_MDS_no_Healthy$target), 'Unique_Interaction'])
	# 	m <- unique(comms_core_MDS_no_Healthy$Unique_Interaction)
	# 	n <- unique(data_elder_flt$Unique_Interaction)
	# 	k <- unique(c(data_elder_flt[grepl(ct, data_elder_flt$source) | grepl(ct, data_elder_flt$target), 'Unique_Interaction'], 
	# 			comms_core_MDS_no_Healthy[grepl(ct, comms_core_MDS_no_Healthy$source) | grepl(ct, comms_core_MDS_no_Healthy$target), 'Unique_Interaction']))
	# 	pval_phyper <- phyper(length(q) -1 , length(m),length(n),length(k), lower.tail=FALSE)
	# 	Enrichment_results <- rbind(Enrichment_results, data.frame(Cluster=ct, pval=pval_phyper))
	# }
	# Enrichment_results$signif <- Enrichment_results$pval<0.05

	# Enrichment_results_2 <- data.frame( source=NULL, target=NULL, pval=NULL)
	# for (ct_source in all_sources){
	# 	print(paste0('FROM:  ',ct_source))
	# 	for (ct_target in all_targets){
	# 		print(paste0('TO:  ', ct_target))
	# 		# q <- communications in MDS involving the ct
	# 		# m <- number of communications in MDS
	# 		# n <- number of communications in elder
	# 		# k <- communications in involving the ct
	# 		# phyper(q, m, n, k)
	# 		q <- unique(comms_core_MDS_no_Healthy[grepl(ct_source, comms_core_MDS_no_Healthy$source) & grepl(ct_target, comms_core_MDS_no_Healthy$target), 'Unique_Interaction'])
	# 		m <- unique(comms_core_MDS_no_Healthy[grepl(ct_source, comms_core_MDS_no_Healthy$source), 'Unique_Interaction'])
	# 		n <- unique(data_elder_flt[grepl(ct_source, data_elder_flt$source), 'Unique_Interaction'])
	# 		k <- unique(c(data_elder_flt[grepl(ct_target, data_elder_flt$target), 'Unique_Interaction'], 
	# 				comms_core_MDS_no_Healthy[grepl(ct_target, comms_core_MDS_no_Healthy$target), 'Unique_Interaction']))
	# 		pval_phyper <- phyper(length(q) -1 , length(m),length(n),length(k), lower.tail=FALSE)
	# 		Enrichment_results_2 <- rbind(Enrichment_results_2, data.frame(source=ct_source, target=ct_target, pval=pval_phyper))
	# 	}
	# }
	# Enrichment_results_2$signif <- Enrichment_results_2$pval<0.05


setEnrichrSite("Enrichr")

get_enrichr_results <- function(genes){
	dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human', 'Reactome_2022')
	all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL)
		enriched <- enrichr(genes, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
		if(nrow(enriched) > 1){
			all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes))
		}
	return(all_results)
}

all_targets <- sort(unique(stringr::str_remove(c(comms_core_MDS_no_Healthy$target, data_elder_flt$target), '&[\\w]+$')))
all_sources <- sort(unique(stringr::str_remove(c(comms_core_MDS_no_Healthy$source, data_elder_flt$source), '&[\\w]+$')))
enrich_results <- data.frame(Comparison=NULL, Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, Comparison=NULL)
for (celltype in unique(all_sources, all_targets)){
	print(celltype)
	comms_core_MDS_no_Healthy_tmp <- comms_core_MDS_no_Healthy[grepl(celltype, comms_core_MDS_no_Healthy$source), ]
	data_elder_flt_tmp <- data_elder_flt[grepl(celltype, data_elder_flt$source), ]
	if(nrow(comms_core_MDS_no_Healthy_tmp) != 0){
		enrich_mds   <- get_enrichr_results(unique(comms_core_MDS_no_Healthy_tmp$ligand_complex, comms_core_MDS_no_Healthy_tmp$receptor_complex))
		enrich_mds$Comparison   <- paste0(celltype, '_MDS')
	}else{
		enrich_mds <- data.frame(Comparison=NULL, Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, Comparison=NULL)
	}
	if(nrow(data_elder_flt_tmp) != 0){
		enrich_elder <- get_enrichr_results(unique(data_elder_flt$ligand_complex, data_elder_flt$receptor_complex))
		enrich_elder$Comparison <- paste0(celltype, '_Elder')
	}else{
		enrich_elder <- data.frame(Comparison=NULL, Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, Comparison=NULL)
	}
	enrich_tmp <- rbind(enrich_mds, enrich_elder)
	enrich_results <- rbind(enrich_results, enrich_tmp)
}


plotter <- enrich_results
plotter[plotter$Combined.Score >210, 'Combined.Score'] <- 210
terms_2_plot <- c(
'DAP12 Signaling', 'HIF\\-1 signaling pathway', 'Degradatiom of extracellular Matrix', 'Activation of Matrix metalloproteinases', 'Extracellular Matrix organization', 'Hemostasis', 'microRNAs in cancer', 'N\\-glycan Trimming in ER and ', 'Oncogenic MAPK Signaling', 'Apoptosis', 'Apoptotic Cleveage of cellular proteins', 'Progammed cell death', 'Apoptotic process', 'RHOBTB1 GTPase cycle', 'RHOBTB GTPase cycle', 'RHO GTPase cycle', '^Antigen', 'hemopoiesis', 'hematopoietic stem cell proliferation', 'negative regulation of cell development', 'GO:1902254', 'GO:1902230', 'GO:1902166', 'GO:0017148', 'GO:2000774', 'GO:2001200', 'GO:1902107', 'GO:0045931')
terms_2_plot <- unique(grep(paste(terms_2_plot, collapse='|'), plotter$Term, value=T))
plotter <- plotter[plotter$Term %in% terms_2_plot, ]
plotter$Term <- stringr::str_remove(plotter$Term,'\\([^)]*\\)')
plotter$Term <- stringr::str_remove(plotter$Term, '[\\s]+R-HSA-[\\d]+')

plotter$Celltype   <- stringr::str_extract(plotter$Comparison, '^[\\w-]+(?=_)')
plotter$Comparison <- stringr::str_extract( plotter$Comparison, '(?<=_)[\\w][^_]+$')

plot_terms <- c(
'apoptotic signaling pathway in response', 
'apoptotic signaling pathway by p53 class', 
'stem cell proliferation', 
'hemopoiesis', 
'leukocyte differentiation', 
'dendritic cell differentiation', 
'DAP12', 
'cellular senescence', 
'regulation of translation', 
'MAPK', 
'HIF-1 signaling pathway')

plot_terms_order <- c("negative regulation of intrinsic apoptotic signaling pathway in response to DNA damage ", 
"negative regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator ", 
"negative regulation of intrinsic apoptotic signaling pathway by p53 class mediator ", 
"hematopoietic stem cell proliferation ", 
"hemopoiesis ", 
"positive regulation of leukocyte differentiation ", 
"positive regulation of dendritic cell differentiation ", 
"DAP12 Signaling", 
"positive regulation of cellular senescence ", 
"negative regulation of translation ", 
"Oncogenic MAPK Signaling", 
"HIF-1 signaling pathway")


plotter <- plotter[plotter$Term %in% plot_terms_order, ]
plotter$Term <- factor(plotter$Term, levels=rev(plot_terms_order))

get_dotplot <- function(data){
ggplot() + 
geom_point(data, mapping= aes(y=Term, x=Celltype, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(name='Combined\nScore', low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 45)) +
facet_wrap(~Comparison, ncol=2)+
scale_alpha(range = c(1, 0.1)) + 
scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(5, 1), breaks=c(0.01, 0.05), limits=c(0, 0.1)) +
theme_minimal() +   
theme(legend.position='right',
	axis.text.y=element_text(size=6, color='black'),
	axis.text.x=element_text(size = 6, angle=45, vjust = 1, hjust = 1, color='black'),
	legend.key.size = unit(0.2, "cm"), 
	legend.title=element_text(size=7), 
	axis.title=element_blank(), 
	legend.text=element_text(size=5)) 
}

# pdf('./Plots/Liana/MDS_and_Elder/Fig6_DP.pdf')
# get_dotplot(plotter)
# dev.off()

# WriteXLS::WriteXLS(enrich_results, './Plots/Liana/MDS_and_Elder/Enrichr_ALL.xlsx')

# generate_evenly_distributed_vector <- function(start, end, num_values) {
#   if (num_values <= 1) {
#     stop("Number of values must be greater than 1.")
#   }
#   sequence_length <- num_values - 1
#   step_size <- (end - start) / sequence_length
#   vector <- seq(start, end, by = step_size)
#   return(vector)
# }



# generate_evenly_distributed_vector(0, max(max(hm_MDS), max(hm_Elder)), 9)
# color_func_HM <- circlize::colorRamp2(generate_evenly_distributed_vector(0, max(max(hm_MDS), max(hm_Elder)), 9),RColorBrewer::brewer.pal(9,"Blues"))
celltypes_mds_cols <- c(
	'HSC del5q',
	'HSC normal', 
	'LMPP del5q', 
	'LMPP normal',
	'GMP del5q',
	'GMP normal', 
	'Granulocyte del5q',
	'Granulocyte normal',
	'Monocytes del5q',
	'Monocytes normal', 
	'DendriticCell del5q', 
	'DendriticCell normal',
	'CLP del5q',
	'CLP normal', 
	'pro-B del5q',
	'pro-B normal', 
	'MEP del5q',
	'MEP normal', 
	'MK_Prog del5q', 
	'MK_Prog normal', 
	'EarlyErythroid del5q',
	'EarlyErythroid normal',
	'Basophil del5q',
	'Basophil normal')

celltypes_mds_rows <- c(
	'HSC del5q',
	'HSC normal', 
	'LMPP del5q', 
	'LMPP normal',
	'GMP del5q',
	'GMP normal', 
	'Granulocyte del5q',
	'Granulocyte normal',
	'Monocytes del5q',
	'Monocytes normal', 
	'DendriticCell del5q', 
	'DendriticCell normal',
	'CLP del5q',
	'CLP normal', 
	'pro-B del5q',
	'pro-B normal', 
	'MEP del5q',
	'MEP normal', 
	'MK_Prog del5q', 
	'MK_Prog normal', 
	'LateErythroid del5q',
	'LateErythroid normal',
	'Basophil del5q',
	'Basophil normal')
hm_MDS <-  as.matrix(table(comms_core_MDS_no_Healthy$source, comms_core_MDS_no_Healthy$target))
colnames(hm_MDS) <- sub('&', ' ', colnames(hm_MDS))
rownames(hm_MDS) <- sub('&', ' ', rownames(hm_MDS))
hm_MDS <- hm_MDS[rownames(hm_MDS) %in% celltypes_mds_rows, colnames(hm_MDS) %in% celltypes_mds_cols]
hm_MDS <- as.data.frame.matrix(hm_MDS)
hm_MDS <- hm_MDS[celltypes_mds_rows,  celltypes_mds_cols]
hm_MDS <- as.matrix(hm_MDS)

row_ha = rowAnnotation(Total_row = anno_barplot(rowSums(hm_MDS), width = unit(0.4, "cm")) , show_annotation_name=FALSE)
column_ha = HeatmapAnnotation(Total_column = anno_barplot(colSums(hm_MDS), height = unit(0.4, "cm")), show_annotation_name=FALSE)

HM_MDS <- Heatmap(hm_MDS, name='Number of\ninteractions', 
col = RColorBrewer::brewer.pal(9,"Blues"), 
rect_gp = gpar(col = "grey", lwd = 1),
top_annotation = column_ha, right_annotation = row_ha, 
cluster_rows = FALSE,
cluster_columns = FALSE, 
row_title = 'Source',
column_title = 'Target',column_names_rot=45,
column_names_gp = grid::gpar(fontsize = 6, angle=45),
row_names_gp = grid::gpar(fontsize = 4.5),
show_heatmap_legend = TRUE
)

cell_type_order <- 	c('HSC', 
	'LMPP', 
	'GMP', 
	'Granulocyte',
	'Monocytes', 
	'DendriticCell', 
	'CLP', 
	'pro-B', 
	'MEP', 
	'MK_Prog', 
	'EarlyErythroid',
	'LateErythroid', 
	'Basophil')


hm_Elder <-  as.matrix(table(data_elder_flt$source, data_elder_flt$target))
colnames(hm_Elder) <- sub('&WT', '', colnames(hm_Elder))
rownames(hm_Elder) <- sub('&WT', '', rownames(hm_Elder))
hm_Elder <- as.data.frame.matrix(hm_Elder)
hm_Elder <- hm_Elder[cell_type_order,  cell_type_order]
hm_Elder <- as.matrix(hm_Elder)
row_ha = rowAnnotation(Total_row = anno_barplot(rowSums(hm_Elder), width = unit(0.4, "cm")), show_annotation_name=FALSE)
column_ha = HeatmapAnnotation(Total_column = anno_barplot(colSums(hm_Elder),  height = unit(0.4, "cm")), show_annotation_name=FALSE)
HM_ELDER <- Heatmap(hm_Elder, name='Number of\ninteractions',
col = RColorBrewer::brewer.pal(9,"Greens"), 
rect_gp = gpar(col = "grey", lwd = 1),
top_annotation = column_ha, right_annotation = row_ha, 
cluster_rows = FALSE,
cluster_columns = FALSE, 
row_title = 'Source',
column_title = 'Target', column_names_rot=45,
column_names_gp = grid::gpar(fontsize = 6, angle=45),
row_names_gp = grid::gpar(fontsize = 4.5),
heatmap_legend_param = list(direction = "vertical"),
show_heatmap_legend=FALSE
) 


pdf('./Plots/PaperFigures/Fig5.pdf', 10, 12)
cowplot::plot_grid(
	cowplot::plot_grid(
		grid::grid.grabExpr(ComplexHeatmap::draw(HM_ELDER, heatmap_legend_side = "left")),
		grid::grid.grabExpr(ComplexHeatmap::draw(HM_MDS, heatmap_legend_side = "right")), 
	ncol=2, rel_widths=c(1,1)), 
	cowplot::plot_grid(get_dotplot(plotter),NULL, ncol=2, rel_widths=c(1,0.3)),
ncol=1, rel_heights=c(0.6,1))
get_gvenn( data_venn)
get_dotplot(plotter)
dev.off()






interactions_elder <- c('B_HSC_HMGB1_CXCR4', 'Basophil_HSC_HMGB1_CXCR4', 'Basophil_HSC_HMGB1_CXCR4', 'Basophil_HSC_HMGB1_CXCR4', 'CLP_HSC_HMGB1_CXCR4', 'CLP_HSC_HMGB1_CXCR4', 'CLP_HSC_HMGB1_CXCR4', 'DendriticCell_HSC_HMGB1_CXCR4', 'DendriticCell_HSC_HMGB1_CXCR4', 'DendriticCell_HSC_HMGB1_CXCR4', 'EarlyErythroid_HSC_HMGB1_CXCR4', 'EarlyErythroid_HSC_HMGB1_CXCR4', 'EarlyErythroid_HSC_HMGB1_CXCR4', 'Granulocyte_HSC_HMGB1_CXCR4', 'Granulocyte_HSC_HMGB1_CXCR4', 'Granulocyte_HSC_HMGB1_CXCR4', 'LateErythroid_HSC_HMGB1_CXCR4', 'LateErythroid_HSC_HMGB1_CXCR4', 'LateErythroid_HSC_HMGB1_CXCR4', 'MK_Prog_HSC_HMGB1_CXCR4', 'MK_Prog_HSC_HMGB1_CXCR4', 'MK_Prog_HSC_HMGB1_CXCR4')
interactions_mds <- c(
'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Granulocyte_Monocytes_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_Basophil_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_CLP_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_EarlyErythroid_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1', 'Monocytes_GMP_AGTRAP_RACK1')

library(liana)

pdf('chords.pdf', 10, 10)
chord_freq(data_elder_flt[data_elder_flt$Unique_Interaction %in% interactions_elder,], cex =0.4)
chord_freq(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Unique_Interaction %in% interactions_mds,], cex =0.4)
cowplot::plot_grid(cowplot::get_legend(chord_freq(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Unique_Interaction %in% interactions_mds,], cex =0.4)))
dev.off()