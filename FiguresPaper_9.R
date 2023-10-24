
library(enrichR)
library(ggplot2)

get_enrichr_results <- function(data_set){
	all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
	for (CT in names(data_set)){
		message(CT)
		tmp <- data_set[[CT]]
		if (!'gene' %in% colnames(tmp)){
			tmp$gene <- rownames(tmp)
		}
		if(!'avg_logFC' %in% colnames(tmp)){
			tmp$avg_logFC <- tmp$avg_log2FC
		}
		# tmp_pos <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >0,'gene']
		# tmp_neg <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC <0,'gene']
		tmp_pos <- tmp[tmp$avg_logFC >0,'gene']
		tmp_neg <- tmp[tmp$avg_logFC <0,'gene']
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

get_dotplot <- function(data){
ggplot() + 
geom_point(data, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(name = 'Combined\nScore', low =scales::muted("blue"),mid = "white",high = '#db3e25',midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(3, 0.5)) +
facet_wrap(~Comparison, nrow=1) + 
theme(axis.text.y=element_text(size=8),
axis.text.x=element_text(size = 8, angle=45, vjust = 1, hjust = 1),
axis.title.x=element_blank(),
legend.position='right') + theme(legend.key.size = unit(0.2, "cm"), 
legend.title=element_text(size=7), legend.text=element_text(size=5)) 
}


get_SimiC_densisites <- function( data_auc=df_auc, regulon){
	phenotypes <- c('Diagnosed MDS', 'Non Responder', 'Partial Responder')
	color_map <- c("#DF8444","#023047","#499957")
	names(color_map) <- phenotypes
	p <- ggplot(data_auc[data_auc$driver ==regulon,], aes(x=value, y=..scaled.., fill=TimePoint)) + #xlim(0,1) + 
	geom_density(alpha = 0.6, adjust = 1/2) + theme_classic() + 
	scale_fill_manual(values=color_map) + scale_y_continuous(breaks=seq(0,1,0.5)) +
	theme(legend.position = 'none') + 
	labs(y=regulon, x='Activity score') + 
	theme(axis.title.x=element_blank(), 
	axis.title.y=element_text(size=8))
	return(p)
}

# get_SimiC_regulon <- function(data, regulon){
# 	regulon_densities <- list()
# 	phenotypes <- c('Elder', 'Diagnosed MDS', 'Non Responder', 'Partial Responder')
# 	color_map <- c("#f94144","#219ebc","#023047","#fb8500")
# 	p <- ggplot(data[data$driver ==regulon,], aes(x=value, fill=TimePoint)) + 
# 				geom_density(alpha = 0.8, adjust = 1/8) + theme_classic() + 
# 				scale_fill_manual(values=color_map) +
# 				theme(legend.position = 'top') + 
# 				ggtitle(paste0(regulon))
# 	return(p)
# }

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human', 'Reactome_2022')


Results_DE_Del5qPRVsDel5qNR_MAST <- readRDS('./Data/Results_DE_Del5qPRVsDel5qNR_MAST.rds')
enrich_results <- get_enrichr_results(Results_DE_Del5qPRVsDel5qNR_MAST)

terms_2_plot <- c('GO:0031624', 'GO:0044390', 'GO:0004842', 'GO:0010498', 'GO:0043161', 'GO:1905037', 'GO:0000045', 'Autophagy', 'Ubiquitin mediated proteolysis', 'Phosphatidylinositol signaling system', 'R-HSA-9006335', 'R-HSA-9027276', 'PD\\-L1')
terms_2_plot <- unique(grep(paste(terms_2_plot, collapse='|'), enrich_results$Term, value=T))

enrich_results[enrich_results$Combined.Score >200, 'Combined.Score'] <- 200
enrich_results$Combined.Score_size <- enrich_results$Combined.Score

terms_2_remove <- grep(paste(c('9612973', '9613829', '9663891'),collapse='|'), terms_2_plot, value=T)
terms_2_plot <- terms_2_plot[!terms_2_plot %in% terms_2_remove]

enrich_results <- enrich_results[enrich_results$Term %in% terms_2_plot, ]

enrich_results$Term <- stringr::str_remove(enrich_results$Term,'\\([^)]*\\)')
enrich_results[enrich_results$is_pos == FALSE, 'Combined.Score'] <- enrich_results[enrich_results$is_pos == FALSE, 'Combined.Score']*-1
# get_order
for (term in unique(enrich_results$Term)){
	enrich_results[enrich_results$Term == term, 'Combined.Score_OVERALL'] <- sum(enrich_results[enrich_results$Term == term, 'Combined.Score'])
}
sort_order <-  unique(enrich_results[,c('Term', 'Combined.Score_OVERALL') ])
sort_order <- sort_order[order(sort_order$Combined.Score_OVERALL, decreasing=FALSE),'Term']
enrich_results$Term <- factor(enrich_results$Term, levels = sort_order)

enrich_results$cell_type <- factor(enrich_results$cell_type, levels=sort(unique(enrich_results$cell_type)))
enrich_results$Comparison <- 'Del(5q) Partial Responder vs Del(5q) Non Responder'
enrich_results <- enrich_results[enrich_results$Combined.Score>0, ]
enrich_results$Term <- stringr::str_remove(enrich_results$Term,'\\([^)]*\\)')
enrich_results$Term <- stringr::str_remove(enrich_results$Term,'(?=R-HSA)[\\w\\d-]+')


SimiC_aucs_PostDel3 <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_df_auc_PostDel3.rds')
SimiC_aucs_PostDel3$TimePoint <- ifelse(SimiC_aucs_PostDel3$TimePoint == 'MDS', 'Diagnosed MDS', 
								 ifelse(SimiC_aucs_PostDel3$TimePoint == 'NR_post', 'Non Responder', 
								 ifelse(SimiC_aucs_PostDel3$TimePoint == 'PR', 'Partial Responder', 'UPS')))

pdf('./Plots/PaperFigures/Fig9.pdf', width=8.3, height=8)
legend <- cowplot::get_legend(get_SimiC_densisites(SimiC_aucs_PostDel3, 'IRF1') + theme(legend.position='bottom'))
cowplot::plot_grid(
	get_dotplot( enrich_results),
	cowplot::plot_grid(
		cowplot::plot_grid( get_SimiC_densisites(SimiC_aucs_PostDel3, 'IRF1') + theme(legend.position='none'),
							get_SimiC_densisites(SimiC_aucs_PostDel3, 'JARID2') + theme(legend.position='none'),
							get_SimiC_densisites(SimiC_aucs_PostDel3, 'NCOR1') + theme(legend.position='none')+ theme(axis.title.x=element_text(size=9)),
							get_SimiC_densisites(SimiC_aucs_PostDel3, 'CUX1') + theme(legend.position='none')+ theme(axis.title.x=element_text(size=9)),
							ncol=2, nrow=2),
		legend, 
		nrow=2, rel_heights = c(0.9, 0.1)),
	nrow=2, rel_heights = c(1, 1)
	)
dev.off()
