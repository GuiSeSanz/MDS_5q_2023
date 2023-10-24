
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
		if(nrow(enriched) > 0){
			all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
		}
		enriched <- enrichr(tmp_neg, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
		if(nrow(enriched) > 0){
			all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))
		}
	}
	return(all_results)
}


get_dotplot <- function(data, Comparison){
	ggplot() + 
	geom_point(data, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
	scale_fill_gradient2(name = 'Combined\nScore', low ="#4882cf",mid = "white",high = "#db3e25",midpoint = 0)+
	scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 50)) +
	scale_alpha(range = c(1, 0.1)) + 
	scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(5, 0.5), breaks=c(0.01, 0.05, 0.1, 0.5), limits=c(0, 0.5)) +
	facet_wrap(~Comparison, nrow=1) + 
	theme_minimal() +   
	theme(legend.position='bottom',
		axis.text.y=element_text(size=6),
		axis.text.x=element_text(size = 6, angle=45, vjust = 1, hjust = 1),
		legend.key.size = unit(0.2, "cm"), 
		legend.title=element_text(size=7), 
		axis.title=element_blank(), 
		legend.text=element_text(size=5)) 
}


setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Biological_Process_2021")



get_dotplot_filteredOrdered_GO <- function(enrichr_results, GO_list, description){
	terms_2_plot <- unique(grep(paste(GO_list, collapse='|'), enrichr_results$Term, value=T))

	enrichr_results$GO_Term <- stringr::str_extract(enrichr_results$Term,'\\(GO:[^)]*\\)')
	enrichr_results$GO_Term <- gsub('\\(|\\)', '', enrichr_results$GO_Term)

	enrichr_results <- enrichr_results[enrichr_results$Term %in% terms_2_plot, ]
	print(unique(enrichr_results$cell_type))

	# enrichr_results <- enrichr_results[match(GO_list, enrichr_results$GO_Term),]

	enrichr_results$Term <- factor(enrichr_results$Term, levels = unique(enrichr_results$Term))
	enrichr_results[enrichr_results$is_pos == FALSE, 'Combined.Score'] <- enrichr_results[enrichr_results$is_pos == FALSE, 'Combined.Score']*-1
	# saturate the color scale
	enrichr_results[enrichr_results$Combined.Score >  200, 'Combined.Score'] <- 200
	enrichr_results[enrichr_results$Combined.Score < -200, 'Combined.Score'] <- -200
	# remove with regex the GO terms
	enrichr_results$Term <- stringr::str_remove(enrichr_results$Term,'\\s\\(GO[^)]*\\)')
	print(unique(enrichr_results$cell_type))
	enrichr_results$Comparison <- description
	p <- get_dotplot(enrichr_results)
	return(p)
}


Results_DE_del5qPR_Vs_del5qMDS    <- readRDS('/mnt/md0/gserranos/MDS/Data/Results_DE_del5qPR_Vs_del5qMDS.rds')
Results_DE_del5qPR_Vs_del5qMDS <- get_enrichr_results(Results_DE_del5qPR_Vs_del5qMDS)


Results_DE_del5qNR_Post_Vs_del5qMDS <- readRDS('/mnt/md0/gserranos/MDS/Data/Results_DE_del5qNR_Post_Vs_del5qMDS.rds')
Results_DE_del5qNR_Post_Vs_del5qMDS <- get_enrichr_results(Results_DE_del5qNR_Post_Vs_del5qMDS)








GO_list_A <- c('GO:0004842', 'GO:0061659', 'GO:0061630', 'GO:0061631', 'GO:0000209', 'GO:0016567', 'GO:0070936', 'GO:0031624', 'GO:0043130', 'GO:0044390', 'GO:0043065', 'GO:0042981', 'GO:0006412', 'GO:0003743', 'GO:0006415', 'GO:0070126', 'GO:0006414', 'GO:0032543', 'GO:0006296', 'GO:0033683', 'GO:0006283', 'GO:0006289', 'GO:0006287', 'GO:0006295', 'GO:0006293', 'GO:0006284', 'GO:0006281', 'GO:0006302', 'GO:0000725', 'GO:0000045', 'GO:1905037', 'GO:0016236', 'GO:0090305', 'GO:0000422', 'GO:0006661', 'GO:0046488', 'GO:0046854', 'GO:0046856', 'GO:0004438', 'GO:0052744', 'GO:0005545', 'GO:0106018', 'GO:0016307')
GO_list_B <- c('GO:0097027', 'GO:1904668', 'GO:1904666', 'GO:0051443', 'GO:0043161', 'GO:0006511', 'GO:0032436', 'GO:1902850', 'GO:0000070', 'GO:0007080', 'GO:1901990', 'GO:0140014', 'GO:0044772', 'GO:0010389', 'GO:0090068', 'GO:0010564', 'GO:0007346', 'GO:0040001', 'GO:0051315', 'GO:0000819', 'GO:0000086', 'GO:0044839', 'GO:0051781', 'GO:1901992', 'GO:0010971', 'GO:1902751', 'GO:0045840')

pdf('./Plots/PaperFigures/FigS5.pdf')
cowplot::plot_grid(
	get_dotplot_filteredOrdered_GO(Results_DE_del5qPR_Vs_del5qMDS, GO_list_A, 'del(5q) PR Vs del(5q) MDS'),
	get_dotplot_filteredOrdered_GO(Results_DE_del5qNR_Post_Vs_del5qMDS, GO_list_B, 'del(5q) NR Vs del(5q) MDS'),
ncol=2, labels = c('A', 'B'))

dev.off()


# table_suppl_Fig_5 




