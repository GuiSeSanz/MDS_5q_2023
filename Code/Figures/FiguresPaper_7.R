
library(enrichR)
library(ggplot2)

get_enrichr_results <- function(data_set, DIR=NULL){
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
		enriched <- enriched[enriched$Adjusted.P.value < 0.5,]
		if(nrow(enriched) > 1){
			if(DIR == 'POS'){
				all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
			}
		}
		enriched <- enrichr(tmp_neg, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.5,]
		if(nrow(enriched) > 1){
			if(DIR == 'NEG'){
				all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))
			}
		}
	}
	return(all_results)
}

get_SimiC_regulon <- function(data, regulon){
	regulon_densities <- list()
	data$TimePoint <- ifelse(data$TimePoint == 'Elder', 'Healthy', 
					  ifelse(data$TimePoint == 'Diagnosed', 'Diagnosed MDS', 
					  ifelse(data$TimePoint == 'PR', 'Partial Responder', 
					  ifelse(data$TimePoint == 'CR', 'Complete Responder', 'UPS'))))
	phenotypes <- c('Healthy', 'Diagnosed MDS', 'Partial Responder', 'Complete Responder')
	color_map <- c("#3172AC","#DF8444","#499957","#C93839")
	names(color_map) <- phenotypes
	p <- ggplot(data[data$driver ==regulon,], aes(x=value, fill=TimePoint)) + 
				geom_density(alpha = 0.6, adjust = 1/8) + 
				ggtitle(paste0(regulon)) +
				theme_classic() + 
				scale_fill_manual(values=color_map) +
				theme(
					axis.title=element_blank(),
					axis.text=element_text(size=6),
					legend.position='none',
					plot.title = element_text(size = 8))
	return(p)
}

get_dotplot <- function(data){
ggplot() + 
geom_point(data, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
scale_fill_gradient2(name = 'Combined\nScore', low ="#4882cf",mid = "white",high = "#db3e25",midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
scale_alpha(range = c(1, 0.1)) + 
scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(5, 1), breaks=c(0.01, 0.05, 0.1, 0.5), limits=c(0, 0.5)) +
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

######################################################################################################
setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human', 'Reactome_2022')

MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')

post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
all_5q_depleted_cells_COPYKAT <- readRDS('./Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')

######################## DOTPLOTS VS MDS ########################
all_DE_pseudo_Non5qCR_Vs_non5qDiagnosed   <- readRDS('./Data/Results_DE_Non5qCR_Vs_non5qDiagnosed.rds')
all_results_Non5qCR_Vs_non5qDiagnosed <- get_enrichr_results(all_DE_pseudo_Non5qCR_Vs_non5qDiagnosed, 'POS')

all_DE_pseudo_Non5qPR_Vs_non5qDiagnosed <- readRDS('./Data/Results_DE_Non5qPR_Vs_non5qDiagnosed.rds')
all_results_Non5qPR_Vs_non5qDiagnosed <- get_enrichr_results(all_DE_pseudo_Non5qPR_Vs_non5qDiagnosed, 'POS')

#  Non5qCR_Vs_non5qDiagnosed -> Complete Responder Vs Diagnosed
#  Non5qPR_Vs_non5qDiagnosed -> Partial Responder Vs Diagnosed

all_results_Non5qCR_Vs_non5qDiagnosed$Comparison <- 'Complete Responder\nVs\nDiagnosed'
all_results_Non5qPR_Vs_non5qDiagnosed$Comparison <- 'Partial Responder\nVs\nDiagnosed'
all_results_Vs_MDS <- rbind(all_results_Non5qCR_Vs_non5qDiagnosed, all_results_Non5qPR_Vs_non5qDiagnosed)

all_results_Vs_MDS[all_results_Vs_MDS$Combined.Score >210, 'Combined.Score'] <- 210
all_results_Vs_MDS$Combined.Score_size <- all_results_Vs_MDS$Combined.Score

terms_2_plot <- c('GO:0031624', 'GO:0044390', 'GO:0004842', 'GO:0010498', 'GO:0043161', 'GO:1905037', 'GO:0000045', 'Autophagy', 'Ubiquitin mediated proteolysis', 'Phosphatidylinositol signaling system', 'R-HSA-9006335', 'R-HSA-9027276', 'PD\\-L1')
terms_2_plot <- unique(grep(paste(terms_2_plot, collapse='|'), all_results_Vs_MDS$Term, value=T))

terms_2_remove <- grep(paste(c('9612973', '9663891'), collapse='|'), terms_2_plot, value=T)
terms_2_plot <- terms_2_plot[!terms_2_plot %in% terms_2_remove]


all_results_Vs_MDS <- all_results_Vs_MDS[all_results_Vs_MDS$Term %in% terms_2_plot, ]
all_results_Vs_MDS$Term <- stringr::str_remove(all_results_Vs_MDS$Term,'\\([^)]*\\)')
all_results_Vs_MDS[all_results_Vs_MDS$is_pos == FALSE, 'Combined.Score'] <- all_results_Vs_MDS[all_results_Vs_MDS$is_pos == FALSE, 'Combined.Score']*-1
# get_order
for (term in unique(all_results_Vs_MDS$Term)){
	all_results_Vs_MDS[all_results_Vs_MDS$Term == term, 'Combined.Score_OVERALL'] <- sum(all_results_Vs_MDS[all_results_Vs_MDS$Term == term, 'Combined.Score'])
}
sort_order <-  unique(all_results_Vs_MDS[,c('Term', 'Combined.Score_OVERALL') ])
sort_order <- sort_order[order(sort_order$Combined.Score_OVERALL, decreasing=FALSE),'Term']
all_results_Vs_MDS$Term <- factor(all_results_Vs_MDS$Term, levels = sort_order)

all_results_Vs_MDS$cell_type <- factor(all_results_Vs_MDS$cell_type, levels=sort(unique(all_results_Vs_MDS$cell_type)))

# dotplot_nondelVsMDS <- ggplot() + 
# geom_point(all_results_Vs_MDS, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
# scale_fill_gradient2(low =scales::muted("blue"),mid = "white",high = scales::muted("red"),midpoint = 0, breaks=c(-200, -100, 0, 100, 200))+
# scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
# scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) +
# facet_wrap(~Comparison) + 
# theme(axis.text.y=element_text(size=8),
# axis.text.x=element_text(size = 8, angle=-90, vjust=0.2, hjust=0.05),
# legend.position='right') 

######################## DOTPLOTS VS ELDER ########################

all_DE_pseudo_Non5qCR_Vs_Elder   <- readRDS('./Data/Results_DE_Non5qCR_Vs_Elder.rds')
all_results_Non5qCR_Vs_Elder <- get_enrichr_results(all_DE_pseudo_Non5qCR_Vs_Elder,'NEG')

all_DE_pseudo_Non5qPR_Vs_Elder <- readRDS('./Data/Results_DE_Non5qPR_VsElder.rds')
all_results_Non5qPR_Vs_Elder <- get_enrichr_results(all_DE_pseudo_Non5qPR_Vs_Elder,'NEG')

all_results_Non5qCR_Vs_Elder$Comparison <- 'Complete Responder\nVs\nHealthy'
all_results_Non5qPR_Vs_Elder$Comparison <- 'Partial Responder\nVs\nHealthy'
all_results_Vs_Elder <- rbind(all_results_Non5qCR_Vs_Elder, all_results_Non5qPR_Vs_Elder)

all_results_Vs_Elder[all_results_Vs_Elder$Combined.Score >210, 'Combined.Score'] <- 210
all_results_Vs_Elder$Combined.Score_size <- all_results_Vs_Elder$Combined.Score

terms_2_plot <- c('GO:0006412','GO:0006364','GO:0043043','GO:0006415','GO:0070125','GO:0006414','GO:0032543','GO:0070126')
terms_2_plot <- unique(grep(paste(terms_2_plot, collapse='|'), all_results_Vs_Elder$Term, value=T))


all_results_Vs_Elder <- all_results_Vs_Elder[all_results_Vs_Elder$Term %in% terms_2_plot, ]

all_results_Vs_Elder$Term <- stringr::str_remove(all_results_Vs_Elder$Term,'\\([^)]*\\)')
all_results_Vs_Elder[all_results_Vs_Elder$is_pos == FALSE, 'Combined.Score'] <- all_results_Vs_Elder[all_results_Vs_Elder$is_pos == FALSE, 'Combined.Score']*-1
# get_order
for (term in unique(all_results_Vs_Elder$Term)){
	all_results_Vs_Elder[all_results_Vs_Elder$Term == term, 'Combined.Score_OVERALL'] <- sum(all_results_Vs_Elder[all_results_Vs_Elder$Term == term, 'Combined.Score'])
}
sort_order <-  unique(all_results_Vs_Elder[,c('Term', 'Combined.Score_OVERALL') ])
sort_order <- sort_order[order(sort_order$Combined.Score_OVERALL, decreasing=FALSE),'Term']
all_results_Vs_Elder$Term <- factor(all_results_Vs_Elder$Term, levels = sort_order)

all_results_Vs_Elder$cell_type <- factor(all_results_Vs_Elder$cell_type, levels=sort(unique(all_results_Vs_Elder$cell_type)))

# dotplot_nondelVsElder <- ggplot() + 
# geom_point(all_results_Vs_Elder, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
# scale_fill_gradient2(low =scales::muted("blue"),mid = "white",high = scales::muted("red"),midpoint = 0,)+
# scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
# scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(3, 0.5)) +
# facet_wrap(~Comparison) + 
# theme(axis.text.y=element_text(size=8),
# axis.text.x=element_text(size = 8, angle=-90, vjust=0.2, hjust=0.05),
# legend.position='right') 



all_results_dotplot <- rbind(all_results_Vs_Elder, all_results_Vs_MDS)
all_results_dotplot$Term <- factor(stringr::str_remove(all_results_dotplot$Term, '[\\s]+R-HSA-[\\d]+'), 
					levels= stringr::str_remove(levels(all_results_dotplot$Term), '[\\s]+R-HSA-[\\d]+'))

all_results_dotplot$Comparison <- factor(all_results_dotplot$Comparison, 
				levels = c('Complete Responder\nVs\nDiagnosed', 'Partial Responder\nVs\nDiagnosed', 'Complete Responder\nVs\nHealthy', 'Partial Responder\nVs\nHealthy'))
##############################################################################

#  There are less cells as some clusters have disapeared based on the annotation

features <- c('IKZF1', 'HBA1', 'CA2',  'HBB')
features <- intersect(intersect(features, rownames(post_data)), rownames(MDS_data))
tmp_post <- FetchData(post_data, vars=c('Sample', 'cell_5q', 'Cluster_names', features))
tmp_mds  <- FetchData(MDS_data, vars=c('Sample', 'cell_5q', 'Cluster_names', features))
plotter <- rbind(tmp_post, tmp_mds)
plotter <- reshape2::melt(plotter)
plotter$contrast <- ifelse(plotter$Sample %in% unique(MDS_data$Sample), 'Diagnosed', 
					ifelse(plotter$Sample == 'FS-0406-post', 'PR', 'CR'))
plotter <- plotter[plotter$Cluster_names %in% c('LateErythroid','MEP', 'EarlyErythroid'), ]
plotter$Cluster_names <- ifelse(plotter$Cluster_names == 'LateErythroid', 'Late\nErythroid', 
						 ifelse(plotter$Cluster_names == 'EarlyErythroid', 'Early\nErythroid', 'MEP'))

phenotypes <- c('Elder', 'Diagnosed\nMDS', 'Partial\nResponder', 'Complete\nResponder')
color_map <- c("#3172AC","#DF8444","#499957","#C93839")
names(color_map) <- phenotypes

plotter$contrast <- ifelse(plotter$contrast == 'Diagnosed', 'Diagnosed\nMDS', 
					ifelse(plotter$contrast == 'PR', 'Partial\nResponder', 'Complete\nResponder'))
boxplot_ery <- ggplot(plotter[plotter$cell_5q == 'normal', ], aes(x=Cluster_names, y=value, fill=contrast, color =contrast)) + 
geom_boxplot(outlier.size=0.4,outlier.alpha = 0.2)+ facet_wrap(~variable, nrow=2, ncol=2, scales='free_y') + 
theme_bw() + ylab('Expression') + 
scale_fill_manual(values=color_map[names(color_map) %in% unique(plotter$contrast)])+ 
scale_color_manual(values=scales::muted(color_map[names(color_map) %in% unique(plotter$contrast)]))+
theme(legend.position='right', axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x=element_blank(), 
strip.background = element_blank()) + 
guides(fill=guide_legend(title="Genotype"), color=guide_legend(title="Genotype"))
# ggpubr::stat_compare_means(aes(group=contrast), label = "p.signif", method="wilcox.test", hide.ns=FALSE)

SimiC_aucs_PostNormal <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_df_auc_PostNormal.rds')

pdf('./Plots/PaperFigures/Fig7.pdf', width=8.3, height=11.7)
legend <- cowplot::get_legend(get_SimiC_regulon(SimiC_aucs_PostNormal, 'IRF1') + 
					theme(legend.position='top', legend.key.size = unit(0.5, 'cm'))+ 
					guides(fill=guide_legend(title="Genotype")))
cowplot::plot_grid(
	cowplot::plot_grid(get_dotplot( all_results_dotplot),
	cowplot::plot_grid(
		cowplot::plot_grid(
					boxplot_ery,
				cowplot::plot_grid(
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'IRF1') ,
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'KAT6B'),
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'CUX1')  + ylab('Activity score'),
								nrow=3, ncol=1),
				nrow=2, rel_heights=c(2, 3))
		,
		cowplot::plot_grid(
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'JARID2'),
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'ZNF451'),
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'NCOR1'),
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'ADNP') ,
								get_SimiC_regulon(SimiC_aucs_PostNormal, 'SMARCE1') + ylab('Activity score'),
		nrow=5, ncol=1),
	ncol=2)
, nrow=2, rel_heights=c(1, 1.3)),
legend, nrow=2, rel_heights=c(5, 0.2))
dev.off()



##############################################################################

