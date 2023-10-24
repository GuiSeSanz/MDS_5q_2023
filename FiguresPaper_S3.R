
library(ggplot2)

all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')
all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')

all_DE_pseudo_Non5qCR_Vs_Elder   <- readRDS('./Data/Results_DE_Non5qCR_Vs_Elder.rds')
all_DE_pseudo_Non5qPR_Vs_Elder <- readRDS('./Data/Results_DE_Non5qPR_VsElder.rds')

all_DE_Non5qCR_Vs_non5qDiagnosed <- readRDS('./Data/Results_DE_Non5qCR_Vs_non5qDiagnosed.rds')
all_DE_Non5qPR_CR_Vs_non5qDiagnosed <- readRDS('./Data/Results_DE_Non5qPR_CR_Vs_non5qDiagnosed.rds')

all_DE_NR_Vs_MDS <- readRDS('./Data/Results_DE_del5qNR_Post_Vs_del5qMDS.rds')

all_DE_PR_Vs_del5qNR <- readRDS('./Data/Results_DE_Del5qPRVsDel5qNR_MAST.rds') # barplots + volcanos


all_DE_Non5qPR_Vs_non5qMDS <- readRDS('./Data/Results_DE_Non5qPR_Vs_non5qDiagnosed.rds')

all_DE_del5qPR_Vs_del5qMDS <- readRDS('./Data/Results_DE_del5qPR_Vs_del5qMDS.rds')

get_barplots_DEG <- function(data_de, comparison_title=NULL, logFC_cutoff=2 ){
	for (cell_type in names(data_de)){
		data_de[[cell_type]]$cell_type <- cell_type
	}
	data_de <- do.call("rbind", data_de)
	if (!'avg_logFC' %in% names(data_de)){
		data_de$avg_logFC <- data_de$avg_log2FC
	}
	plotter <- data_de[abs(data_de$avg_logFC) >logFC_cutoff,]
	plotter$is_pos <- ifelse(plotter$avg_logFC >=0, 'UP', 'DOWN')
	plotter <- as.data.frame(table(plotter$is_pos, plotter$cell_type))
	plotter$Freq <- ifelse(plotter$Var1 == 'UP', plotter$Freq, plotter$Freq*-1)
	p <- ggplot(plotter, aes(x=Var2, y=Freq,  fill=Var1)) + 
	geom_hline(yintercept = 0, color='black') + 
	geom_bar(stat = "identity", position = "stack", color='black') + theme_classic() + 
		scale_fill_manual(values=c(DOWN='#5074AF', UP='#C64032')) + 
	labs(y='Number of DEG') + 
	guides(fill=guide_legend(title=paste0("FDR<0.05\n|logFC|>=", logFC_cutoff)))+ 
	scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) +
	labs(title=comparison_title) +
	theme(
		plot.title = element_text(hjust=0.5, size=10, face='italic'), 
		axis.text.x=element_text(angle=90, hjust=1), 
		axis.title=element_blank(), 
		legend.position='none') 
	return(p)
}

get_volcano <- function(results, CType, sub= "5q Vs Elder"){
	p <-  EnhancedVolcano::EnhancedVolcano(results,
				lab = results$gene_name,
				x = 'avg_logFC',
				y = 'p_val_adj',col=c('black', 'black', 'black', 'red3'),
				subtitle = NULL, titleLabSize=2,
				pointSize=0.3, colAlpha = 1, 
				caption = NULL, xlim=c(-3,3),
				axisLabSize = 10, title= CType,
				pCutoff = 0.05, borderWidth= 0.5, 
				FCcutoff = 2) +
				theme(legend.position='none', 
				text = element_text(family = "Helvetica", size = 5),
				axis.title= element_text(family = "Helvetica", size = 12))
	return(p)
}



noaxis_title <- theme(axis.title=element_blank())
noaxis_y_title <- theme(axis.title.y=element_blank())
noaxis_x_title <- theme(axis.title.x=element_blank())

volcanos_list <- list()
for (ct in names(all_DE_PR_Vs_del5qNR)){
	tt <- all_DE_PR_Vs_del5qNR[[ct]]
	tt$gene_name <- rownames(tt)
	tt$avg_logFC <- tt$avg_log2FC
	volcanos_list[[ct]] <- get_volcano(tt, ct, sub='')
}




# pdf('./Plots/PaperFigures/FigS3.pdf', width=8.3, height=11.7)
# cowplot::plot_grid(

# 	cowplot::plot_grid(
# 		get_barplots_DEG(all_results_psuedo_5qVsElder), 
# 		get_barplots_DEG(all_results_psuedo_Non5qVsElder),
# 		get_barplots_DEG(all_DE_PR_Vs_del5qNR, logFC_cutoff=1)
# 	,nrow=3, labels=c('A', 'D', 'G')), 

# 	cowplot::plot_grid(
# 		cowplot::plot_grid( 
# 							get_barplots_DEG(all_DE_pseudo_Non5qCR_Vs_Elder), 
# 							get_barplots_DEG(all_DE_pseudo_Non5qPR_Vs_Elder),
# 							get_barplots_DEG(all_DE_Non5qCR_Vs_non5qDiagnosed),
# 							get_barplots_DEG(all_DE_Non5qPR_CR_Vs_non5qDiagnosed),
# 		ncol=2, labels=c('B', 'C', 'E', 'F')),
# 		cowplot::plot_grid(
# 					volcanos_list[[1]] + noaxis_title,  volcanos_list[[2]] + noaxis_title, volcanos_list[[3]] + noaxis_title, volcanos_list[[4]] + noaxis_title,
# 					volcanos_list[[5]]+ noaxis_title,  volcanos_list[[6]] + noaxis_title, volcanos_list[[7]] + noaxis_title, volcanos_list[[8]] + noaxis_title,
# 					volcanos_list[[9]] + noaxis_title, volcanos_list[[10]] + noaxis_title, volcanos_list[[11]] + noaxis_title,
# 					ncol=4,labels=c('H'))
# 		, nrow=2, rel_heights=c(2,1))

# ,ncol=2, rel_widths=c(1,2))
# dev.off()

# pdf('./Plots/PaperFigures/FigS3.pdf', width=8.3, height=11.7)
# cowplot::plot_grid(

# 	cowplot::plot_grid(
# 		get_barplots_DEG(all_results_psuedo_5qVsElder), 
# 		get_barplots_DEG(all_results_psuedo_Non5qVsElder),
# 		get_barplots_DEG(all_DE_PR_Vs_del5qNR, logFC_cutoff=1)
# 	,nrow=3, labels=c('A', 'D', 'G')), 

# 	cowplot::plot_grid(
# 		cowplot::plot_grid( 
# 							get_barplots_DEG(all_DE_pseudo_Non5qCR_Vs_Elder), 
# 							get_barplots_DEG(all_DE_pseudo_Non5qPR_Vs_Elder),
# 							get_barplots_DEG(all_DE_Non5qCR_Vs_non5qDiagnosed),
# 							get_barplots_DEG(all_DE_Non5qPR_CR_Vs_non5qDiagnosed),
# 		ncol=2, labels=c('B', 'C', 'E', 'F')),
# 		cowplot::plot_grid(
# 					volcanos_list[[1]] + noaxis_title,  volcanos_list[[2]] + noaxis_title, volcanos_list[[3]] + noaxis_title, volcanos_list[[4]] + noaxis_title,
# 					volcanos_list[[5]]+ noaxis_title,  volcanos_list[[6]] + noaxis_title, volcanos_list[[7]] + noaxis_title, volcanos_list[[8]] + noaxis_title,
# 					volcanos_list[[9]] + noaxis_title, volcanos_list[[10]] + noaxis_title, volcanos_list[[11]] + noaxis_title,
# 					ncol=4,labels=c('H'))
# 		, nrow=2, rel_heights=c(2,1))

# ,ncol=2, rel_widths=c(1,2))
# dev.off()







# pdf('./Plots/PaperFigures/FigS3.pdf', width=8.3, height=11.7)
# cowplot::plot_grid(

# 	cowplot::plot_grid(
# 		get_barplots_DEG(all_results_psuedo_5qVsElder), 
# 		get_barplots_DEG(all_results_psuedo_Non5qVsElder),
# 		get_barplots_DEG(all_DE_PR_Vs_del5qNR, logFC_cutoff=1)
# 	,nrow=3, labels=c('A', 'D', 'G')), 

# 	cowplot::plot_grid(
# 		cowplot::plot_grid( 
# 							get_barplots_DEG(all_DE_pseudo_Non5qCR_Vs_Elder), 
# 							get_barplots_DEG(all_DE_pseudo_Non5qPR_Vs_Elder),
# 							get_barplots_DEG(all_DE_Non5qCR_Vs_non5qDiagnosed),
# 							get_barplots_DEG(all_DE_Non5qPR_CR_Vs_non5qDiagnosed),
# 		ncol=2, labels=c('B', 'C', 'E', 'F')),
# 		cowplot::plot_grid(
# 					volcanos_list[[1]] + noaxis_title,  volcanos_list[[2]] + noaxis_title, volcanos_list[[3]] + noaxis_title, volcanos_list[[4]] + noaxis_title,
# 					volcanos_list[[5]]+ noaxis_title,  volcanos_list[[6]] + noaxis_title, volcanos_list[[7]] + noaxis_title, volcanos_list[[8]] + noaxis_title,
# 					volcanos_list[[9]] + noaxis_title, volcanos_list[[10]] + noaxis_title, volcanos_list[[11]] + noaxis_title,
# 					ncol=4,labels=c('H'))
# 		, nrow=2, rel_heights=c(2,1))

# ,ncol=2, rel_widths=c(1,2))
# dev.off()




pdf('./Plots/PaperFigures/FigS3.pdf', width=8.3, height=11.7)

cowplot::plot_grid(
	cowplot::plot_grid(
		get_barplots_DEG(all_results_psuedo_5qVsElder,     'del(5q) diagnosed\nVs\nHealthy'), 
		get_barplots_DEG(all_results_psuedo_Non5qVsElder,  'non-del(5q) diagnosed\nVs\nHealthy'), 
		get_barplots_DEG(all_DE_Non5qCR_Vs_non5qDiagnosed, 'non-del(5q) Complete-Responder\nVs\nnon-del(5q) diagnosed'),
		get_barplots_DEG(all_DE_Non5qPR_Vs_non5qMDS,       'non-del(5q) Partial-Responder\nVs\nnon-del(5q) diagnosed'),
		get_barplots_DEG(all_DE_pseudo_Non5qCR_Vs_Elder,   'non-del(5q) Complete-Responder\nVs\nHealthy'), 
		get_barplots_DEG(all_DE_pseudo_Non5qPR_Vs_Elder,   'non-del(5q) Partial-Responder\nVs\nHealthy'),
		get_barplots_DEG(all_DE_del5qPR_Vs_del5qMDS,       'del(5q) Partial-Responder\nVs\ndel(5q) diagnosed'),
		get_barplots_DEG(all_DE_PR_Vs_del5qNR,             'del(5q) Partial-Responder\nVs\ndel(5q) Non-Responder', logFC_cutoff=0),
		get_barplots_DEG(all_DE_NR_Vs_MDS,                 'del(5q) Non-Responder\nVs\ndel(5q) diagnosed'),
	labels='AUTO', ncol=3)
	, 
	cowplot::get_legend(get_barplots_DEG(all_DE_NR_Vs_MDS)+ theme(legend.position='bottom')),
	nrow=2, rel_heights=c(2,0.1))
dev.off()



data_de <- all_DE_PR_Vs_del5qNR
max_tmp <- c()
for (ct in names(data_de)){
	tmp <- data_de[[ct]]
	# tmp <- tmp[tmp$p_val_adj<0.05, ]
	if('avg_logFC' %in% names(tmp)){
		tmp <- tmp[tmp$p_val_adj<0.05 , ]
	}else{
		tmp <- tmp[tmp$p_val_adj<0.05 , ]
	}
	max_tmp <- c(max_tmp,nrow(tmp))
}
summary(max_tmp)


for (ct in names(all_DE_PR_Vs_del5qNR)){
	tmp <- all_DE_PR_Vs_del5qNR[[ct]]
	tmp$gene_name <- rownames(tmp)
	#sort the colnames having gene_name at the beginning
	all_DE_PR_Vs_del5qNR[[ct]] <- tmp 
}

all_data <- list(
	do.call('rbind', all_results_psuedo_5qVsElder),
	do.call('rbind', all_results_psuedo_Non5qVsElder),
	do.call('rbind', all_DE_Non5qCR_Vs_non5qDiagnosed),
	do.call('rbind', all_DE_Non5qPR_Vs_non5qMDS),
	do.call('rbind', all_DE_pseudo_Non5qCR_Vs_Elder),
	do.call('rbind', all_DE_pseudo_Non5qPR_Vs_Elder),
	do.call('rbind', all_DE_del5qPR_Vs_del5qMDS),
	do.call('rbind', all_DE_PR_Vs_del5qNR),
	do.call('rbind', all_DE_NR_Vs_MDS)
)

names(all_data) <- c(
	'del(5q)diagnosedVsHealthy',
	'non-del(5q)diagnosedVsHealthy',
	'non-del(5q)CRVsnon-del(5q)diagn',
	'non-del(5q)PCRVsnon-del(5q)diag',
	'non-del(5q)CRVsHealthy',
	'non-del(5q)PCRVsHealthy',
	'del(5q)PCRVsdel(5q)diagnosed',
	'del(5q)PRVsdel(5q)NR',
	'del(5q)NRVsdel(5q)diagnosed'
)

for (name in names(all_data)){
	print(nchar(name))
}
WriteXLS::WriteXLS(all_data, './Plots/PaperFigures/Suppl_Table_3.xlsx', col.names=TRUE, BoldHeaderRow=TRUE, row.names=FALSE, SheetNames=names(all_data))


