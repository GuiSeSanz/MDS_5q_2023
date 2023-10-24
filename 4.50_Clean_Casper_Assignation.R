
library(Seurat)
library(ggplot2)

HL_signature <- function(data, control_type){
	signature_df <- get_genes(control_type)
	data <- data[rownames(data) %in% signature_df$gene_name,]
    logpluspointfive <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(data, 2, logpluspointfive))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))
	zscore$gene_name <- rownames(zscore)
    ann_markers <- merge(zscore, signature_df, by='gene_name')
	ann_markers <- ann_markers[ !duplicated(ann_markers$gene_name),]
	res <- ann_markers[, !colnames(ann_markers) %in% c('gene_name', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
	rownames(res) <- ann_markers$gene_name
	res <- setNames(as.data.frame(colSums(res[, !colnames(res) %in% c('gene_name', 'log2FoldChange', 'padj'), drop=FALSE], )), c('High_pondered'))
} 

get_genes <- function(sample_name_filter){
	de_files <- c(  'HSC5qvsElderly.tsv',
					'HSCs5qvsnon5q.tsv',
					'MEPs5qvsElderly.csv',
					'MEPs5qvsnon5q.tsv',
					'CMPs5qvsElderly.tsv',
					'CMPs5qvsnon5q.tsv', 
					'GMPs5qvsElderly.tsv',
					'GMPs5qvsnon5q.tsv')

	data_folder <- '/home/sevastopol/data/gserranos/MDS/Data/Annotation'
	# sample_name_filter <- elder or non5q
	genes <- data.frame(X=NULL, log2FoldChange=NULL)
	for (fl in de_files){
		if (stringr::str_detect(fl, stringr::fixed(sample_name_filter, ignore_case=TRUE))){
			tmp <- read.table(paste(data_folder, fl, sep = '/'), sep=',', header=TRUE)
			tmp <- tmp[!is.na(tmp$padj),]
			tmp <- tmp[tmp$padj<0.05, c('X', 'log2FoldChange')]
			genes <- rbind(genes, tmp)
		}
	}
	# print(dim(genes))
	genes <- setNames(genes, c('gene_name', 'log2FoldChange'))
	return(genes)
}


get_annotation_sofia_04 <- function(){
	annotation_Sofia_04 <- list()
	annotation_Sofia_04[['0' ]] <- 'HSC'
	annotation_Sofia_04[['1' ]] <- 'HSC'
	annotation_Sofia_04[['2' ]] <- 'EL'
	annotation_Sofia_04[['3' ]] <- 'E-E'
	annotation_Sofia_04[['4' ]] <- 'CFU-MK/MEP'
	annotation_Sofia_04[['5' ]] <- 'LMPP/CLP/T'
	annotation_Sofia_04[['6' ]] <- 'Mono/DC'
	annotation_Sofia_04[['7' ]] <- 'CLP/ProB'
	annotation_Sofia_04[['8' ]] <- 'Unknown'
	annotation_Sofia_04[['9' ]] <- 'CLP'
	annotation_Sofia_04[['10']] <- 'LMPP'
	annotation_Sofia_04[['11']] <- 'GMP/Granul'
	annotation_Sofia_04[['12']] <- 'E-E'
	annotation_Sofia_04[['13']] <- 'Unknown'
	annotation_Sofia_04[['14']] <- 'CLP/ProB'
	annotation_Sofia_04[['15']] <- 'Baso'
	annotation_Sofia_04[['16']] <- 'Unknown'
	return(annotation_Sofia_04)
}


annotate_clusters <- function(data, annotation_list){
	cluster_col <- grep('Cluster', colnames(data))
	cluster_names <- c()
	pb <- progress::progress_bar$new(total = nrow(data), format="[:bar] :percent eta: :eta")
	for (cl in data[[cluster_col]]){
		cluster_names <- c(cluster_names , annotation_list[[as.character(cl)]])
		pb$tick()
	}
	data$ClusterName <- cluster_names
	return(data)
	
}

get_umap_annot <- function(plotter){
	colors <- ggthemes::tableau_color_pal('Classic 20')(length(unique(plotter$ClusterName)))
	p<- ggplot(plotter, aes(x=UMAP_1, y=UMAP_2, color = ClusterName)) + geom_point(alpha=0.6 , size =1) + 
		theme_void() + scale_color_manual(name='',values=colors) + theme(legend.position = 'bottom')
	return(p)
}


get_dist <- function(plotter){
	p <- ggplot(plotter, aes(x=ClusterName,fill=ClusterName)) + geom_bar() + theme_classic() +
	theme(legend.position='none', axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),  axis.title.x= element_blank()) +
	scale_fill_manual(values=ggthemes::tableau_color_pal('Classic 20')(length(unique(plotter$ClusterName)))) + 
	facet_wrap(~State, nrow=2)
	return(p)
}

combined.sct_geneset <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
data <- as.data.frame(combined.sct_geneset@assays$SCT@data)
signature_non_5q <- HL_signature(data,'non5q')

signature_non_5q$Sample <- stringr::str_extract(rownames(signature_non_5q), '^[A-Z0-9]+')
signature_non_5q$Case <- ifelse(stringr::str_detect(signature_non_5q$Sample, '^GSM') , 'Elderly', '5q')


all_seurat_integrated_sct <- readRDS('./Data/all_seurat_integrated_sct.rds')
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
# annotation resolution is 0.4
coords <- merge(coords, setNames(as.data.frame(all_seurat_integrated_sct$integrated_snn_res.0.4), c('Cluster_04')), by=0)




# # get the overlapping cells between the signature and casper
# selected_cells_consensus <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL) 
# # get the cells selected by casper but no in the signature
# discarded_cells_casper <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL)
# # get the cells selected by signature but no in the casper
# discarded_cells_signature <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL)
# percentages_2_keep = c(5, 10, 20, 30)
# for (n in c(percentages_2_keep)){
# 	for (sample_name in unique(all_seurat_integrated_sct$Sample)){
# 		print(paste0('Processing the ', n , '% for  ',sample_name))
# 		coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]
# 		chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'.rds'))
# 		casper_results <- reshape::melt(chrMat)
# 		casper_results$value2 <- "neutral"
# 		casper_results$value2[casper_results$value > 0] <- "amplification"
# 		casper_results$value2[casper_results$value < 0] <- "deletion"
# 		casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
# 			"deletion", "neutral"))
# 		casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
# 		casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')


# 		signature_non_5q_tmp <- signature_non_5q[rownames(signature_non_5q) %in% coords_smp$Row.names ,]
# 		signature_non_5q_tmp_top    <-  signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),]
# 		signature_non_5q_tmp_bottom <-  signature_non_5q_tmp[signature_non_5q_tmp$High_pondered < quantile(signature_non_5q_tmp$High_pondered,prob=n/100),]
# 		# These cells come from the counts 
# 		casper_results_del <- casper_results[casper_results$value2 == 'deletion' & casper_results$X2 == '5q' & casper_results$phenotype == sample_name,'X1'] 
# 		casper_results_amp <- casper_results[casper_results$value2 == 'amplification' & casper_results$X2 == '5q' & casper_results$phenotype == sample_name,'X1'] 
# 		casper_results_neu <- casper_results[casper_results$value2 == 'neutral'  & casper_results$X2 == '5q' & casper_results$phenotype == sample_name,'X1']
# 		casper_results_del <- casper_results_del[casper_results_del %in% coords_smp$Row.names]
# 		casper_results_neu <- casper_results_neu[casper_results_neu %in% coords_smp$Row.names]
# 		if(length(casper_results_del) + length(casper_results_neu) + length(casper_results_amp) != nrow(signature_non_5q_tmp)){
# 			print('Some dimensions are not correct... :c')
# 		}
# 		# get the cells overlaping results casper and sign
# 		true_deleted_cells <- intersect(casper_results_del, rownames(signature_non_5q_tmp_top))
# 		true_neutral_cells <- intersect(casper_results_neu, rownames(signature_non_5q_tmp_bottom))
# 		selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_deleted_cells, State='deleted', Sample = sample_name, Selected = n))
# 		selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_neutral_cells, State='neutral', Sample = sample_name, Selected = n))
# 		# get the cells predicted by casper but not the signature
# 		tmp_disc_casper <- setdiff(casper_results_del, rownames(signature_non_5q_tmp_top))
# 		discarded_cells_casper <- rbind(discarded_cells_casper, data.frame(cell_id =tmp_disc_casper, State='deleted', Sample = sample_name, Selected = n ))
# 		tmp_disc_casper <- setdiff(casper_results_neu, rownames(signature_non_5q_tmp_bottom))
# 		discarded_cells_casper <- rbind(discarded_cells_casper, data.frame(cell_id =tmp_disc_casper, State='neutral', Sample = sample_name, Selected = n ))
# 		# get the cells predicted by signature but not casper
# 		tmp_disc_sig <- setdiff(rownames(signature_non_5q_tmp_top), casper_results_del)
# 		discarded_cells_signature <- rbind(discarded_cells_signature,  data.frame(cell_id =tmp_disc_sig, State='deleted', Sample = sample_name, Selected = n ))
# 		tmp_disc_sig <- setdiff(rownames(signature_non_5q_tmp_bottom), casper_results_neu)
# 		discarded_cells_signature <- rbind(discarded_cells_signature,  data.frame(cell_id =tmp_disc_sig, State='neutral', Sample = sample_name, Selected = n ))

# 	}
# }




casper_results <- function(samples){
	all_results <- data.frame(X1=NULL, X2=NULL, value=NULL, value2=NULL, phenotype=NULL)
	for (sample_name in samples){
		print(paste0('Processing ',sample_name))
		coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]
		chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'.rds'))
		casper_results <- reshape::melt(chrMat)
		casper_results$value2 <- "neutral"
		casper_results$value2[casper_results$value > 0] <- "amplification"
		casper_results$value2[casper_results$value < 0] <- "deletion"
		casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
			"deletion", "neutral"))
		casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
		casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')
		all_results <- rbind(all_results, casper_results)
	}
	return(all_results)
}


casper_results <- casper_results(unique(all_seurat_integrated_sct$Sample))
# get the overlapping cells between the signature and casper
selected_cells_consensus <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL) 
# get the cells selected by casper but no in the signature
discarded_cells_casper <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL)
# get the cells selected by signature but no in the casper
discarded_cells_signature <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL)
percentages_2_keep = c(5, 10, 20, 30)
for (n in c(percentages_2_keep)){
	signature_non_5q_mds <- signature_non_5q[signature_non_5q$Sample %in% c(unique(all_seurat_integrated_sct$Sample)),]
	signature_non_5q_tmp <- signature_non_5q_mds[rownames(signature_non_5q_mds) %in% coords$Row.names ,]
	signature_non_5q_tmp_top    <-  signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),]
	signature_non_5q_tmp_bottom <-  signature_non_5q_tmp[signature_non_5q_tmp$High_pondered < quantile(signature_non_5q_tmp$High_pondered,prob=n/100),]
	# These cells come from the counts 
	casper_results_del <- casper_results[casper_results$value2 == 'deletion' & casper_results$X2 == '5q','X1'] 
	casper_results_amp <- casper_results[casper_results$value2 == 'amplification' & casper_results$X2 == '5q','X1'] 
	casper_results_neu <- casper_results[casper_results$value2 == 'neutral'  & casper_results$X2 == '5q','X1']
	# casper_results_del <- casper_results_del[casper_results_del %in% coords_smp$Row.names]
	# casper_results_neu <- casper_results_neu[casper_results_neu %in% coords_smp$Row.names]
	if(length(casper_results_del) + length(casper_results_neu) + length(casper_results_amp) != nrow(signature_non_5q_tmp)){
		print('Some dimensions are not correct... :c')
	}
	# get the cells overlaping results casper and sign
	true_deleted_cells <- intersect(casper_results_del, rownames(signature_non_5q_tmp_top))
	true_neutral_cells <- intersect(casper_results_neu, rownames(signature_non_5q_tmp_bottom))
	selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_deleted_cells, State='deleted', Sample = sample_name, Selected = n))
	selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_neutral_cells, State='neutral', Sample = sample_name, Selected = n))
	# get the cells predicted by casper but not the signature
	tmp_disc_casper <- setdiff(casper_results_del, rownames(signature_non_5q_tmp_top))
	discarded_cells_casper <- rbind(discarded_cells_casper, data.frame(cell_id =tmp_disc_casper, State='deleted', Sample = sample_name, Selected = n ))
	tmp_disc_casper <- setdiff(casper_results_neu, rownames(signature_non_5q_tmp_bottom))
	discarded_cells_casper <- rbind(discarded_cells_casper, data.frame(cell_id =tmp_disc_casper, State='neutral', Sample = sample_name, Selected = n ))
	# get the cells predicted by signature but not casper
	tmp_disc_sig <- setdiff(rownames(signature_non_5q_tmp_top), casper_results_del)
	discarded_cells_signature <- rbind(discarded_cells_signature,  data.frame(cell_id =tmp_disc_sig, State='deleted', Sample = sample_name, Selected = n ))
	tmp_disc_sig <- setdiff(rownames(signature_non_5q_tmp_bottom), casper_results_neu)
	discarded_cells_signature <- rbind(discarded_cells_signature,  data.frame(cell_id =tmp_disc_sig, State='neutral', Sample = sample_name, Selected = n ))

}

pdf('./Plots/CASPER/CASPER_results_consensus_signature.pdf')
# print the UMAP with all the cells to compare the results
	tmp <- annotate_clusters(coords, get_annotation_sofia_04())
	get_umap_annot(tmp) +  ggtitle(paste0('Original n_cells:', nrow(tmp)))

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 30,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 30,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))
	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_30 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_30') + facet_wrap(~State) + theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 20,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 20,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_20 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_20') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 10,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 10,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_10 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_10') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 5,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 5,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_5 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_5') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

dev.off()



pdf('./Plots/CASPER/CASPER_results_casper_no_signature.pdf')
# print the UMAP with all the cells to compare the results
	tmp <- annotate_clusters(coords, get_annotation_sofia_04())
	get_umap_annot(tmp) +  ggtitle(paste0('Original n_cells:', nrow(tmp)))

	tmp <- coords[coords$Row.names %in% discarded_cells_casper[discarded_cells_casper$Selected == 30,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_casper[discarded_cells_casper$Selected == 30,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))
	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_30 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_30') + facet_wrap(~State) + theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_casper[discarded_cells_casper$Selected == 20,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_casper[discarded_cells_casper$Selected == 20,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_20 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_20') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_casper[discarded_cells_casper$Selected == 10,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_casper[discarded_cells_casper$Selected == 10,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_10 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_10') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_casper[discarded_cells_casper$Selected == 5,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_casper[discarded_cells_casper$Selected == 5,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_5 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_5') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

dev.off()



pdf('./Plots/CASPER/CASPER_results_signature_no_casper.pdf')
# print the UMAP with all the cells to compare the results
	tmp <- annotate_clusters(coords, get_annotation_sofia_04())
	get_umap_annot(tmp) +  ggtitle(paste0('Original n_cells:', nrow(tmp)))

	tmp <- coords[coords$Row.names %in% discarded_cells_signature[discarded_cells_signature$Selected == 30,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_signature[discarded_cells_signature$Selected == 30,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))
	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_30 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_30') + facet_wrap(~State) + theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_signature[discarded_cells_signature$Selected == 20,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_signature[discarded_cells_signature$Selected == 20,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_20 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_20') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_signature[discarded_cells_signature$Selected == 10,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_signature[discarded_cells_signature$Selected == 10,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_10 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_10') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% discarded_cells_signature[discarded_cells_signature$Selected == 5,'cell_id'],]
	tmp <- merge(tmp, discarded_cells_signature[discarded_cells_signature$Selected == 5,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_5 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_5') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

dev.off()






pdf('./Plots/CASPER/CASPER_results_signature_percluster.pdf')
tmp <- merge(coords, signature_non_5q, by.x='Row.names', by.y=0)
tmp <- annotate_clusters(tmp, get_annotation_sofia_04())
ggplot(tmp, aes(x=ClusterName, y= High_pondered, fill = ClusterName)) + geom_boxplot() + theme_classic() +
scale_fill_manual(values =ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$ClusterName)))) + 
theme(legend.position='none', axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggplot(tmp, aes(x=High_pondered, y= ClusterName, fill = ClusterName)) +  ggridges::geom_density_ridges2() + theme_classic() +
scale_fill_manual(values =ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$ClusterName)))) + 
theme(legend.position='none', axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
ggplot(tmp, aes(x=High_pondered, y= ClusterName, fill = ClusterName)) +  ggridges::geom_density_ridges2() + theme_classic() +
scale_fill_manual(values =ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$ClusterName)))) + 
theme(legend.position='none', axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) + facet_wrap(~Sample)
dev.off()







casper_results <- function(samples){
	all_results <- data.frame(X1=NULL, X2=NULL, value=NULL, value2=NULL, phenotype=NULL)
	for (sample_name in samples){
		print(paste0('Processing the',sample_name))
		coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]
		chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'.rds'))
		casper_results <- reshape::melt(chrMat)
		casper_results$value2 <- "neutral"
		casper_results$value2[casper_results$value > 0] <- "amplification"
		casper_results$value2[casper_results$value < 0] <- "deletion"
		casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
			"deletion", "neutral"))
		casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
		casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')
		all_results <- rbind(all_results, casper_results)
	}
	return(all_results)
}

casper_predictions <- casper_results(unique(all_seurat_integrated_sct$Sample))
ann_04 <- get_annotation_sofia_04()

pdf('./Plots/CASPER/CASPER_results.pdf')
	tmp <- merge(coords, casper_predictions[casper_predictions$X2 == '5q' & casper_predictions$value2 %in% c('neutral', 'deletion'),], by.x='Row.names', by.y='X1')
	tmp$State <- tmp$value2
	tmp$ClusterName <- apply(tmp, 1, FUN=function(x) ann_04[[as.character(x[['Cluster_04']])]])
	print(dim(tmp))
	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) , ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + facet_wrap(~State) + theme(legend.position='none'),
		nrow=2
	)
dev.off()







#### New Signatures Battery


casper_results <- function(samples){
	all_results <- data.frame(X1=NULL, X2=NULL, value=NULL, value2=NULL, phenotype=NULL)
	for (sample_name in samples){
		print(paste0('Processing ',sample_name))
		coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]
		chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'.rds'))
		casper_results <- reshape::melt(chrMat)
		casper_results$value2 <- "neutral"
		casper_results$value2[casper_results$value > 0] <- "amplification"
		casper_results$value2[casper_results$value < 0] <- "deletion"
		casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
			"deletion", "neutral"))
		casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
		casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')
		all_results <- rbind(all_results, casper_results)
	}
	return(all_results)
}


casper_results <- casper_results(unique(all_seurat_integrated_sct$Sample))
# get the overlapping cells between the signature and casper

selected_cells_consensus <- data.frame(cell_id = NULL, State=NULL, Sample = NULL, Selected=NULL) 

percentages_2_keep = c(10, 20, 30, 50)

for (n in c(percentages_2_keep)){
	print(paste0(n, '%'))
	Signatures_mds <- Signatures[Signatures$Sample %in% c(unique(all_seurat_integrated_sct$Sample)),]
	Signatures_mds <- Signatures_mds[rownames(Signatures_mds) %in% coords$Row.names ,]
	Signatures_mds    <-  Signatures_mds[Signatures_mds$Pearson_Rank > quantile(Signatures_mds$Pearson_Rank,prob=1-n/100),]
	# These cells come from the counts 
	casper_results_smp <- casper_results[casper_results$X1 %in% coords$Row.names,]
	casper_results_del <- casper_results[casper_results$value2 == 'deletion' & casper_results$X2 == '5q','X1'] 
	casper_results_amp <- casper_results[casper_results$value2 == 'amplification' & casper_results$X2 == '5q','X1'] 
	casper_results_neu <- casper_results[casper_results$value2 == 'neutral'  & casper_results$X2 == '5q','X1']

	# get the cells overlaping results casper and sign
	true_deleted_cells <- intersect(casper_results_del, rownames(Signatures_mds))
	true_neutral_cells <- casper_results_neu
	selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_deleted_cells, State='deleted',  Selected = n))
	selected_cells_consensus <- rbind(selected_cells_consensus, data.frame(cell_id = true_neutral_cells, State='neutral',  Selected = n))
	
}





pdf('./Plots/CASPER/CASPER_results_consensus_signaturePearson.pdf')
# print the UMAP with all the cells to compare the results
	get_umap_annot(coords) +  ggtitle(paste0('Original n_cells:', nrow(coords)))


	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 50,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 50,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_50 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_50') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)


	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 30,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 30,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	print(dim(tmp))
	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_30 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_30') + facet_wrap(~State) + theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 20,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 20,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_20 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_20') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

	tmp <- coords[coords$Row.names %in% selected_cells_consensus[selected_cells_consensus$Selected == 10,'cell_id'],]
	tmp <- merge(tmp, selected_cells_consensus[selected_cells_consensus$Selected == 10,c('cell_id', 'State')], by.x = 'Row.names', by.y = 'cell_id')
	print(dim(tmp))

	cowplot::plot_grid(
		cowplot::plot_grid(get_dist(tmp), NULL, get_umap_annot(tmp) + ggtitle(paste0('Selected_10 n_cells:', nrow(tmp))), ncol=3, rel_widths = c(0.3,0.1,0.6)),
		get_umap_annot(tmp) + ggtitle('Selected_10') + facet_wrap(~State)+ theme(legend.position='none'),
		nrow=2
	)

dev.off()

