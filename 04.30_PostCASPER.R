

library(Seurat)
library(CaSpER)
library(ggplot2)
data("hg38_cytoband")




get_deletions_chromosomes <- function(data){
	control <- 'SMD35344'
	sample <-  setdiff(unique(casper_results$phenotype), control)
	names5q <- ifelse(levels(data$X2) %in% c('5q', '5p'), 'red', 'black')
	data$phenotype <- factor(data$phenotype, levels = c(sample, control))
	p <- ggplot(data, aes(x=X2, fill = value2)) + geom_bar(position='fill') + theme_classic() + facet_wrap(.~phenotype, nrow=2) +
		scale_fill_manual(values=c('neutral'="#999999", 'amplification'="#E69F00", 'deletion'="#56B4E9"))+ 
		theme(axis.text.x = element_text(angle = 45, hjust=1, size=6, color = names5q)) +
		ggtitle('Distribution of the CNV per chromosome')
	return(p)
}

get_umap <- function(coordinates, data){
	plotter <- merge(coordinates, data, by.x='Row.names', by.y='X1')
	p <- ggplot(plotter, aes(x = UMAP_1, y = UMAP_2))  + 
	geom_point(data = plotter[plotter$value2 == 'neutral',],  color='#999999', alpha=0.2, size=1) + 
	geom_point(data = plotter[plotter$value2 == 'amplification',],  color='#E69F00', alpha=0.5,  size=1) + 
	geom_point(data = plotter[plotter$value2 == 'deletion',],  color='#56B4E9', alpha=0.5,  size=1) + 
	theme_classic() + facet_wrap(~X2) + theme(legend.position="bottom") +
	ggtitle('Umap of the chr5 CNV affected cells')
	return(p)
}

get_cluster_deletion <- function(coordinates, data){
	plotter <- merge(coordinates, data, by.x='Row.names', by.y='X1')
	p <- ggplot(plotter, aes(x = Cluster, fill = value2)) + geom_bar(position='fill') + theme_classic() +
		scale_fill_manual(values=c('neutral'="#999999", 'amplification'="#E69F00", 'deletion'="#56B4E9"))+ 
		theme(axis.text.x = element_text(angle = 45, hjust=1, size=6)) + 
		ggtitle('CNV per cluster, resolution 0.4')
	return(p)
}

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

get_venn <- function(data, percent = TRUE){
	p <- ggvenn::ggvenn(data,
			fill_color = destiny::cube_helix(length(data)),
			stroke_size = 0.4,
			show_percentage = percent,
			fill_alpha = 0.4,
			stroke_color = 'white',
			stroke_alpha = 1,
			stroke_linetype = 'solid',
			text_color = 'black',
			set_name_size = 4,
			text_size = 3.5,
			label_sep = ','
			)+ theme(plot.title = element_text(hjust = 0.5))
	return(p)
}


get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)), 
	order.by = "freq", empty.intersections = "on")
	plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
			nrow=2, align='hv', rel_heights = c(3,1),
			rel_widths = c(2,3))
	return(plot)
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

combined.sct_geneset <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
data <- as.data.frame(combined.sct_geneset@assays$SCT@data)
signature_non_5q <- HL_signature(data,'non5q')
signature_elder  <- HL_signature(data,'elder')

all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
# annotation resolution is 0.4
coords <- merge(coords, setNames(as.data.frame(all_seurat_integrated_sct$integrated_snn_res.0.4), c('Cluster')), by=0)

deletion_percentages <- data.frame(Samples = c('SMD34459',	'SMD35109',	'SMD35303',	'SMD37209', 'SMD29402'), 
								Karyotype = c(43,75,90,30, 0), 
								CASPER= c(0,0,0,0, 0))

get_pval <- function(casper_num, signature_num, all_res){
	q <- length(intersect(casper_num,signature_num))
	k <- length(signature_num)
	m <- length(casper_num)
	n <- length(all_res[all_res$X2 %in% c('5q'),'X1']) - m
	pval <- phyper(q, m, n, k, lower.tail=FALSE)
	return(pval)
}


all_5q_depleted_cells <- c()
for (sample_name in c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303', 'SMD29402')){	print(sample_name)
	# sample_name <- 'SMD35109'
	coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]

	chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'_proper_control.rds'))
	casper_results <- melt(chrMat)
	casper_results$value2 <- "neutral"
	casper_results$value2[casper_results$value > 0] <- "amplification"
	casper_results$value2[casper_results$value < 0] <- "deletion"
	casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
		"deletion", "neutral"))
	casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
	casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')
	deletion_percentages[deletion_percentages$Samples == sample_name,'CASPER'] <- (prop.table(table(casper_results[casper_results$phenotype != 'SMD35344' & casper_results$ X2 == '5q', 'value2']))*100)[['deletion']]
	all_5q_depleted_cells <- c(all_5q_depleted_cells, as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1']))
	pdf(paste0('./Plots/CASPER/Test_Results_CASPER_',sample_name,'.pdf'), height=5)
	print(get_deletions_chromosomes(casper_results))
	if (sample_name != 'SMD29402'){
		print(get_umap(coords_smp, casper_results[stringr::str_detect(casper_results$X1, sample_name) & casper_results$X2 %in% c('5p', '5q'),]))
		print(get_cluster_deletion(coords_smp, casper_results))
		casper_results_5q <- casper_results[casper_results$X2 %in% c('5q'),]
		print(get_cluster_deletion(coords_smp[coords_smp$Row.names %in% casper_results_5q$X1,], casper_results_5q) + ggtitle('Chr5q CNV per cluster, resolution 0.4'))
	} 
	dev.off()
	
	# create the 
	signature_non_5q_tmp <- signature_non_5q[stringr::str_detect(rownames(signature_non_5q),sample_name), , drop=FALSE]
	signature_non_5q_tmp$X1 <- rownames(signature_non_5q_tmp)
	pdf(paste0('./Plots/CASPER/Venn_Results_CASPER_',sample_name,'.pdf'))
	n <- 50
	data_venn <- list()
	data_venn[['casper_5q_del']] <- as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1'])
	data_venn[['top_50_signature']] <- signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),'X1']
	pval <- get_pval(data_venn[['casper_5q_del']], data_venn[['top_50_signature']] , casper_results)
	print(get_venn(data_venn)+ggtitle(paste0('top50% pval::', pval)))
	n <- 20
	data_venn <- list()
	data_venn[['casper_5q_del']] <- as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1'])
	data_venn[['top_20_signature']] <- signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),'X1']
	print(get_venn(data_venn)+ggtitle(paste0('top20% pval::', pval)))
	n <- 10
	data_venn <- list()
	data_venn[['casper_5q_del']] <- as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1'])
	data_venn[['top_10_signature']] <- signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),'X1']
	print(get_venn(data_venn)+ggtitle(paste0('top10% pval::', pval)))
	n <- 5
	data_venn <- list()
	data_venn[['casper_5q_del']] <- as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1'])
	data_venn[['top_5_signature']] <- signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),'X1']
	print(get_venn(data_venn)+ggtitle(paste0('top5% pval::', pval)))
	n <- 1
	data_venn <- list()
	data_venn[['casper_5q_del']] <- as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1'])
	data_venn[['top_1_signature']] <- signature_non_5q_tmp[signature_non_5q_tmp$High_pondered > quantile(signature_non_5q_tmp$High_pondered,prob=1-n/100),'X1']
	print(get_venn(data_venn)+ggtitle(paste0('top1% pval::', pval)))
	dev.off()	

}

pdf(paste0('./Plots/CASPER/CASPER_percentages.pdf'))
gridExtra::grid.table(deletion_percentages, rows=NULL)
dev.off()

pdf('./Plots/CASPER/CASPER_percentage_corr.pdf')
ggplot(deletion_percentages, aes(x=Karyotype, y= CASPER, color=Samples)) + geom_point(size = 6)  + 
ggprism::theme_prism() + scale_color_manual(values=ggthemes::tableau_color_pal('Classic 10 Medium')(5)) + 
theme(legend.position='bottom') + geom_abline(intercept = 0)
dev.off()
# table(casper_results[casper_results$phenotype != 'SMD35344' &casper_results$ X2 == '5q', 'value2'])

# venn diagram for all the DE genes
de_files <- c(  'HSC5qvsElderly.tsv',
				'HSCs5qvsnon5q.tsv',
				'MEPs5qvsElderly.csv',
				'MEPs5qvsnon5q.tsv',
				'CMPs5qvsElderly.tsv',
				'CMPs5qvsnon5q.tsv', 
				'GMPs5qvsElderly.tsv',
				'GMPs5qvsnon5q.tsv')
data_folder <- '/home/sevastopol/data/gserranos/MDS/Data/Annotation'

de_genes_vsNon5q <- list()
de_genes_vsElder <- list()
for (fl in de_files){
	cell_type <- stringr::str_extract(fl, '^[A-Z0-9^]+')
	tmp <- read.table(paste(data_folder, fl, sep = '/'), sep=',', header=TRUE)
	tmp <- tmp[!is.na(tmp$padj),]
	tmp <- as.character(tmp[tmp$padj<0.05 & tmp$log2FoldChange !=0, 'X'])
	if(stringr::str_detect(fl, 'Elder')){
		de_genes_vsElder[[cell_type]] <- tmp
	}else if( stringr::str_detect(fl, 'non5q')){
		de_genes_vsNon5q[[cell_type]] <-  tmp
	}

}


pdf(paste0('./Plots/CASPER/Venn_Results_DE_celltype.pdf'))

cowplot::plot_grid(
	get_venn(de_genes_vsNon5q, FALSE)+ggtitle('non5q'),
	get_upset(de_genes_vsNon5q), nrow=2
)
uu <- get_upset(de_genes_vsElder)
upset_plot <- cowplot::plot_grid(NULL, uu$Main_bar, uu$Sizes, uu$Matrix,
                            nrow=2, align='hv', rel_heights = c(3,1),
                           rel_widths = c(2,3))
cowplot::plot_grid(
	get_venn(de_genes_vsElder, FALSE)+ggtitle('elder'),
	get_upset(de_genes_vsElder), nrow=2
)
dev.off()








all_seurat_integrated_sct$depleted_5q <- ifelse(colnames(all_seurat_integrated_sct) %in% all_5q_depleted_cells, '5q', 'mds')
Idents(all_seurat_integrated_sct) <- 'depleted_5q'
 
marker_5q <- FindMarkers(all_seurat_integrated_sct, ident.1 = '5q', ident.2='mds')


# pdf('./Test.pdf')
# UpSetR::upset(movies,attribute.plots=list(gridrows=60,plots=list(list(plot=UpSetR::scatter_plot, x="ReleaseDate", y="AvgRating"),
# list(plot=UpSetR::scatter_plot, x="ReleaseDate", y="Watches"),list(plot=UpSetR::scatter_plot, x="Watches", y="AvgRating"),
# list(plot=UpSetR::histogram, x="ReleaseDate")), ncols = 2))
# dev.off()