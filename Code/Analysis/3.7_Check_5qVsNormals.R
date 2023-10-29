

library(Seurat)
library(ggplot2)
library(argparse)


lib_size_normalization <-  function(counts){
	# remove the genes with all zeroes
	counts <- counts[, colSums(counts)>0]
	# normalize the counts
	libsize <- colSums(counts)
	size_factors <- libsize/mean(libsize)
	norm_counts <- counts/size_factors
	normalized_counts <- log(norm_counts+1)
	return(normalized_counts)
}

get_data <- function(sample_5q_norm_name, data_normal, data_5q = seurat_per_sample){
	data_5q_norm <- data_5q[[sample_5q_norm_name]]
	data_normal_norm <- data_normal[,common_genes]
	data_normal_norm <- data_normal[,common_genes]
	plotter <- as.data.frame(colMeans(data_normal))
	data_5q_norm <- as.data.frame(colMeans(data_5q_norm))
	plotter <- merge(plotter, data_5q_norm, by= 0)
	plotter <- setNames(plotter, c('gene_name', 'Normal', 'sample_5q'))
	plotter$Sample <- sample_5q_norm_name
	plotter$Pos <- ifelse(plotter$gene_name %in% genes_5q, '5q',
						  ifelse( plotter$gene_name %in% genes_chr5, 'chr5', 'other'))
	return(plotter)
}

plot_correlation <- function(plotter, title){
	wt <- wilcox.test(plotter$Normal, plotter$sample_5q)
	p <- ggplot() +
		geom_point(plotter[plotter$Pos == 'other', ], mapping = aes(x = sample_5q, y = Normal, color = Pos), alpha=0.3) + 
		geom_point(plotter[plotter$Pos == 'chr5', ], mapping = aes(x = sample_5q, y = Normal, color = Pos), alpha=0.6) + 
		geom_point(plotter[plotter$Pos == '5q', ], mapping = aes(x = sample_5q, y = Normal, color = Pos), alpha=0.9) + 
		theme_classic() + viridis::scale_color_viridis(discrete=TRUE) + 
		labs(title=title, subtitle = paste0('  Wilcox test W:',formatC(wt$statistic, format = "e", digits = 2), '  p.val: ',formatC(wt$p.val, format = "e", digits = 2))) + 
		geom_abline(intercept = 0, slope = 1) + 
		theme(legend.position = "bottom", axis.text.x = element_blank(), axis.text.y = element_blank(), 
			   axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
	return(p)
}

get_umap_signature <- function(df, title){
	plot <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color= Signature)) + 
	geom_point(alpha=0.9, size = 0.8) + 
	# scale_color_gradientn(low="grey90", high ="blue", name = 'Expression') + 
	# guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
	ggprism::theme_prism() + ggtitle(title)+  labs(x = "UMAP 1", y = 'UMAP 2') +
	theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 12),  axis.text=element_blank(), axis.ticks=element_blank())#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
	if(is.numeric(df$Signature)){
		  plot <- plot + viridis::scale_color_viridis()
	}else{
	  plot <- plot + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))
	}
	return(plot)
}

get_expression_signature <- function(gene_list, df, coords = coords, upplim=NULL, split_HL=FALSE, noPlot=FALSE){
		check_genes <- function(gene_list, df){
			not_found_genes <- gene_list[which(!gene_list %in% colnames(df))]
			if(length(not_found_genes)!=0){
				message(paste0('Left behind ', length(not_found_genes), ': ', paste0(not_found_genes, collapse='; ')))
				gene_list <- gene_list[which(gene_list %in% colnames(df))]
			}else{
				gene_list <- gene_list
			}

			return(gene_list)
		}
	gene_list <- check_genes(gene_list, df)
	tmp <- df[, colnames(df) %in% gene_list]
	tmp$Signature <- rowMeans(tmp)
	# tmp$cell_id <- sub('-', '\\.',rownames(tmp))
	tmp <- merge(tmp, coords, by=0)
	if(!is.null(upplim)){
		upp_outliers <- which(tmp$Signature>upplim) 
		tmp[upp_outliers , 'Signature'] <- upplim
	}
	plot <- get_umap_signature(tmp, 'title')
	if(split_HL){
		plot <- plot + facet_wrap(~BinScore)
	}
	if(noPlot){
		return(tmp)
	}else{
		return(plot)
	}
}

get_top_bottom <- function(){
	normData <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$SCT@data)))
	coords_mds <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
	coords_mds$Sample <- stringr::str_extract(rownames(coords_mds), '^[\\w]+(?=_)')
	signature <-  get_expression_signature(unique(genes_5q), normData, coords = coords_mds,noPlot=TRUE)
	sub_signature <- split(signature, signature$Sample)

	top_list <- list()
	bottom_list <- list()
	for (sample in names(sub_signature)){
		tmp <- sub_signature[[sample]]
		tmp <- tmp[, c('Row.names', 'Signature')]
		tmp <- tmp[order(-tmp$Signature),]
		ten_percent <- floor(nrow(tmp)/10)
		top <- tmp[1:ten_percent,]
		a <- nrow(tmp) - ten_percent
		bottom <- tmp[a:nrow(tmp),]
		top_list[[sample]] <- top
		bottom_list[[sample]] <- bottom
	}
	return(list(Top=top_list, Bottom=bottom_list))
}

####### Annotations #########
	# gene_list <-  read.table('./Data/Annotation/5q_signature_5q13-33.txt', fill=TRUE, header=TRUE)
	# genes_5q <- c()
	# for (name in names(gene_list)){
	# 	genes_5q <- c(genes_5q, as.character(levels(gene_list[[name]])))
	# }
	genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
	genes_5q <- as.character(unique(genes_expanded$V5))

	ann <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/Homo_sapiens.GRCh38.105_CHR5.gtf', sep='\t')
	ann <- ann[, c('V1', 'V3', 'V4', 'V5', 'V9'),]
	ann <- ann[ann$V3 == 'gene',]
	ann$gene_id <- stringr::str_extract(ann$V9,'(?<=gene_id\\s)[\\w]+(?=;)')
	ann$gene_name <- stringr::str_extract(ann$V9,'(?<=gene_name\\s)[\\w]+(?=;)')
	ann <- ann[, c('V1', 'V4', 'V5', 'gene_id', 'gene_name'),]
	ann <- ann[!is.na(ann$gene_name),]
	genes_chr5 <- unique(ann$gene_name)
#############################





parser <- ArgumentParser(description='Create plots for the normals from young or senior samples')
parser$add_argument('-S', '--Sample', default="young", type="character", help='young or senior')
args <- parser$parse_args()
SAMPLE <- toupper(args$Sample)


print(paste0('Sample: ', SAMPLE))

all_seurat_integrated_sct <- readRDS('/mnt/md0/gserranos/MDS/Data/all_seurat_integrated_sct.rds')
# Idents(all_seurat_integrated_sct) <- 'Sample'
seurat_per_sample <- SplitObject(all_seurat_integrated_sct, split.by = "Sample")

for (sample in names(seurat_per_sample)){
	counts <- t(as.data.frame(seurat_per_sample[[sample]]@assays$RNA@counts))
	sample_norm <- lib_size_normalization(counts)
	seurat_per_sample[[sample]] <- sample_norm
}

data_path <- '/home/tereshkova/data/gserranos/MDS/Data/Normal_Data/RawData/'
if (SAMPLE == "YOUNG"){
	normal_list <-c( 'GSM5460406_young1_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460407_young2_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460408_young3_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460409_young4_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460410_young5_filtered_feature_bc_matrix_counts.tsv')
}else if(SAMPLE == "SENIOR"){
	normal_list <-c('GSM5460411_elderly1_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460412_elderly2_filtered_feature_bc_matrix_counts.tsv',
					'GSM5460413_elderly3_filtered_feature_bc_matrix_counts.tsv')
}


normal_list2 <- normal_list
normal_plot_list <- list()
for (normal_sample1 in normal_list){
	message(normal_sample1)
	normal_name_1 <- stringr::str_extract(normal_sample1,'(?<=_)[a-z]+[\\d]{1}')
	data_norm_1 <- read.table(paste0(data_path, normal_sample1), sep='\t', header=TRUE)
	rownames(data_norm_1) <- data_norm_1$X
	data_norm_1 <- data_norm_1[,-1]
	data_norm_1 <- lib_size_normalization(data_norm_1)
	plot_list <- list()
	common_genes <- intersect(colnames(data_norm_1), unique(unlist(lapply(seurat_per_sample, colnames))))

	for (sample in names(seurat_per_sample)){
		plot_list[[sample]] <- plot_correlation(get_data(sample, data_norm_1), sample) + theme(legend.position='none')
	}
	# Normal Vs 5q
	legend <- cowplot::get_legend( plot_correlation(get_data(sample, data_norm_1), sample))
	pdf(paste0('./Plots/Integrated/5q_vs_',normal_name_1,'_', SAMPLE,'_mean.pdf'))
	print(cowplot::plot_grid(
		cowplot::plot_grid(
			plotlist = plot_list, 
			ncol = 2
		),
		legend,
		nrow=2, rel_heights=c(1, 0.1)))
	dev.off()
	# Normal Vs Normal
	normal_list2 <- normal_list2[!normal_list2 %in% normal_sample1]
	for (normal_sample2 in normal_list2){
		if(normal_sample1 != normal_sample2){
			print(paste0(normal_sample1, '_Vs_',  normal_sample2))

			normal_name_2 <- stringr::str_extract(normal_sample2,'(?<=_)[a-z]+[\\d]{1}')
			young2 <- read.table(paste0(data_path, normal_sample2), sep='\t', header=TRUE)
			rownames(young2) <- young2$X 
			young2 <- young2[,-1]
			data_norm_2 <- lib_size_normalization(young2)

			data_norm_2_plot <- as.data.frame(colMeans(data_norm_2))
			young_norm_plot <- as.data.frame(colMeans(data_norm_1))
			gene_2_keep <- intersect(rownames(data_norm_2_plot), rownames(young_norm_plot))

			data_norm_2_plot <- data_norm_2_plot[rownames(data_norm_2_plot) %in% gene_2_keep,, drop=FALSE]
			data_norm_2_plot <- data_norm_2_plot[rownames(data_norm_2_plot) %in% gene_2_keep,, drop=FALSE]
			control_plotter <- merge(young_norm_plot, data_norm_2_plot, by= 0)
			control_plotter <- setNames(control_plotter, c('gene_name', 'Normal', 'Normal_2'))
			control_plotter$Pos <- ifelse(control_plotter$gene_name %in% genes_5q, '5q',
									ifelse( control_plotter$gene_name %in% genes_chr5, 'chr5', 'other'))
			
			p <- ggplot() +
				geom_point(control_plotter[control_plotter$Pos == 'other', ], mapping = aes(x = Normal_2, y = Normal, color = Pos), alpha=0.3) + 
				geom_point(control_plotter[control_plotter$Pos == 'chr5', ] , mapping = aes(x = Normal_2, y = Normal, color = Pos), alpha=0.6) + 
				geom_point(control_plotter[control_plotter$Pos == '5q', ]   , mapping = aes(x = Normal_2, y = Normal, color = Pos), alpha=0.9) + 
				theme_classic() + viridis::scale_color_viridis(discrete=TRUE) + ggtitle('Control') + 
				geom_abline(intercept = 0, slope = 1) +labs(title="Control", x =normal_name_2, y = normal_name_1) +
				theme(legend.position = "bottom", axis.text.x = element_blank(), axis.text.y = element_blank(), 
						axis.ticks.x = element_blank(), axis.ticks.y = element_blank())
			normal_plot_list[[paste0(normal_name_1, normal_name_2)]] <- p
			
			
		}
	}
}


# test
pdf(paste0('./Plots/Integrated/normal_vs_SMD_mean.pdf'))
print(plot_list[1])
dev.off()

normal_plot_list_2 <- list() 
for (name in names(normal_plot_list)){
	normal_plot_list_2[[name]] <- normal_plot_list[[name]] + theme(legend.position='none') + labs(title="")
}

pdf(paste0('./Plots/Integrated/normal_vs_normal_',SAMPLE,'_mean.pdf'))
cowplot::plot_grid(plotlist = normal_plot_list_2 , nrow=3)
dev.off()

top_bottom <- get_top_bottom()

for (normal_sample1 in normal_list){
	print(normal_sample1)
	normal_name_1 <- stringr::str_extract(normal_sample1,'(?<=_)[a-z]+[\\d]{1}')
	data_norm_1 <- read.table(paste0(data_path, normal_sample1), sep='\t', header=TRUE)
	rownames(data_norm_1) <- data_norm_1$X
	data_norm_1 <- data_norm_1[,-1]
	data_norm_1 <- lib_size_normalization(data_norm_1)
	plot_list_TB <- list()
	common_genes <- intersect(colnames(data_norm_1), unique(unlist(lapply(seurat_per_sample, colnames))))

	pdf(paste0('./Plots/Integrated/',normal_name_1,'_top_bottom.pdf'))
	for (sample in names(seurat_per_sample)){
		print(sample)
		top <- top_bottom$Top[[sample]]
		bottom <- top_bottom$Bottom[[sample]]
		data <- seurat_per_sample[[sample]]
		data_top <- data[top$Row.names,colnames(data) %in% common_genes]
		data_bottom <- data[bottom$Row.names,colnames(data) %in% common_genes]
		data_top <- as.data.frame(colMeans(data_top))
		data_bottom <- as.data.frame(colMeans(data_bottom))
		normal <- data_norm_1[,colnames(data_norm_1) %in% common_genes]
		normal <- as.data.frame(colMeans(normal))
		plotter_top     <- setNames(merge(normal, data_top, by= 0), c('gene_name', 'Normal', 'sample_5q'))
		plotter_bottom	<- setNames(merge(normal, data_bottom, by= 0), c('gene_name', 'Normal', 'sample_5q'))

		plotter_top$Pos <- ifelse(plotter_top$gene_name %in% genes_5q, '5q',
						  ifelse( plotter_top$gene_name %in% genes_chr5, 'chr5', 'other'))
		plotter_bottom$Pos <- ifelse(plotter_bottom$gene_name %in% genes_5q, '5q',
						  ifelse( plotter_bottom$gene_name %in% genes_chr5, 'chr5', 'other'))
		
		print(cowplot::plot_grid(plot_correlation(plotter_top, paste0(sample, '_top')),
							plot_correlation(plotter_bottom, paste0(sample, '_bottom')),
							ncol=2))
		
	}
	dev.off()
}




