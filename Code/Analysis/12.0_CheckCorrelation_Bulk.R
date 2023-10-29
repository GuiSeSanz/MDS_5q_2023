
library(Seurat)
library(ggplot2)

normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

load('./Data/5Q-ANALYSIS.RData')
data_5q    <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

# select common genes
common_genes <- intersect(rownames(data_5q),rownames(countData))

sc_list <- SplitObject(data_5q, split.by = "Sample")

sc_list <- lapply(X = sc_list, FUN = function(x){
					x <- as.data.frame(x@assays$RNA@counts)
					x <- x[common_genes, ]
})

countData <- countData[common_genes,]


names <- colnames(countData)
names <- names[stringr::str_detect(names, '^SMD')]

cluster_cells <- FetchData(data_5q, vars= 'Cluster_names')
unique(stringr::str_extract(names, '^SMD[0-9]+'))
unique(stringr::str_extract(colnames(sc_counts), '^SMD[0-9]+'))



CELLTYPES <- intersect(unique(data_5q$Cluster_names), unique(stringr::str_extract(names, '[A-Z]+$')))




for (ct in CELLTYPES){
	message(ct)
	samples_2_keep_sc <- rownames(cluster_cells[cluster_cells$Cluster_names == ct,, drop=FALSE])
	for (smp in names(sc_list)){
		sc_tmp <- sc_list[[smp]]
		sc_tmp <- sc_tmp[,colnames(sc_tmp) %in% samples_2_keep_sc]
		sc_tmp <- as.data.frame(rowSums(sc_tmp))
		bulk_samples <- grep(ct, names, value=TRUE)
		countData_tmp <- countData[, bulk_samples]
		plot_list <- list()
		for (i in 1:ncol(countData_tmp)){
			x <- normalize_TPM(countData_tmp[, i])
			y <- normalize_TPM(sc_tmp[, 1])
			corr <- round(cor(x, y, method = "spearman"), 2)
			plot_list[[paste0(i, smp)]] <- ggplot(data = data.frame(x = log10(x), y = log10(y)), aes(x = x, y = y), alpha=0.8) +
			geom_point() + theme_classic() +
			geom_abline(intercept = 0, slope = 1, color = "red") +
			labs(x = "log10 Bulk TPM", y = "log10 Single cell TPM",
				title = paste0(colnames(countData_tmp)[i], " Vs " , smp," ", corr)) + 
			theme(plot.title = element_text(size=5))

		}
		pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/Correlation_Bulk/Correlation_', smp, '_', ct, '.pdf'))
			print(cowplot::plot_grid(plotlist = plot_list, ncol = 3))
		dev.off()
	}
}




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


viridis::scale_color_viridis(discrete=TRUE) 



for (ct in CELLTYPES){
	message(ct)
	samples_2_keep_sc <- rownames(cluster_cells[cluster_cells$Cluster_names == ct,, drop=FALSE])
	for (smp in names(sc_list)){
		sc_tmp <- sc_list[[smp]]
		sc_tmp <- sc_tmp[,colnames(sc_tmp) %in% samples_2_keep_sc]
		sc_tmp <- as.data.frame(rowSums(sc_tmp))
		bulk_samples <- grep(ct, names, value=TRUE)
		countData_tmp <- countData[, bulk_samples]
		plot_list <- list()
		for (i in 1:ncol(countData_tmp)){
			x <- normalize_TPM(countData_tmp[, i, drop=FALSE])
			y <- normalize_TPM(sc_tmp[, 1, drop=FALSE])
			plotter <- setNames(data.frame(x = log10(x), y = log10(y)), c('bulk', 'sc'))
			plotter$Pos <- ifelse(rownames(plotter) %in% genes_5q, '5q',
									ifelse( rownames(plotter) %in% genes_chr5, 'chr5', 'other'))
			corr <- round(cor(x, y, method = "spearman"), 2)
			plot_list[[paste0(i, smp)]] <- ggplot(data = plotter[plotter$Pos != 'other', ], aes(x = bulk, y = sc, color = Pos), alpha=0.8) +
			geom_point() + theme_classic() + viridis::scale_color_viridis(discrete=TRUE) + ggtitle('Control') +
			geom_abline(intercept = 0, slope = 1, color = "red") +
			labs(x = "log10 Bulk TPM", y = "log10 Single cell TPM",
				title = paste0(colnames(countData_tmp)[i], " Vs " , smp," ", corr)) + 
			theme(plot.title = element_text(size=5))

		}
		pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/Correlation_Bulk/Correlation_', smp, '_', ct, '_color_chr.pdf'))
			print(cowplot::plot_grid(plotlist = plot_list, ncol = 3))
		dev.off()
	}
}



library(Seurat)
library(ggplot2)

normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

# load('./Data/5Q-ANALYSIS.RData')
data_5q    <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_elder <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')



Create_pseudoBulk_by_Sample_Celltype <- function(data, split, cell_name){
	message(paste0('Splitting Object per: ', split))
	obj_lits <- SplitObject(data, split.by=split)
	message(paste0('Creating pseudobulk by: ', split))
	obj_lits_pseudoBulk <- lapply(X=obj_lits, FUN=function(x){
		smp_name <- unique(x$Sample)
		message(smp_name)
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudoBulk <- NULL
		for(ct in sort(unique(x[[cell_name]][[1]]))){
			message(paste0('--', ct,'--'))
			samples_2_keep <- colnames(x)[which(x[[cell_name]] == ct)]
			if (length(samples_2_keep) > 2){
				tmp <- counts[, colnames(counts) %in% samples_2_keep]
				tmp <- setNames(as.data.frame(rowSums(tmp)), paste0(smp_name, '_', ct))
				pseudoBulk[[paste0(smp_name, '_', ct)]]
				if (is.null(pseudoBulk)){
					pseudoBulk <- tmp
				}
				else{
					pseudoBulk <- cbind(pseudoBulk, tmp)
				}
			}else{message(paste0('Not enough cells to create pseudoBulk for this celltype: ', ct))}
		}
		return(pseudoBulk)
	})
	message('Yay!')
	return(obj_lits_pseudoBulk)
}


data_5q_list <-  SplitObject(data_5q, split.by='Sample')
data_5q_list <- lapply(X=data_5q_list, FUN=function(x){
	counts <- as.data.frame(x@assays$RNA@counts)
	pseudo <- as.data.frame(rowSums(counts))
	return(pseudo)
})

data_elder_list <-  SplitObject(data_elder, split.by='Sample')



pseudobulk_elder <- dplyr::bind_cols(elder_integrated_list_pseudobulk)
pseudobulk_5q <- dplyr::bind_cols(data_5q_annotated_list_pseudobulk)

pseudobulk_elder <- apply(pseudobulk_elder, 2, normalize_TPM)
pseudobulk_5q <- apply(pseudobulk_5q, 2, normalize_TPM)








normalize_TPM <- function(x){
	x <- x/sum(x)*1e6
	return(x)
}

data_5q    <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_elder <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

genes_2_check <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/genes_in_the_deleted_region/13-33-Table_1.tsv', sep='\t', header=T)
genes_2_check <- genes_2_check$Gene.name
genes_2_check <- genes_2_check[genes_2_check!= '']

data_5q_list <- SplitObject(data_5q, split.by='Sample')
data_5q_list <- lapply(data_5q_list, FUN=function(x){
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudo <- setNames(as.data.frame(rowSums(counts)), unique(x$Sample))
		return(pseudo)
})

data_elder_list <- SplitObject(data_elder, split.by='Sample')
data_elder_list <- lapply(data_elder_list, FUN=function(x){
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudo <- setNames(as.data.frame(rowSums(counts)), unique(x$Sample))
		return(pseudo)
})
data_5q_list <-  dplyr::bind_cols(data_5q_list)
data_elder_list <- dplyr::bind_cols(data_elder_list)

pdf('./Plots/Test.pdf')
# TMP normalization
data_5q_list_tpm <- apply(data_5q_list, 2, normalize_TPM)
data_elder_list_tpm <- apply(data_elder_list, 2, normalize_TPM)

data_5q_list_tpm <- data_5q_list_tpm[rownames(data_5q_list_tpm) %in% genes_2_check, ]
data_elder_list_tpm <- data_elder_list_tpm[rownames(data_elder_list_tpm) %in% genes_2_check, ]

data_5q_list_tpm <- colSums(data_5q_list_tpm)
data_elder_list_tpm <- colSums(data_elder_list_tpm)


plotter2 <- rbind(reshape2::melt(data_5q_list_tpm), reshape2::melt(data_elder_list_tpm))
plotter2$genotype <- ifelse(stringr::str_detect(rownames(plotter2), '^SMD+'), 'MDS', 'Control')
plotter2$Sample <- rownames(plotter2)

ggplot(plotter2, aes(x=genotype, y=value, fill=genotype)) +
	labs(y ='Sum of counts in region', title= 'TPM')+
	geom_boxplot() + geom_point() + 
	ggrepel::geom_text_repel(aes(label=Sample), size=3) +
	scale_fill_manual(values=c(MDS='#E69F00',  Control='#56B4E9')) +
	theme_classic() +
	theme(axis.title.x=element_blank())
dev.off()




data_5q_list <- SplitObject(data_5q, split.by='Sample')

pseudoBulk <- NULL
for (sample in names(data_5q_list)){
	message(sample)
	data_5q <- data_5q_list[[sample]]
	tmp_counts <- as.data.frame(data_5q@assays$RNA@counts)
	is_5q <- data_5q[['cell_5q']]
	is_5q <- rownames(is_5q[is_5q$cell_5q=='del5q', , drop=FALSE])
	tmp_counts_del5q  <- tmp_counts[, colnames(tmp_counts) %in% is_5q]
	tmp_counts_normal <- tmp_counts[, !colnames(tmp_counts) %in% is_5q]
	tmp_counts_del5q  <- setNames(as.data.frame(rowSums(tmp_counts_del5q)), 
								paste0(unique(data_5q$Sample), '_del5q'))
	tmp_counts_normal <- setNames(as.data.frame(rowSums(tmp_counts_normal)), 
								paste0(unique(data_5q$Sample), '_normal'))
	tmp_counts_del5q  <- normalize_TPM(tmp_counts_del5q)
	tmp_counts_normal <- normalize_TPM(tmp_counts_normal)
	if(is.null(pseudoBulk)){
		pseudoBulk <- cbind(tmp_counts_del5q, tmp_counts_normal)
	}
	else{
		pseudoBulk <- cbind(pseudoBulk, tmp_counts_del5q, tmp_counts_normal)
	}
}

genes_2_check2 <- c('CD74','BTF3','COX7C','HINT1','RPS23','RPS14')
pseudoBulk <- pseudoBulk[rownames(pseudoBulk) %in% genes_2_check2, ]
pseudoBulk <- colSumxxs(pseudoBulk)
pseudoBulk <- reshape2::melt(pseudoBulk)
pseudoBulk$Sample <- stringr::str_extract(rownames(pseudoBulk), '^[\\w\\d]+(?=_)')
pseudoBulk$genotype <- ifelse(stringr::str_detect(rownames(pseudoBulk), '_del5q'), 'del5q', 'normal')

pdf('./Plots/Pseudobulk_6_genes_delVsNon.pdf')
# ggplot(pseudoBulk, aes(x=genotype, y=value, fill=genotype)) +
# 	labs(y ='Sum of counts in region', title= 'TPM')+
# 	geom_boxplot() + geom_point() + 
# 	ggrepel::geom_text_repel(aes(label=Sample), size=3) +
# 	scale_fill_manual(values=c(del5q='#E69F00',  normal='#56B4E9')) +
# 	theme_classic() +
# 	theme(axis.title.x=element_blank())

	ggplot(pseudoBulk, aes(x=genotype, y=value, fill=Sample, group=Sample)) +
	labs(y ='Sum of counts in region', title= 'TPM')+
	geom_point() + geom_line() + 
	ggrepel::geom_text_repel(aes(label=Sample), size=3) +
	scale_fill_manual(values=c(del5q='#E69F00',  normal='#56B4E9')) +
	theme_classic() +
	theme(axis.title.x=element_blank())
dev.off()








data_5q_list <- lapply(data_5q_list, FUN=function(x){
		counts <- as.data.frame(x@assays$RNA@counts)
		pseudo <- setNames(as.data.frame(rowSums(counts)), unique(x$Sample))
		return(pseudo)
})








pseudoBulk_5q    <- data_5q_list[rownames(data_5q_list) %in% genes_2_check, ]
pseudoBulk_elder <-  data_elder_list[rownames(data_elder_list) %in% genes_2_check, ]
plotter_HM <- cbind(pseudoBulk_5q[rownames(pseudoBulk_5q) %in% rownames(pseudoBulk_elder), ], 
					pseudoBulk_elder[rownames(pseudoBulk_elder) %in% rownames(pseudoBulk_5q), ])

annotation_col <- data.frame(Genotype=ifelse(stringr::str_detect(colnames(plotter_HM), '^SMD+'), 'MDS', 'Control'))
rownames(annotation_col) <- colnames(plotter_HM)
annotation_color <- list(Genotype=c(MDS='#E69F00', Control='#56B4E9'))

pdf('./Plots/Test.pdf')
pheatmap::pheatmap(plotter_HM, annotation_col= annotation_col, annotation_colors=annotation_color, fontsize=3, main='PseudoBulks of 5q samples and elders') 
dev.off()


table(data_5q@assays$integrated@var.features %in% genes_2_check)
data_5q_scaled_data <- as.data.frame(data_5q@assays$SCT@scale.data)
data_5q_data <- as.data.frame(data_5q@assays$SCT@data)

plooter_hm <- data_5q_data[rownames(data_5q_data) %in% genes_2_check,]
plooter_hm_scale <- data_5q_scaled_data[rownames(data_5q_scaled_data) %in% genes_2_check,]
pdf('./Plots/Test.pdf')
pheatmap::pheatmap(plooter_hm[,sample(ncol(plooter_hm_scale), 1e3)],  scale='none', 
	show_rownames=TRUE, fontsize= 2,show_colnames=FALSE, color=viridis::viridis(50),
	cell_width=5, cluster_cols = TRUE)
dev.off()

