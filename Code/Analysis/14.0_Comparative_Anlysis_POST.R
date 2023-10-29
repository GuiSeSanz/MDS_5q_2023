
library(Seurat)
library(Libra)
library(ggplot2)

get_volcano <- function(results, CType, sub= "5q Vs Elder"){
	p <-  EnhancedVolcano::EnhancedVolcano(results,
				lab = results$gene_name,
				x = 'avg_logFC', captionLabSize=7,
				y = 'p_val_adj',col=c('black', 'black', 'black', 'red3'),
				subtitle =sub, axisLabSize = 5, title= CType,
				pCutoff = 0.05, labSize=2,
				FCcutoff = 0.2) +
				theme(legend.position='none', text = element_text(family = "Helvetica", size = 2),
				axis.title= element_text(family = "Helvetica", size = 7))
	return(p)
}

get_tornado <- function(results, title){
	p <- ggplot(results, aes(x = NES, y = reorder(pathway, NES), fill = is_pos)) + 
			geom_bar(stat = "identity", position = "identity") + theme_classic() +
			ylab('Pathways in C2')+
			theme(legend.position="none", text = element_text(family = "Helvetica", size = 7), axis.text.y= element_text(size = 5),) + 
			scale_fill_manual(values=c(Neg='#006e90', Pos='#cc2936')) + ggtitle(title)
	return(p)
}

get_prop_tables_per_sample <- function(Sobj){
	print(
		prop.table(table(Sobj$Sample, Sobj$cell_5q), margin=1)*100
	)
}

get_HM_TopBottom <- function(DE_data, data_sc, meta, name){
	all_results_psuedo <- DE_data 
	pdf(paste0('./Plots/DE_POST/Heatmap_TopBottom_',name, '.pdf'))
	for (ct in sort(names(all_results_psuedo))){
		message(ct)
		tmp <- all_results_psuedo[[ct]]
		tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
		gene_names <- c(tmp[1:10,]$gene, tmp[(nrow(tmp)-9):nrow(tmp),]$gene)
		meta_tmp <- meta[meta$cell_type == ct,]
		contrast_stats <- as.data.frame(table(meta_tmp$label))
		if(any(contrast_stats[,2]<20)){
			message('Not enough samples for ', ct)
			next
		}
		meta_tmp <- split(meta_tmp, meta_tmp$label)
		sampled_rows <- lapply(meta_tmp, function(df) sample(rownames(df), 20))
		annotation_column <- data.frame(samples=sampled_rows[[1]], label=names(sampled_rows)[[1]])
		annotation_column <- rbind(annotation_column, 
							 data.frame(samples=sampled_rows[[2]], label=names(sampled_rows)[[2]]))
		tmp <- data_sc[gene_names, annotation_column$samples]
		rownames(annotation_column) <- annotation_column$samples
		annotation_column <- annotation_column[,2, drop=FALSE]
		anno_colors <- list(label = c('#006e90', '#cc2936'))
		names(anno_colors[['label']]) <- unique(annotation_column$label)
		# pdf('./Plots/DE_POST/Heatmap_TopBottom_TEST.pdf')
		pheatmap::pheatmap(tmp, annotation_col = annotation_column,
					annotation_colors = anno_colors, 
					show_rownames = TRUE, show_colnames = FALSE, 
					cluster_rows = FALSE, cluster_cols = FALSE, 
					scale = 'row', fontsize_row = 5, fontsize_col = 5,
					cellwidth = 5, cellheight = 5, border_color = NA,
					main = ct)
	}
	dev.off()
}

# get_HM_TopBottom_MAST(results, contrast_data, ct )
get_HM_TopBottom_MAST <- function(de_results, data, ct ){
	data_tmp <- as.data.frame(data@assays$RNA@counts)
	cells <- FetchData(data, vars=c('Status', 'Cluster_names'))
	de_results <- de_results[order(de_results$avg_log2FC), ]
	de_results$gene_name <- rownames(de_results)
	gene_names <- c(de_results[1:10,]$gene_name, 
					de_results[(nrow(de_results)-9):nrow(de_results),]$gene_name)
	cells_tmp <- cells[cells$Cluster_names == ct & cells$Status != 'None',]
	cells_tmp <- split(cells_tmp, cells_tmp$Status)
	sampled_rows <- lapply(cells_tmp, function(df) {
							if (nrow(df) < 20) {
								sample(rownames(df), nrow(df))
							} else {
								sample(rownames(df), 20)
							}
							})
	annotation_column <- data.frame(samples=sampled_rows[[1]], label=names(sampled_rows)[[1]])
	annotation_column <- rbind(annotation_column, 
							data.frame(samples=sampled_rows[[2]], label=names(sampled_rows)[[2]]))
	tmp <- data_tmp[gene_names,annotation_column$samples]
	rownames(annotation_column) <- annotation_column$samples
	annotation_column <- annotation_column[,2, drop=FALSE]
	anno_colors <- list(label = c('#006e90', '#cc2936'))
	names(anno_colors[['label']]) <- unique(annotation_column$label)
	# pdf('./Plots/DE_POST/Heatmap_TopBottom_TEST.pdf')
	p <- pheatmap::pheatmap(tmp, annotation_col = annotation_column,
				annotation_colors = anno_colors, 
				show_rownames = TRUE, show_colnames = FALSE, 
				cluster_rows = FALSE, cluster_cols = FALSE, 
				scale = 'row', fontsize_row = 5, fontsize_col = 5,
				cellwidth = 5, cellheight = 5, border_color = NA,
				main = ct, silent=TRUE)
	return(p)
}
 

saveDE_results <- function(results, name){
		
	results <- as.data.frame(results)
	DE_flt <- results[results$p_val_adj < 0.05,]
	all_results_psuedo <- split( DE_flt , f = DE_flt$cell_type )
	saveRDS(all_results_psuedo ,paste0('./Data/Results_DE_',name, '.rds'))
	pdf(paste0('./Plots/DE_POST/Volcano_And_tornado_',name, '.pdf'))
	for (cell_type in sort(names(all_results_psuedo))){
		print(cell_type)
		tmp <- all_results_psuedo[[cell_type]]
		tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
		rownames(tt) <- tt$gene_name
		print(get_volcano(tt, cell_type, sub=name))
		# tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
		# ranks <- setNames(tmp$avg_logFC, tmp$gene)
		# results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
		# results_c2 <- results_c2[results_c2$padj < 0.05,]
		# results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
		# results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
		# results_c2$pathway <- factor(results_c2$pathway)
		# if(nrow(results_c2)>0){
		# print(get_tornado(results_c2, cell_type))
		# }
	}
	dev.off()
	WriteXLS::WriteXLS(all_results_psuedo, ExcelFileName=paste0('./Plots/DE_POST/Results_DE_', name,'.xlsx'), 
	SheetNames = names(all_results_psuedo),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

}

get_cells_by_condition_by_cluster <- function(obj){
	stat_list <- list()
	meta <- FetchData(obj, vars = c('cell_5q', 'Cluster_names', 'Sample'))
	for (smp in unique(meta$Sample)){
		tmp <- meta[meta$Sample == smp,]
		print(smp)
		print(table(tmp$Cluster_names, tmp$cell_5q))
		stat_list[[smp]] <- as.data.frame(table(tmp$Cluster_names, tmp$cell_5q))
	}
	return(stat_list)
}


pathways_c2 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c2.cp.v7.5.1.symbols.gmt')

post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
all_5q_depleted_cells_COPYKAT <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')
#  There are less cells as some clusters have disapeared based on the annotation

# elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')

elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

get_prop_tables_per_sample(post_data)
get_prop_tables_per_sample(MDS_data)
get_prop_tables_per_sample(pre_post_data)

get_cells_by_condition_by_cluster(post_data)
get_cells_by_condition_by_cluster(MDS_data)
get_cells_by_condition_by_cluster(pre_post_data)

############################
# 1 Non del(5q) analysis
############################

# 1.1 Non del(5q) CR Vs non del(5q) MDS 

post_data_normal <- subset(x=post_data, subset = Sample == 'FS-0634-post')
post_data_normal <- subset(x = post_data_normal, subset = cell_5q == "normal")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'CR'

MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "normal")
expression_MDS <- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'Diagnosed'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_MDS))

all_data <- merge(expression_post[genes_2_keep,], expression_MDS[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_MDS)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
						   c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('CR', 'Diagnosed'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qCR_Vs_non5qDiagnosed')
saveDE_results(DE, 'Non5qCR_Vs_non5qDiagnosed')


# 1.2 Non del(5q) CR Vs non del(5q) PR
get_prop_tables_per_sample(post_data)
post_data_normal <- subset(x = post_data, subset = cell_5q == "normal")
post_data_normal <- PrepSCTFindMarkers(post_data_normal)
results_DE <- list()
cell_types <- sort(unique(post_data_normal$Cluster_names))
HM_LIST <- list()
pdf('./Plots/DE_POST/Results_DE_Non5qCR_Vs_non5qPR_MAST.pdf')
for(ct in cell_types){
	metadata <- FetchData(post_data_normal, vars = c('Sample', 'Cluster_names'))
	cell_1 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'FS-0634-post',])
	cell_2 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'FS-0406-post',])
	post_data_normal$Status <- ifelse(colnames(post_data_normal) %in% cell_1, 'CompleteResponder', 
									 ifelse(colnames(post_data_normal) %in% cell_2, 'PartialResponder', 'None'))
	Idents(post_data_normal)  <- 'Status'
	results <- FindMarkers(post_data_normal, 
							ident.1 = 'CompleteResponder',
							ident.2 = 'PartialResponder',
							assay='RNA',
							test.use = 'MAST', 
							latent.vars='Sample')
	results_DE[[ct]] <- results
	HM_LIST[[ct]] <- get_HM_TopBottom_MAST(results, post_data_normal, ct )

	tt <- results
	tt$gene_name <- rownames(tt)
	tt$avg_logFC <- tt$avg_log2FC
	tt <- tt[tt$p_val_adj <= 0.05,]
	print(get_volcano(tt, ct, sub= paste0("CompleteResponder Vs PartialResponder")))
	tt <- tt[order(tt$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tt$avg_logFC, tt$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
		print(get_tornado(results_c2, ct))
	}
}
dev.off()


pdf('./Plots/DE_POST/HM_DE_Non5qCR_Vs_non5qPR_MAST.pdf')
for (ct in cell_types){
	 grid::grid.draw(HM_LIST[[ct]]$gtable)
	 grid::grid.newpage()
}
dev.off()

results_DE <- lapply(results_DE, function(x) {
				x <- x[x$p_val_adj <= 0.05,]
				print(dim(x))
				return(x)
})

WriteXLS::WriteXLS(results_DE, ExcelFileName=paste0('./Plots/DE_POST/Results_DE_Non5qCR_Vs_non5qPR_MAST.xlsx'), SheetNames = names(results_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

# 1.3 nonDel(5q) CR Vs Elder
get_prop_tables_per_sample(post_data)

post_data_normal <- subset(x=post_data, subset = Sample == 'FS-0634-post')
post_data_normal <- subset(x = post_data_normal, subset = cell_5q == "normal")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'CR'

elder_data$Status <- 'Healthy'
expression_elder <- as.data.frame(elder_data@assays$RNA@counts)

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_elder))

all_data <- merge(expression_post[genes_2_keep,], expression_elder[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_elder)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   elder_data[[c('Status', 'Sample', 'Cluster_names')]]), 
					c('label', 'replicate', 'cell_type'))

table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

metadata$label <- factor(metadata$label, levels=c('CR', 'Healthy'))
DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qCR_Vs_Elder')
saveDE_results(DE, 'Non5qCR_Vs_Elder')

# 1.4 nonDel(5q) PR vs nonDel(5q) MDS
get_prop_tables_per_sample(MDS_data)
get_prop_tables_per_sample(post_data)

post_data_normal <- subset(x=post_data, subset = Sample == 'FS-0406-post')
post_data_normal <- subset(x = post_data_normal, subset = cell_5q == "normal")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'PR'

MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "normal")
expression_MDS <- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'Diagnosed'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_MDS))

all_data <- merge(expression_post[genes_2_keep,], expression_MDS[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_MDS)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
						   c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PR', 'Diagnosed'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qPR_Vs_non5qDiagnosed')

saveDE_results(DE, 'Non5qPR_Vs_non5qDiagnosed')

# 1.5 nonDel(5q) PR + CR vs nonDel(5q) MDS
get_prop_tables_per_sample(MDS_data)
get_prop_tables_per_sample(post_data)

post_data_normal <- subset(x = post_data, subset = cell_5q == "normal")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'PR_CR'

MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "normal")
expression_MDS <- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'Diagnosed'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_MDS))

all_data <- merge(expression_post[genes_2_keep,], expression_MDS[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_MDS)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
					c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PR_CR', 'Diagnosed'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qPR_CR_Vs_non5qDiagnosed')

saveDE_results(DE, 'Non5qPR_CR_Vs_non5qDiagnosed')


# 1.6 nonDel(5q) PR + CR vs Elderly
table(elder_data$Sample, elder_data$Cluster_names)
get_prop_tables_per_sample(post_data)

post_data_normal <- subset(x = post_data, subset = cell_5q == "normal")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'PR_CR'

expression_elder <- as.data.frame(elder_data@assays$RNA@counts)
elder_data$Status <- 'Elder'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_elder))

all_data <- merge(expression_post[genes_2_keep,], expression_elder[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_elder)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   elder_data[[c('Status', 'Sample', 'Cluster_names')]]), 
					c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PR_CR', 'Elder'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=2)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qPR_CRVsElder')

saveDE_results(DE, 'Non5qPR_CRVsElder')

# 1.7 nonDel(5q) PR vs Elderly
table(elder_data$Sample, elder_data$Cluster_names)
get_prop_tables_per_sample(post_data)

post_data_normal <- subset(x = post_data, subset = cell_5q == "normal")
post_data_normal <- subset(x = post_data_normal, subset = Sample == "FS-0406-post")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'PR'

expression_elder <- as.data.frame(elder_data@assays$RNA@counts)
elder_data$Status <- 'Elder'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_elder))

all_data <- merge(expression_post[genes_2_keep,], expression_elder[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_elder)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   elder_data[[c('Status', 'Sample', 'Cluster_names')]]), 
					c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PR', 'Elder'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qPR_VsElder')

saveDE_results(DE, 'Non5qPR_VsElder')


############################
# 2 Del(5q) analysis
############################
# 2.1 del(5q) NonResponder Vs del(5q) partial Responder

get_prop_tables_per_sample(post_data)
get_prop_tables_per_sample(pre_post_data)

post_data_del5q <- subset(x=post_data, subset = Sample == 'FS-0406-post')
post_data_del5q <- subset(x = post_data_del5q, subset = cell_5q == "del5q")
post_data_del5q$Status <- 'PR'

pre_post_data_del5q <- subset(x=pre_post_data, subset = Sample == 'SMD132114579')
pre_post_data_del5q <- subset(x = pre_post_data_del5q, subset = cell_5q == "del5q")
pre_post_data_del5q$Status <- 'NoR'

contrast_data <- merge(x=post_data_del5q, y=pre_post_data_del5q)
# contrast_data <- PrepSCTFindMarkers(contrast_data)


results_DE <- list()
cell_types <- sort(unique(contrast_data$Cluster_names))
metadata <- FetchData(contrast_data, vars = c('Sample', 'Cluster_names'))
HM_LIST <- list()
pdf('./Plots/DE_POST/Results_DE_Del5qPRVsDel5qNR_MAST.pdf')
for(ct in cell_types){
	message(ct)
	cell_1 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'FS-0406-post',])
	cell_2 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'SMD132114579',])
	if(length(cell_1) == 0 || length(cell_2) == 0){
		next
	}
	contrast_data$Status <- ifelse(colnames(contrast_data) %in% cell_1, 'PartialResponder', 
									 ifelse(colnames(contrast_data) %in% cell_2, 'NonResponder', 'None'))
	Idents(contrast_data)  <- 'Status'
	results <- FindMarkers(contrast_data, 
							ident.1 = 'PartialResponder',
							ident.2 = 'NonResponder',
							assay='RNA',
							test.use = 'MAST', 
							latent.vars='Sample')
	results_DE[[ct]] <- results
	HM_LIST[[ct]] <- get_HM_TopBottom_MAST(results, contrast_data, ct )
	tt <- results
	tt$gene_name <- rownames(tt)
	tt$avg_logFC <- tt$avg_log2FC
	tt <- tt[tt$p_val_adj <= 0.05,]
	print(get_volcano(tt, ct, sub= paste0("PartialResponder Vs NonResponder")))
	tt <- tt[order(tt$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tt$avg_logFC, tt$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
		print(get_tornado(results_c2, ct))
	}
}
dev.off()

results_DE <- lapply(results_DE, function(x) {
				x <- x[x$p_val_adj <= 0.05,]
				print(dim(x))
				return(x)
})

saveRDS(results_DE, './Data/Results_DE_Del5qPRVsDel5qNR_MAST.rds')

pdf('./Plots/DE_POST/HM_DE_Del5qPRVsDel5qNR_MAST.pdf')
for (ct in cell_types){
	 grid::grid.draw(HM_LIST[[ct]]$gtable)
	 grid::grid.newpage()
}
dev.off()

WriteXLS::WriteXLS(results_DE, ExcelFileName=paste0('./Plots/DE_POST/Results_DE_Del5qPRVsDel5qNR_MAST.xlsx'), 
SheetNames = names(results_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

# 2.2 del(5q) NonResponderPost Vs del(5q) NonResponderPre
get_prop_tables_per_sample(pre_post_data)
pre_post_data_del5q <- subset(x = pre_post_data, subset = cell_5q == "del5q")

results_DE <- list()
cell_types <- sort(unique(pre_post_data_del5q$Cluster_names))
metadata <- FetchData(pre_post_data_del5q, vars = c('Sample', 'Cluster_names'))
HM_LIST <- list()
pdf('./Plots/DE_POST/Results_DE_Del5qNR_preVsDel5qNR_post_MAST.pdf')
for(ct in cell_types){
	cell_1 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'SMD211420',])
	cell_2 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$Sample == 'SMD132114579',])
	if(length(cell_1) == 0 || length(cell_2) == 0){
		next
	}
	pre_post_data_del5q$Status <- ifelse(colnames(pre_post_data_del5q) %in% cell_1, 'Pre', 
									 ifelse(colnames(pre_post_data_del5q) %in% cell_2, 'Post', 'None'))
	Idents(pre_post_data_del5q)  <- 'Status'
	results <- FindMarkers(pre_post_data_del5q, 
								ident.1 = 'Pre',
								ident.2 = 'Post',
								assay='RNA',
								test.use = 'MAST', 
								latent.vars='Sample')
	results_DE[[ct]] <- results
	tt <- results
	tt$gene_name <- rownames(tt)
	tt$avg_logFC <- tt$avg_log2FC
	tt <- tt[tt$p_val_adj <= 0.05,]
	HM_LIST[[ct]] <- get_HM_TopBottom_MAST(results, pre_post_data_del5q, ct )

	print(get_volcano(tt, ct, sub= paste0("PartialResponderPre Vs PartialResponderPost")))
	tt <- tt[order(tt$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tt$avg_logFC, tt$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
		print(get_tornado(results_c2, ct))
	}
}
dev.off()

pdf('./Plots/DE_POST/HM_DE_Del5qNR_preVsDel5qNR_post_MAST.pdf')
for (ct in cell_types){
	 grid::grid.draw(HM_LIST[[ct]]$gtable)
	 grid::grid.newpage()
}
dev.off()

results_DE <- lapply(results_DE, function(x) {
				x <- x[x$p_val_adj <= 0.05,]
				print(dim(x))
				return(x)
})


WriteXLS::WriteXLS(results_DE, ExcelFileName=paste0('./Plots/DE_POST/Results_DE_Del5qNR_preVsDel5qNR_post_MAST.xlsx'), SheetNames = names(results_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

# 2.3 del(5q) partial Responder Vs del(5q) MDS

get_prop_tables_per_sample(post_data)
get_prop_tables_per_sample(MDS_data)

post_data_normal <- subset(x=post_data, subset = Sample == 'FS-0406-post')
post_data_normal <- subset(x = post_data_normal, subset = cell_5q == "del5q")
expression_post  <- as.data.frame(post_data_normal@assays$RNA@counts)
post_data_normal$Status <- 'PartialResponder'

MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "del5q")
expression_MDS <- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'Diagnosed'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_MDS))

all_data <- merge(expression_post[genes_2_keep,], expression_MDS[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_MDS)

metadata <- setNames(rbind(post_data_normal[[c('Status', 'Sample', 'Cluster_names')]], 
						   MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
							c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PartialResponder', 'Diagnosed'))

table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'del5qPR_Vs_del5qMDS')
saveDE_results(DE, 'del5qPR_Vs_del5qMDS')


# 2.4 Del(5q)NR_Post VS del(5q) MDS

get_prop_tables_per_sample(pre_post_data)
get_prop_tables_per_sample(MDS_data)

pre_post_data_del <- subset(x=pre_post_data, subset = Sample == 'SMD132114579')
pre_post_data_del <- subset(x = pre_post_data_del, subset = cell_5q == "del5q")
expression_pre_post  <- as.data.frame(pre_post_data_del@assays$RNA@counts)
pre_post_data_del$Status <- 'NR_post'

MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "del5q")
expression_MDS <- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'Diagnosed'

genes_2_keep <- intersect(rownames(expression_pre_post), rownames(expression_MDS))

all_data <- merge(expression_pre_post[genes_2_keep,], expression_MDS[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_pre_post) + ncol(expression_MDS)

metadata <- setNames(rbind(pre_post_data_del[[c('Status', 'Sample', 'Cluster_names')]], 
						   MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
						   c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('NR_post', 'Diagnosed'))

table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'del5qNR_Post_Vs_del5qMDS')
saveDE_results(DE, 'del5qNR_Post_Vs_del5qMDS')


############################
# 3 Mixed analysis
############################

# 3.1 del(5q) PR Vs non-del(5q) PR
get_prop_tables_per_sample(post_data)
post_data_pr <- subset(x=post_data, subset = Sample == 'FS-0406-post')

results_DE <- list()
cell_types <- sort(unique(post_data_pr$Cluster_names))
metadata <- FetchData(post_data_pr, vars = c('cell_5q', 'Cluster_names', 'Sample'))
HM_LIST <- list()
pdf('./Plots/DE_POST/Results_DE_Del5qPRVsNonDel5qPR_MAST.pdf')
for(ct in cell_types){
	cell_1 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$cell_5q == 'del5q',])
	cell_2 <-  rownames(metadata[metadata$Cluster_names == ct & metadata$cell_5q == 'normal',])
	if(length(cell_1) == 0 || length(cell_2) == 0){
		next
	}
	post_data_pr$Status <-  ifelse(colnames(post_data_pr) %in% cell_1, 'del5qPR', 
							ifelse(colnames(post_data_pr) %in% cell_2, 'normalPR', 'None'))
	Idents(post_data_pr)  <- 'Status'
	results <- FindMarkers(post_data_pr, 
								ident.1 = 'del5qPR',
								ident.2 = 'normalPR',
								assay='RNA') 
								# No need to correct by sample as is within patient and otherwise it fails
								# test.use = 'MAST', 
								# latent.vars='Sample')
	results_DE[[ct]] <- results
	tt <- results
	tt$gene_name <- rownames(tt)
	tt$avg_logFC <- tt$avg_log2FC
	tt <- tt[tt$p_val_adj <= 0.05,]
	if(nrow(tt)>10){
		HM_LIST[[ct]] <- get_HM_TopBottom_MAST(results, post_data_pr, ct )
		print(get_volcano(tt, ct, sub= paste0("del(5q)PR Vs Non-del(5q)PR")))
		tt <- tt[order(tt$avg_logFC, decreasing=TRUE),]
		ranks <- setNames(tt$avg_logFC, tt$gene)
		results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
		results_c2 <- results_c2[results_c2$padj < 0.05,]
		results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
		results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
		results_c2$pathway <- factor(results_c2$pathway)
		if(nrow(results_c2)>0){
			print(get_tornado(results_c2, ct))
		}
	}

}
dev.off()


pdf('./Plots/DE_POST/HM_DE_Del5qPRVsNonDel5qPR_MAST.pdf')
for (ct in cell_types){
	 grid::grid.draw(HM_LIST[[ct]]$gtable)
	 grid::grid.newpage()
}
dev.off()

results_DE <- lapply(results_DE, function(x) {
				x <- x[x$p_val_adj <= 0.05,]
				print(dim(x))
				return(x)
})

WriteXLS::WriteXLS(results_DE, ExcelFileName=paste0('./Plots/DE_POST/Results_DE_Del5qPRVsNonDel5qPR_MAST.xlsx'), SheetNames = names(results_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)


# 3.2 Del(5q)PR VS Elderly

table(elder_data$Sample, elder_data$Cluster_names)
get_prop_tables_per_sample(post_data)

post_data_del <- subset(x = post_data, subset = cell_5q == "del5q")
post_data_del <- subset(x = post_data_del, subset = Sample == "FS-0406-post")
expression_post  <- as.data.frame(post_data_del@assays$RNA@counts)
post_data_del$Status <- 'PR'

expression_elder <- as.data.frame(elder_data@assays$RNA@counts)
elder_data$Status <- 'Elder'

genes_2_keep <- intersect(rownames(expression_post), rownames(expression_elder))

all_data <- merge(expression_post[genes_2_keep,], expression_elder[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_post) + ncol(expression_elder)

metadata <- setNames(rbind(post_data_del[[c('Status', 'Sample', 'Cluster_names')]], 
						   elder_data[[c('Status', 'Sample', 'Cluster_names')]]), 
					c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('PR', 'Elder'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)
table(colnames(all_data) %in% rownames(metadata))

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results <- de_results[de_results$p_val_adj < 0.05,]
de_results <- split( de_results , f = de_results$cell_type )
get_HM_TopBottom(de_results, all_data, metadata, 'Non5qPR_VsElder')

saveDE_results(DE, 'Non5qPR_VsElder')


cell_by_sample <- c(get_cells_by_condition_by_cluster(post_data),
					get_cells_by_condition_by_cluster(MDS_data),
					get_cells_by_condition_by_cluster(pre_post_data)
					)

WriteXLS::WriteXLS(cell_by_sample, ExcelFileName=paste0('./Plots/DE_POST/Cells_by_Sample_by_Condition.xlsx'), 
SheetNames = names(cell_by_sample),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)





get_prop_tables_per_sample(post_data)
get_prop_tables_per_sample(MDS_data)
get_prop_tables_per_sample(pre_post_data)


'CD274' %in% rownames(post_data)
# 'PDL-1' %in% rownames(MDS_data)
'CD274' %in% rownames(pre_post_data)

DefaultAssay(post_data) <- 'SCT'
DefaultAssay(pre_post_data) <- 'SCT'
post_data_meta <- FetchData(post_data, vars=c('Sample', 'Cluster_names', 'cell_5q','CD34'))
pre_post_data_meta <- FetchData(pre_post_data, vars=c('Sample', 'Cluster_names', 'cell_5q','CD34'))

all_data <- rbind(post_data_meta, pre_post_data_meta)
all_data <- reshape2::melt(all_data)

pdf('./Plots/Test_cd274.pdf')
ggplot(all_data, aes(y=value, x=cell_5q, fill=cell_5q)) + 
geom_boxplot() + facet_grid(~Cluster_names) + 
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
facet_wrap(~Sample, scales = 'free')
dev.off()