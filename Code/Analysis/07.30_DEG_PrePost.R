library(Seurat)
library(future)
plan("multicore", workers = 64)


Samples_PrePost <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
PLOT_PATH <- '/home/tereshkova/data/gserranos/MDS/Data/DE_results'

cells_selected_CASPER <- readRDS('./Data/CASPER/SelectedCells5q_PrePost_CASPER.rds')
all_5q_selected_CopyKat <- readRDS('./Data/CopyKat/all_5q_selected_cellsPrePost.rds')
real_5q_cells <- readRDS('./Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')

# maybe_5q <- c(setdiff(all_5q_selected_CopyKat$Cell_id, real_5q_cells), setdiff(cells_selected_CASPER$X1, real_5q_cells))
# table(cells_selected_CASPER$X1 %in% real_5q_cells)
# tmp <- Samples_PrePost[[c('Cluster_names', 'Sample')]]
# tmp$Cell_id <- rownames(tmp)
# tmp$is_5q <- ifelse(tmp$Cell_id %in% real_5q_cells, '5q', 
# 					ifelse(tmp$Cell_id %in% maybe_5q, 'maybe_5q', 'normal'))
Samples_PrePost$cell_id <- colnames(Samples_PrePost)
Samples_PrePost$cell_5q <- ifelse(Samples_PrePost$cell_id %in% real_5q_cells, '5q','normal')

############################################################################################################
# 5qPre Vs 5qPost
############################################################################################################

DefaultAssay(Samples_PrePost) <- 'SCT'
res <- 'integrated_snn_res.0.8_sub'
message(res)
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster in sort(unique(Samples_PrePost[['Cluster_names']][[1]]))){
	message(cluster)
	data <- subset(Samples_PrePost, subset = Cluster_names == cluster)
	data <- subset(data, subset = cell_5q == "5q")
	data <- PrepSCTFindMarkers(data)
	Idents(data) <- 'Sample'
	tmp <- FindMarkers(data,
					ident.1 = 'SMD211420',
					ident.2 = 'SMD132114579')
	tmp <- tmp[order(tmp$p_val_adj),]
	tmp <- tmp[tmp$p_val_adj < 0.05 ,]
	tmp$gene_name <- rownames(tmp)
	tmp$cluster <- cluster
	tmp <- tmp[, c('gene_name', 'cluster',  'p_val', 'avg_log2FC', 'p_val_adj')]
	Results <- rbind(Results, tmp)
	
}

saveRDS(Results, paste0(PLOT_PATH,'/PrePost_DE_5qPreVs5qPost.rds'))
Results_sp <- split(Results, Results$cluster)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'/PrePost_DE_5qPreVs5qPost.xlsx'), 
SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)



############################################################################################################
# Non5qPre Vs Non5qPost
############################################################################################################

DefaultAssay(Samples_PrePost) <- 'SCT'
res <- 'integrated_snn_res.0.8_sub'
message(res)
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster in sort(unique(Samples_PrePost[['Cluster_names']][[1]]))){
	message(cluster)
	data <- subset(Samples_PrePost, subset = Cluster_names == cluster)
	data <- subset(data, subset = cell_5q == "normal")
	data <- PrepSCTFindMarkers(data)
	Idents(data) <- 'Sample'
	tmp <- FindMarkers(data,
					ident.1 = 'SMD211420',
					ident.2 = 'SMD132114579')
	tmp <- tmp[order(tmp$p_val_adj),]
	tmp <- tmp[tmp$p_val_adj < 0.05 ,]
	tmp$gene_name <- rownames(tmp)
	tmp$cluster <- cluster
	tmp <- tmp[, c('gene_name', 'cluster',  'p_val', 'avg_log2FC', 'p_val_adj')]
	Results <- rbind(Results, tmp)
	
}

saveRDS(Results, paste0(PLOT_PATH,'/PrePost_DE_Non5qPreVsNon5qPost.rds'))
Results_sp <- split(Results, Results$cluster)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'/PrePost_DE_Non5qPreVsNon5qPost.xlsx'), 
SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

