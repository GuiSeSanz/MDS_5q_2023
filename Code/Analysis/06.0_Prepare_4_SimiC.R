
library(reticulate)
reticulate::use_python("/usr/bin/python3")
py_discover_config("magic")
library(Rmagic)

library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
 



#################################################################
# SimiC 5q Vs elder
#################################################################

data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_integrated_sct_subcluster.rds')
data_5q <- subset(data_5q, idents = c( "17_0", "23", "25", "28", "29"), invert = TRUE)

data_elder <-  readRDS(paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_subcluster.rds'))
data_elder <- subset(data_elder, idents = c('5_0', '6_0', '6_1', '12_0', '13_0', '13_1', '14_1', '20_0', '22'), invert = TRUE)

all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') 
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')
real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
real_5q_cells <- all_5q_selected_cells[all_5q_selected_cells$Cell_id %in% real_5q_cells,]



if(!file.exists('/home/sevastopol/data/gserranos/MDS/Data/imputed_MAGIC_5q_and_elder.rds')){
	# We just need an object with the 5q samples and the elder count's, no need the integrated array.
	data_5q <- subset(data_5q, cells=real_5q_cells$Cell_id)
	combined.sct <- merge(data_5q, y= data_elder, add.cell.ids=c('real_5q', 'elder'), project='Combined_5q_elder')
	ncol(combined.sct) == ncol(data_elder) + ncol(data_5q)
	combined.sct_raw <- as.data.frame(combined.sct@assays$RNA@counts)
	unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
	combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
	### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
	combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

	# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
	combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes') 

	saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_5q_and_elder.rds')

}else{
	combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_5q_and_elder.rds')
}

MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df), sample = stringr::str_extract(rownames(data_MAGIC_df), "(^[a-z_5]+)(?!S|G)"))
cell_info$assignment <- ifelse(cell_info$sample == 'real_5q', 0, 1)
# The numerical assignments are:
# 0 <- real_5q
# 1 <- elder
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/5qVsElder.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/5qVsElder_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/5qVsElder_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/5qVsElder_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/5qVsElder_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)







#################################################################
# SimiC 5q Vs Non5q
#################################################################



sc_data_MDS5q_subseted <- readRDS('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')

if(!file.exists('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_5qVsnon5q.rds')){
	combined.sct_raw <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
	unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
	combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
	### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
	combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

	# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
	combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes') 

	saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_5qVsnon5q.rds')
}else{
 combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_5qVsnon5q.rds')
}
MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df), sample = stringr::str_extract(rownames(data_MAGIC_df), "(SMD[0-9]+)(?=_)"))

info_5q <- setNames(as.data.frame(sc_data_MDS5q_subseted$cell_5q), 'assignment')
cell_info <- merge(cell_info, info_5q, by.x='cell_id', by.y=0)
cell_info$assignment <- ifelse(cell_info$assignment == 'del5q', 0, 1)

# The numerical assignments are:
# 0 <- real_5q
# 1 <- normal
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/5qVsnon5q.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/5qVsnon5q_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/5qVsnon5q_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/5qVsnon5q_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t', row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/5qVsnon5q_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)






#################################################################
# Pre Vs Post
#################################################################



sc_data_PrePost <- readRDS('./Data/pre_post_5q_Annotated_final.rds')

if(!file.exists('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_PreVsPost.rds')){
	combined.sct_raw <- as.data.frame(sc_data_PrePost@assays$RNA@counts)
	unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
	combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
	### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
	combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

	# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
	combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes') 
	saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_PreVsPost.rds')
}else{
 combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_PreVsPost.rds')
}

MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)

cell_info <- data.frame(cell_id =rownames(data_MAGIC_df), sample = stringr::str_extract(rownames(data_MAGIC_df), "(SMD[0-9]+)(?=_)"))

# Keep only the cells with 5q
assignment <- readRDS('./Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
cell_info <- cell_info[cell_info$cell_id %in% assignment, ]
cell_info$assignment <- ifelse(cell_info$sample == 'SMD211420', 0, 1)

# The numerical assignments are:
# 0 <- Pre
# 1 <- Post
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/simicprepost.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/PreVsPost_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/PreVsPost_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/PreVsPost_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t', row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/PreVsPost_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)





#################################################################
# Non5qVsElder
#################################################################

data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_elder <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')


if(!file.exists('/home/sevastopol/data/gserranos/MDS/Data/imputed_MAGIC_5q_and_elder.rds')){
	all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') 
	cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')
	# remove the 5qcells and the maybes
	real_5q_cells <- union(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
	data_5q <- as.data.frame(data_5q@assays$RNA@counts)
	data_5q <- data_5q[, !colnames(data_5q) %in% real_5q_cells]
	data_5q_gene_names <-  sapply(rownames(data_5q) , FUN=function(gene_name) gsub('\\.', '_', gene_name))
	data_elder <-as.data.frame(data_elder@assays$RNA@counts)
	data_elder_gene_names <-  sapply(rownames(data_elder) , FUN=function(gene_name) gsub('\\.', '_', gene_name))
	rownames(data_5q) <- data_5q_gene_names
	rownames(data_elder) <- data_elder_gene_names
	combined.sct_raw  <- merge(data_5q, data_elder, by=0)
	rownames(combined.sct_raw) <- combined.sct_raw$Row.names
	combined.sct_raw <- combined.sct_raw[,-1]
	unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
	combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
	### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
	combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

	# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
	combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes') 

	saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_Non5qVsElder.rds')

}else{
	combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_Non5qVsElder.rds')
}

MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df), sample = stringr::str_extract(rownames(data_MAGIC_df), "(^[A-Z0-9]+)"))
cell_info$assignment <- ifelse(grepl('SMD', cell_info$sample), 0, 1)
# The numerical assignments are:
# 0 <- Non_5q
# 1 <- elder
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/Non5qVsElder.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/Non5qVsElder_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/Non5qVsElder_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/Non5qVsElder_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/Non5qVsElder_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)







#################################################################
# 5qVsNon5qVsElder
#################################################################


data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
data_elder <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

data_5q_counts    <- as.data.frame(data_5q@assays$RNA@counts)
data_elder_counts <- as.data.frame(data_elder@assays$RNA@counts)
data_5q_gene_names    <-  sapply(rownames(data_5q_counts) , FUN=function(gene_name) gsub('\\.', '_', gene_name))
data_elder_gene_names <-  sapply(rownames(data_elder_counts) , FUN=function(gene_name) gsub('\\.', '_', gene_name))
rownames(data_5q_counts) <- data_5q_gene_names
rownames(data_elder_counts) <- data_elder_gene_names
combined.sct_raw  <- merge(data_5q_counts, data_elder_counts, by=0)
rownames(combined.sct_raw) <- combined.sct_raw$Row.names
combined.sct_raw <- combined.sct_raw[,-1]
unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes')

saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_All_5qVsElder.rds')
combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_All_5qVsElder.rds')

MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df), sample = stringr::str_extract(rownames(data_MAGIC_df), "(^[A-Z0-9]+)"))
# cell_info$assignment <- ifelse(grepl('SMD', cell_info$sample), 0, 1)
cell_5q <- as.data.frame(data_5q$cell_5q)
cell_5q <- rownames(cell_5q[cell_5q[[1]] == 'del5q',,drop=F])
cell_info$assignment <- ifelse(cell_info$cell_id %in% cell_5q, 0, 
						ifelse(grepl('SMD', cell_info$sample), 1, 2))
# The numerical assignments are:
# 0 <- 5q
# 1 <- non5q
# 2 <- elder
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/ThreeWayComparison.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/ThreeWayComparison_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/ThreeWayComparison_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/ThreeWayComparison_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/ThreeWayComparison_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)



