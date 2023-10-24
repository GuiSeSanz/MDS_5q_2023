
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
#################################################################

post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
all_5q_depleted_cells_COPYKAT <- readRDS('./Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')
#  There are less cells as some clusters have disapeared based on the annotation

elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('./Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')


if(!file.exists('/home/sevastopol/data/gserranos/MDS/Data/imputed_MAGIC_ALL_COMBINED.rds')){
	# We just need an object with the 5q samples and the elder count's, no need the integrated array.
	combined.sct <- merge(post_data, y=c( elder_data, pre_post_data, MDS_data) , 
						add.cell.ids=c('POST', 'ELDER', 'PREPOST', 'MDS'), project='ALL_COMBINED')
	ncol(combined.sct) == ncol(elder_data) + ncol(pre_post_data) + ncol(MDS_data) + ncol(post_data)
	combined.sct_raw <- as.data.frame(combined.sct@assays$RNA@counts)
	unexpresed_genes <- names(which(rowSums(abs(combined.sct_raw ))<1e-6))
	combined.sct_raw  <- combined.sct_raw [ !rownames(combined.sct_raw ) %in% unexpresed_genes, ]
	### RUN MAGIC #### #MAGIC works on a cells x genes matrix, seurat gives a genes x cells matrix
	combined.sct_raw_imputedMAGIC <-Rmagic::library.size.normalize(t(combined.sct_raw ))

	# combined.sct_raw_imputedMAGIC <- sqrt(combined.sct_raw_imputedMAGIC)
	combined.sct_raw_imputedMAGIC <- magic(combined.sct_raw_imputedMAGIC,genes='all_genes') 

	saveRDS(combined.sct_raw_imputedMAGIC, '/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_ALL_COMBINED.rds')

}else{
	combined.sct_raw_imputedMAGIC <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/imputed_MAGIC_ALL_COMBINED.rds')
}

MAX_NUM_TFs=100
MAX_NUM_TARGETS=1000
data_MAGIC_df <- as.data.frame(combined.sct_raw_imputedMAGIC)
rownames(data_MAGIC_df) <- gsub('_', '-', rownames(data_MAGIC_df))
rownames(data_MAGIC_df) <- stringr::str_remove(rownames(data_MAGIC_df), '^(?:POST|PREPOST|MDS|ELDER)-')



################################################################################################
# RUN 1 :: Normal cells   ELDER -> MDS -> PR -> CR
################################################################################################

normal_cells <- c( colnames(elder_data), 
					colnames(MDS_data)[MDS_data$cell_5q == 'normal'], 
					colnames(post_data)[post_data$cell_5q == 'normal'])
normal_cells <- gsub('_', '-', normal_cells)

data_MAGIC_df_Normal <- data_MAGIC_df[rownames(data_MAGIC_df) %in% normal_cells, ]
nrow(data_MAGIC_df_Normal) == length(normal_cells)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df_Normal), sample = stringr::str_remove(rownames(data_MAGIC_df_Normal), '[-A-Z]+-1$'))
cell_info$assignment <- ifelse(grepl('^GSM', cell_info$sample), 0, 
						ifelse(grepl('^SMD', cell_info$sample), 1,
						ifelse(grepl('^FS-0406', cell_info$sample), 2, 3)))
table(cell_info$assignment, cell_info$sample)
# # The numerical assignments are:
# # 0 <- Elder
# # 1 <- MDS
# # 2 <- PR
# # 3 <- CR
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/POST_Run_Normal.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df_Normal),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df_Normal[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df_Normal),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df_Normal[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df_Normal[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/POST_Run_Normal_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/POST_Run_Normal_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/POST_Run_Normal_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/POST_Run_Normal_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)




################################################################################################
# RUN 2 :: del(5q) cells   MDS -> PR -> PRE -> Post
################################################################################################

del5q <- c( colnames(MDS_data)[MDS_data$cell_5q == 'del5q'], 
					colnames(post_data)[post_data$cell_5q == 'del5q' & post_data$Sample == 'FS-0406-post'], 
					colnames(pre_post_data)[pre_post_data$cell_5q == 'del5q'])
del5q <- gsub('_', '-', del5q)

data_MAGIC_df_Del<- data_MAGIC_df[rownames(data_MAGIC_df) %in% del5q, ]
nrow(data_MAGIC_df_Del) == length(del5q)
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df_Del), sample = stringr::str_remove(rownames(data_MAGIC_df_Del), '[-A-Z]+-1$'))
cell_info$assignment <- ifelse(grepl('^SMD3', cell_info$sample), 0,
						ifelse(grepl('^FS-0406', cell_info$sample), 1,
						ifelse(grepl('^SMD211420', cell_info$sample), 2, 
						ifelse(grepl('^SMD132114579', cell_info$sample), 3, 4))))
#        FS-0406-post SMD132114579 SMD211420 SMD34459 SMD35109 SMD35303 SMD37209
#   0            0            0         0     1607    10604    10667     4239
#   1         1939            0         0        0        0        0        0
#   2            0            0     12409        0        0        0        0
#   3            0         8522         0        0        0        0        0
# # # The numerical assignments are:
# # 0 <- MDS
# # 1 <- PR
# # 2 <- NR_Pre
# # 3 <- NR_Post
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/POST_Run_Del1.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df_Del),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df_Del[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df_Del),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df_Del[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df_Del[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/POST_Run_Del1_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/POST_Run_Del1_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/POST_Run_Del1_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/POST_Run_Del1_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)

################################################################################################
# Run 3 ::  MDS(x4) -> NR(pre) -> PR ->  -> NR(post)
################################################################################################

del5q <- c( colnames(MDS_data)[MDS_data$cell_5q == 'del5q'], 
					colnames(post_data)[post_data$cell_5q == 'del5q' & post_data$Sample == 'FS-0406-post'], 
					colnames(pre_post_data)[pre_post_data$cell_5q == 'del5q'])
del5q <- gsub('_', '-', del5q)

data_MAGIC_df_Del<- data_MAGIC_df[rownames(data_MAGIC_df) %in% del5q, ]
# Only changes the assignation order, so the data should stay the same BUT I will create new files to avoid confusion
cell_info <- data.frame(cell_id =rownames(data_MAGIC_df_Del), sample = stringr::str_remove(rownames(data_MAGIC_df_Del), '[-A-Z]+-1$'))
cell_info$assignment <- ifelse(grepl('^SMD3', cell_info$sample), 0,
						ifelse(grepl('^SMD211420', cell_info$sample), 1, 
						ifelse(grepl('^FS-0406', cell_info$sample), 2,
						ifelse(grepl('^SMD132114579', cell_info$sample), 3, 4))))
#        FS-0406-post SMD132114579 SMD211420 SMD34459 SMD35109 SMD35303 SMD37209
#   0            0            0         0     1607    10604    10667     4239
#   1            0            0     12409        0        0        0        0
#   2         1939            0         0        0        0        0        0
#   3            0         8522         0        0        0        0        0
# # # The numerical assignments are:
# # 0 <- MDS
# # 1 <- NR_pre
# # 2 <- PR
# # 3 <- NR_post
assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/POST_Run_Del2.clustAssign.txt', col.names=FALSE, row.names=FALSE)



reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/POST_Run_Del2_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/POST_Run_Del2_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/POST_Run_Del2_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/POST_Run_Del2_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)



################################################################################################
# Run 4 ::  MDS(x4) -> NR(post) -> PR 
################################################################################################
del5q <- c( colnames(MDS_data)[MDS_data$cell_5q == 'del5q'], 
			colnames(post_data)[post_data$cell_5q == 'del5q' & post_data$Sample == 'FS-0406-post'], 
			colnames(pre_post_data)[pre_post_data$cell_5q == 'del5q' & pre_post_data$Sample == 'SMD132114579'])

del5q <- gsub('_', '-', del5q)

data_MAGIC_df_Del<- data_MAGIC_df[rownames(data_MAGIC_df) %in% del5q, ]

cell_info <- data.frame(cell_id =rownames(data_MAGIC_df_Del), sample = stringr::str_remove(rownames(data_MAGIC_df_Del), '[-A-Z]+-1$'))
cell_info$assignment <- ifelse(grepl('^SMD3', cell_info$sample), 0,
						ifelse(grepl('^SMD132114579', cell_info$sample), 1, 
						ifelse(grepl('^FS-0406', cell_info$sample), 2, 9))) # Add 9 as security check


#     FS-0406-post SMD132114579 SMD34459 SMD35109 SMD35303 SMD37209
#   0            0            0     1607    10604    10667     4239
#   1            0         8522        0        0        0        0
#   2         1939            0        0        0        0        0
# # # The numerical assignments are:
# # 0 <- MDS
# # 1 <- NR_post
# # 2 <- PR

assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/POST_Run_Del3.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df_Del),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df_Del[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]])

target_genes <- setdiff(colnames(data_MAGIC_df_Del),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df_Del[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df_Del[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/POST_Run_Del3_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/POST_Run_Del3_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/POST_Run_Del3_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/POST_Run_Del3_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)






################################################################################################
# Run 5 ::  del(5q) NR_Pre -> del(5q) NR_Post
################################################################################################


del5q <- c(colnames(pre_post_data)[pre_post_data$cell_5q == 'del5q'])
del5q <- gsub('_', '-', del5q)

data_MAGIC_df_Del <- data_MAGIC_df[rownames(data_MAGIC_df) %in% del5q, ]
nrow(data_MAGIC_df_Del) == length(del5q)

cell_info <- data.frame(cell_id =rownames(data_MAGIC_df_Del), sample = stringr::str_remove(rownames(data_MAGIC_df_Del), '[-A-Z]+-1$'))
cell_info$assignment <- ifelse(grepl('^SMD211420', cell_info$sample), 0,
						ifelse(grepl('^SMD132114579', cell_info$sample), 1, 9)) # Add 9 as security check

# table(cell_info$assignment, cell_info$sample)
#                    0     1
#   SMD132114579     0  8522
#   SMD211420    12409     0
# # # The numerical assignments are:
# # 0 <- Pre
# # 1 <- Pots

assignment <- cell_info$assignment
write.table(assignment,'./Data/SimiC/POST_Run_Del_PREPOST.clustAssign.txt', col.names=FALSE, row.names=FALSE)

TFs_list <- py_load_object('/home/tereshkova/data/gserranos/MDS/Data/SimiC/human_TF.pickle')
TFs <- intersect(colnames(data_MAGIC_df_Del),TFs_list)

MAD_TFs <- order(apply(data_MAGIC_df_Del[,TFs],2,mad), decreasing = TRUE)
top_MAD_tfs <- c('TP53', na.omit(TFs[MAD_TFs[1:MAX_NUM_TFs]]))

target_genes <- setdiff(colnames(data_MAGIC_df_Del),TFs_list)

MAD_targets <- order(apply(data_MAGIC_df_Del[,target_genes],2,mad), decreasing = TRUE)
top_MAD_targets <- na.omit(target_genes[MAD_targets[1:MAX_NUM_TARGETS]])

input_data <- data_MAGIC_df_Del[,c(top_MAD_tfs,top_MAD_targets)]
TFs <- top_MAD_tfs

reticulate::py_save_object(as.data.frame(input_data), filename = paste0('./Data/SimiC/POST_Run_Del_PREPOST_',MAX_NUM_TARGETS, ".DF.pickle"))
reticulate::py_save_object(TFs, filename = paste0('./Data/SimiC/POST_Run_Del_PREPOST_',MAX_NUM_TARGETS, ".TF.pickle"))

write.table(input_data, file=  paste0('./Data/SimiC/POST_Run_Del_PREPOST_',MAX_NUM_TARGETS, ".DF.csv"), sep='\t',  row.names= TRUE, col.names=TRUE, quote=FALSE)
write.table(TFs, file= paste0('./Data/SimiC/POST_Run_Del_PREPOST_',MAX_NUM_TARGETS, ".TF.csv"), sep='\t', row.names= FALSE, col.names=FALSE, quote=FALSE)

