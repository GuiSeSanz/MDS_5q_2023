



library(reticulate)
reticulate::use_python("/usr/bin/python3")
library(Seurat)
library(cowplot)
library(dplyr)
library(future)
library(viridis)
library(ggplot2)
library(plyr)
library(reshape2)
library(gridExtra)
library(ggridges)
library(stringr)
library(hues)
library(pheatmap)
library(Rtsne)


get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}

plasma <- viridis(50, direction = 1, option = "C")

#### NEED THE INPUT OF THE SEURAT DATA WITH THE NAME OF THE PHENOTYPE

# post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
# all_5q_depleted_cells_COPYKAT <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
# post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')
# #  There are less cells as some clusters have disapeared based on the annotation

# # elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
# MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('./Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')

# prop.table(table(post_data$Sample, post_data$cell_5q), margin=1)*100
prop.table(table(pre_post_data$Sample, pre_post_data$cell_5q), margin=1)*100
# table(MDS_data$Sample, MDS_data$cell_5q)


tmp <- FetchData(pre_post_data , vars=c('Cluster_names', 'cell_5q', 'Sample'))
tmp$cell_id <- rownames(tmp)

tmp$Case <- ifelse(tmp$Sample == 'SMD132114579', 'NR_post', 'NR_Pre')

seurat_data <- tmp
seurat_data$cell_id <- gsub('_', '-', seurat_data$cell_id)
rownames(seurat_data) <- seurat_data$cell_id
# # # The numerical assignments are:
# # 0 <- Pre
# # 1 <- Pots

phenotypes <- c('NR_Pre', 'NR_post')
color_map <- c("#5f0f40","#0f4c5c")
names(color_map) <- phenotypes

# binaryze the phenotypes with base 0
assignment <- as.character(seq(0,length(phenotypes)-1))
names(assignment) <- phenotypes

# Load the weigths and AUCs from SimiC
weights_file <- './Data/SimiC/Results/POST_Run_Del_PREPOST_1000_L10.01_L20.01_Ws_filtered_BIC.pickle'
SimiC_weights <- py_load_object(filename =weights_file)
pd <- import("pandas")
AUCs_file <-'./Data/SimiC/Results/POST_Run_Del_PREPOST_1000_L10.01_L20.01_AUCs.pickle'
aucs <- py_load_object(filename =AUCs_file)
AUCs <- list()
for(i in seq_along(phenotypes)){
  AUCs[[i]] <- aucs[[i]]
}


# filter low R2 targets
unselected_targets <- list()
for (phenotype in phenotypes){
    unselected_targets[[phenotype]] <- SimiC_weights$query_targets[which(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]] < 0.7)]
}

pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_Filtered_R2_targets_hist_PostPREPOST.pdf'))
for (phenotype in phenotypes){
    selectedTargets <- which(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]] > 0.7)
    hist(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]], col='grey', breaks=100, 
    xlab = 'Adjusted R2',
    main = paste0('Phenotype: ',phenotype ,'\n Targets selected: ', length(selectedTargets), ', mean R2: ', mean(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]][selectedTargets])))
}
dev.off()


SimiC_weights_df <- data.frame(driver=NULL, target = NULL, value = NULL, .id = NULL, stringsAsFactors=F)
for(phenotype in phenotypes){
    tmp <-SimiC_weights$weight_dic[[assignment[[phenotype]]]][-nrow(SimiC_weights$weight_dic[[assignment[[phenotype]]]]),]
    rownames(tmp)<-SimiC_weights$TF_ids
    colnames(tmp)<-SimiC_weights$query_targets

    tmp <- setNames(melt(tmp), c('driver', 'target', 'value'))
    tmp$.id <- phenotype
    SimiC_weights_df <- rbind(SimiC_weights_df, tmp)
}




pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_TF_weigths_PostPREPOST.pdf'), onefile = TRUE, width=20)
TF_2_remove <- c()
for(drv in unique(SimiC_weights_df$driver)){
  message(drv)
  tmp_plotter <- SimiC_weights_df[SimiC_weights_df$driver == drv,]
  # remove the unselected targets for each phenotype
  tmp_plotter_filtered <- data.frame(driver=NULL,  target=NULL,  value=NULL,  .id=NULL , stringsAsFactors = F)
  for (phenotype in phenotypes){
      tmp_plotter_filtered <- rbind(tmp_plotter_filtered, tmp_plotter[tmp_plotter$.id == phenotype & !tmp_plotter$target %in% c(unselected_targets[[phenotype]]) ,])
  }
  #
  tmp_plotter_filtered <- tmp_plotter_filtered[order(abs(tmp_plotter_filtered$value), decreasing=T),]
  bests <- unique(tmp_plotter_filtered$target)[1:100]
  tmp_plotter_filtered <- tmp_plotter_filtered[tmp_plotter_filtered$target %in% bests,]
  tmp_plotter_filtered$target <- factor(tmp_plotter_filtered$target, levels = unique(tmp_plotter_filtered$target))

  if(sum(abs(tmp_plotter_filtered$value)) == 0){
    TF_2_remove <- c(TF_2_remove, drv)
    drv <- paste0(drv, '(deleted)')
  }
  p <- ggplot(tmp_plotter_filtered, aes(x=target, y=value, fill=.id)) + 
    geom_bar(stat='identity', position='dodge', color='black') + 
    scale_fill_manual(values=color_map) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(drv) 
  print(p)
}
dev.off()


pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_Targets_weigths_PostPREPOST.pdf'), onefile = TRUE, width=20)
plot_counter <- 1
plot_list <- list()
for(tgt in sort(as.character(unique(SimiC_weights_df$target)))){
  # get the target genes for each TF
  # if the target is discarded on all the phenotypes
  message(tgt)
  if(tgt %in% Reduce(intersect, unselected_targets)){
    next
  # if the target is discarded on at least one of the phenotypes
  }else if(tgt %in% unique(unlist(unselected_targets))){
    tmp_plotter <- data.frame(driver=NULL, target =NULL, value=NULL, .id=NULL)
    for (phenotype in phenotypes){
      if ( !tgt %in% unselected_targets[[phenotype]]){
         tmp_plotter <- rbind(tmp_plotter, SimiC_weights_df[SimiC_weights_df$target == tgt & SimiC_weights_df$.id == phenotype,])
      }
    }  
  }else{
    tmp_plotter <- SimiC_weights_df[SimiC_weights_df$target == tgt,]
  }
  # tmp_plotter$value <- scale(tmp_plotter$value, center=FALSE)
  tmp_plotter <- tmp_plotter[order(abs(tmp_plotter$value), decreasing=T),]
  assign(paste0('p', plot_counter), ggplot(tmp_plotter, aes(x=reorder(driver, -abs(value)), y=value, fill=.id,  palette = "jco")) + 
    geom_bar(stat='identity', position='dodge', color='black') + 
    scale_fill_manual(values=color_map) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(tgt) )
  if(plot_counter == 4){
    grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
    plot_list <- list()
    plot_counter <- 1
  }else{
    plot_counter <- plot_counter +1
  }
}
dev.off()


# AUCs <- list()
# for(i in seq_along(phenotypes)){
# 	AUCs[[i]] <- read.table(paste0('./Data/SimiC/Results//PostPREPOST_1000_L10.01_L20.01_AUCs_', assignment[[i]], '_BIS.csv'), header=T, sep='\t')
# }

# set the phenotype and cluster per TF
# seurat_data$Normal_5q <- ifelse(seurat_data$cell_5q == 'normal', 'non5q', 'real_5q')
# Idents(seurat_data) <- seurat_data$Normal_5q

AUCs_by_state<-list()
state_specific_AUCs<-NULL

for(phenotype in phenotypes){
	print(phenotype)
	index <- as.integer(assignment[phenotype])+1
	# cell_names_phenotype <- colnames(subset(seurat_data, idents = phenotype))
  cell_names_phenotype <- rownames(seurat_data[seurat_data$Case ==phenotype,])
	AUCs_cell_i <- AUCs[[index]]
	# rownames(AUCs_cell_i) <- AUCs_cell_i$X
	# AUCs_cell_i <- AUCs_cell_i[ , -1]
	AUCs_cell_i$TimePoint <- phenotype
  # (stringr::str_extract(rownames(AUCs_cell_i), '(?<=_)[A-Z0-9-_]+$') %in% cell_names_phenotype)
	AUCs_by_state[[phenotype]] <- AUCs_cell_i[rownames(AUCs_cell_i) %in% cell_names_phenotype, ]
	state_specific_AUCs <- rbind(state_specific_AUCs, AUCs_cell_i[rownames(AUCs_cell_i) %in% cell_names_phenotype, ])
}

# df <- do.call('rbind', AUCs_by_state)

for (phenotype in phenotypes){
  if (phenotype == phenotypes[1]){
    df <- AUCs_by_state[[phenotype]]
  }else{
    df <- rbind(df, AUCs_by_state[[phenotype]])
  }
}
df$cell_id <- rownames(df)
# df$cell_id <- stringr::str_extract(rownames(df), '(?<=\\.)[A-Z0-9-_]+$')
# rownames(df) <- df$cell_id

clusters_df <- setNames( seurat_data[, 'Cluster_names', drop=FALSE] , c('cluster_id'))
clusters_df$cell_id <- rownames(clusters_df)
df_w_cluster <- merge(df,clusters_df[, c('cell_id', 'cluster_id')],by="cell_id")

df_auc <- melt(df_w_cluster , id.vars = c('TimePoint','cell_id', 'cluster_id'), variable.name = 'driver', stringsAsFactors =F)



# check the population per phenotype and cluster
tmp <- with(df_w_cluster, table(cluster_id, TimePoint))
tmp <- reshape2::melt(tmp)
 
pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_Population_phenotype_PostPREPOST.pdf'))
  ggplot(tmp, aes(x=cluster_id, y=value, fill=TimePoint)) + 
  geom_bar(stat='identity' , position='dodge',  color='black') + 
  geom_text(aes(label=value), position=position_dodge(width=0.9), angle=90, vjust=0.4, hjust=-0.1) +
  scale_fill_manual(values=color_map) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
dev.off()


##### keep the clusters with >10 representatnts in both populations

clusters_2_keep <- as.data.frame.matrix(with(df_auc, table(cluster_id, TimePoint)))
clusters_2_keep$cluster_id <- rownames(clusters_2_keep)
clusters_2_keep <- clusters_2_keep[rowSums(clusters_2_keep > 20) == ncol(clusters_2_keep), 'cluster_id']


# Keep all cluster or not
df_auc_common <- df_auc[df_auc$cluster_id %in% clusters_2_keep, ]
# df_auc_common <- df_auc
# clusters_2_keep <- paste(unique(df_auc$cluster_id))
# calculate the regulation dissimilarity by the total variation score
MinMax_clust<-NULL
for(cluster in clusters_2_keep){
	print(cluster)
  tmp <- df_auc_common[df_auc_common$cluster_id==cluster,]
  MinMax_val <- NULL
  for (tf in unique(df_auc_common$driver)){
    n_breaks = 100
    Wauc_dist <- list()
    for (phenotype in phenotypes){
      Wauc_dist[[phenotype]] <- hist(tmp[tmp$TimePoint ==phenotype & tmp$driver ==tf, 'value'], breaks = c(seq(0,1, 1/n_breaks)), plot=FALSE)$density
    }
    mat <-  do.call('rbind', Wauc_dist)
    mat <- mat[complete.cases(mat),]
    minmax_diff <- apply(na.omit(mat), 2, max) - apply(na.omit(mat), 2, min)
    variant <- sum(abs(minmax_diff)) / n_breaks
    variant <- variant/ sum(rowSums(mat)!=0)
    MinMax_val <- append(MinMax_val,variant)
  }
  MinMax_val <- setNames(as.data.frame(MinMax_val), c(cluster))
  rownames(MinMax_val) <-  unique(df_auc_common$driver)
  if (is.null(MinMax_clust)){
    MinMax_clust <- MinMax_val
  }else{
    MinMax_clust<-cbind(MinMax_clust,MinMax_val)
  }
}



# plot the densities of the AUC and the score per TF

for (clust in clusters_2_keep){
  print(clust)
  plotter2 <- df_auc[df_auc$cluster_id == clust,]
  clust_name <- sub(' ', '_', clust)
  pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_Auc_Cluster_',clust,'_PostPREPOST.pdf'), width = 15, onefile = TRUE)
  plot_counter <- 1
  for (tf in unique(plotter2$driver)){
    assign( paste0('p', plot_counter), 
    ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=TimePoint)) + 
    geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
    scale_fill_manual(values=color_map) +
    theme(legend.position = 'top')+ geom_rug() + 
    ggtitle(paste0(tf))) #, '   ', MinMax_clust[rownames(MinMax_clust) == tf, clust])) )
    if(plot_counter == 2){
      grid.arrange(p1, p2, ncol=2)
      plot_counter <- 1
    }else{
      plot_counter <- plot_counter +1
    }
  }
  grid.arrange(p1, p2 , ncol=2)
  dev.off()
}

saveRDS(df_auc, '/home/tereshkova/data/gserranos/MDS/Data/SimiC_df_auc_PostPREPOST.rds')

pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_Auc_all_cells_PostPREPOST.pdf'), width = 15, onefile = TRUE)
plot_counter <- 1
plotter2 <- df_auc
for (tf in unique(plotter2$driver)){
	message(tf)
	assign( paste0('p', plot_counter), 
			ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=TimePoint)) + 
			geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
			scale_fill_manual(values=color_map) +
			theme(legend.position = 'top') + 
			ggtitle(paste0(tf))
		) #, '   ', MinMax_clust[rownames(MinMax_clust) == tf, clust])) )
	if(plot_counter == 2){
		grid.arrange(p1, p2, ncol=2)
		plot_counter <- 1
	}else{
		plot_counter <- plot_counter +1
	}
}
grid.arrange(p1, p2 , ncol=2)
dev.off()


# pdf(paste0('./Plots/SimiC/Simic_Auc_all_cells_perPatient_PostPREPOST.pdf'), width = 15, onefile = TRUE)
# plot_counter <- 1
# for (tf in unique(plotter2$driver)){
# 	assign( paste0('p', plot_counter), 
# 	ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=TimePoint)) + 
# 	geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
# 	scale_fill_iwanthue() +
# 	theme(legend.position = 'top')+ geom_rug() + 
# 	ggtitle(paste0(tf))) #, '   ', MinMax_clust[rownames(MinMax_clust) == tf, clust])) )
# 	if(plot_counter == 2){
# 		grid.arrange(p1, p2, ncol=2)
# 		plot_counter <- 1
# 	}else{
# 		plot_counter <- plot_counter +1
# 	}
# }
# grid.arrange(p1, p2 , ncol=2)
# dev.off()



clust_order_asc<-names(sort(apply(MinMax_clust, 2, mean)))
MinMax_clust<-MinMax_clust[,clust_order_asc]
MinMax_df <- as.data.frame(MinMax_clust)
MinMax_df$driver <- rownames(MinMax_df)
MinMax_df <- melt(MinMax_df,variable.name = "cluster_id")
MinMax_clust <- MinMax_clust[!rownames(MinMax_clust) %in% TF_2_remove,]

saveRDS(MinMax_clust, '/home/tereshkova/data/gserranos/MDS/Data/SimiC_MinMax_clust_PostPREPOST.rds')

# Plot the heatmap of the regulatory dissimilarity score
p <- pheatmap::pheatmap(MinMax_clust,color=plasma, fontsize=5, angle_col =45, cellwidth=40, silent=TRUE)
plot_dims <- get_plot_dims(p)
pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_HeatMaps_RegDissScore_PostPREPOST.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()




# all_data <- as.data.frame( seurat_data@assays$RNA@scale.data)
# all_data <- all_data[, colnames(all_data) %in% unique(df_auc[df_auc$TimePoint %in% phenotypes, 'cell_id'])]
# tsne_out <- Rtsne(as.matrix(t(all_data)))
# plotter <-  as.data.frame(tsne_out$Y)
# rownames(plotter) <- colnames(all_data)
# plotter$cell_id <- colnames(all_data)


# genes_2_plot <- c('MYC', 'RUNX3', 'EOMES')
# pdf(paste0('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Simic_tSNE_selection_PostPREPOST.pdf'), height =8, width=18)
# for (gene in genes_2_plot){
#   plotter_tmp <- plotter[rownames(plotter) %in% df_auc[df_auc$TimePoint %in% phenotypes, 'cell_id'],]
#   plotter_tmp <- merge(plotter_tmp, df_auc[df_auc$driver == gene, c('cell_id', 'value', 'TimePoint', 'cluster_id')], by='cell_id')
#   plotter_tmp$cluster_id <- factor(plotter_tmp$cluster_id, levels=c('CD8 Memory', 'CD8 Effector', 'CD8 Exhausted'))

#   g <- ggplot() +
#   geom_point(plotter_tmp[,-6], mapping=aes(x=V1, y=V2, shape = TimePoint), color = "grey", alpha = 0.4) + 
#   geom_point(plotter_tmp, mapping=aes(x=V1, y=V2, color = value, shape = TimePoint)) + 
#   theme_classic() + geom_point(size = 1.5) + 
#   scale_color_viridis(option='inferno') + 
#   facet_wrap(~cluster_id, ncol=3, nrow=1) + labs(x= '', y = gene )
#   print(g)
# }

# dev.off()





saveRDS(SimiC_weights_df,  '/home/tereshkova/data/gserranos/MDS/Data/SimiC_Weigths_PostPREPOST.rds')
saveRDS(unselected_targets, '/home/tereshkova/data/gserranos/MDS/Data/SimiC_unselected_targets_PostPREPOST.rds')

pdf('/home/tereshkova/data/gserranos/MDS/Plots/SimiC/Post/PREPOST/Enrich_Regulons_PostPREPOST.pdf')
theme_mod <-  theme(axis.text.y=element_text(angle=45, vjust=0.5, size=8), legend.position='none')
for(drv in sort(as.character(unique(SimiC_weights_df$driver)))){
  print(drv)
  targets <- SimiC_weights_df[SimiC_weights_df$driver == drv & SimiC_weights_df$value != 0,'target']
  names = clusterProfiler::bitr(sub('\\.','-', targets ), fromType="SYMBOL", toType= "ENTREZID", OrgDb="org.Hs.eg.db")
  results_BP = clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "BP",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE)
  results_MF = clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "MF",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE)
  results_CC = clusterProfiler::enrichGO(names[, 'ENTREZID'], OrgDb= 'org.Hs.eg.db', ont= "CC",  pAdjustMethod = "BH", pvalueCutoff  = 0.05,  qvalueCutoff  = 0.05,readable= TRUE)
  results_KG = clusterProfiler::enrichKEGG(names[, 'ENTREZID'], organism = "hsa",pAdjustMethod = "BH",pvalueCutoff  = 0.05)
  title_theme <- ggdraw() +draw_label(drv)

  plot_list <- list()
  if( results_BP@readable == TRUE & nrow(head(as.data.frame(results_BP)>0 ))){ plot_list[['BP']] <- clusterProfiler::dotplot(results_BP) + theme_mod}else{NULL}
  if( results_MF@readable == TRUE & nrow(head(as.data.frame(results_MF)>0 ))){ plot_list[['MF']] <- clusterProfiler::dotplot(results_MF) + theme_mod}else{NULL}
  if( results_CC@readable == TRUE & nrow(head(as.data.frame(results_CC)>0 ))){ plot_list[['CC']] <- clusterProfiler::dotplot(results_CC) + theme_mod}else{NULL}
  # if( results_KG@readable == TRUE & nrow(head(as.data.frame(results_KG)>0 ))){ plot_list[['KG']] <- clusterProfiler::dotplot(results_KG) + theme_mod}else{NULL}
  
  print(cowplot::plot_grid(title_theme, 
        cowplot::plot_grid(plotlist=plot_list,
                           ncol=2, labels=c('BP', 'MF', 'CC')), 
                           nrow=2, rel_heights=c(0.1,0.9)))
}
dev.off()