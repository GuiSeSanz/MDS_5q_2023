
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


# MDS 5q
annotation_Sofia_08_MDS5q <- list()
annotation_Sofia_08_MDS5q[['0_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['0_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['0_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['0_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['0_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['0_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['0_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['1_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['1_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['1_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['1_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['1_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['1_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['1_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['2_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['2_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['2_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['2_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['2_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['2_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['2_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['3_0'   ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['3_1'   ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['3_2'   ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['3_3'   ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['3_4'   ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['4'   ]] <- 'pro-B'
annotation_Sofia_08_MDS5q[['5_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['5_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['5_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['5_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['5_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['5_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['5_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['6'   ]] <- 'Monocytes'
annotation_Sofia_08_MDS5q[['7_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['7_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['7_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['7_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['7_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['7_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['7_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['8_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['8_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['8_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['8_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['8_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['8_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['8_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['9'   ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['10'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['11'  ]] <- 'Granulocyte'
annotation_Sofia_08_MDS5q[['12_0']] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['12_1']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['12_2']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['13'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['14'  ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['15_0'  ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['15_1'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['15_2'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['15_3'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['15_4'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['15_5'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['15_6'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['16'  ]] <- 'MK_Prog'
# annotation_Sofia_08_MDS5q[['17_0']] <- '??'
annotation_Sofia_08_MDS5q[['17_1_0']] <- 'LMPP'
annotation_Sofia_08_MDS5q[['17_1_1']] <- 'HSC'
annotation_Sofia_08_MDS5q[['17_1_2']] <- 'HSC'
annotation_Sofia_08_MDS5q[['17_1_3']] <- 'GMP'
annotation_Sofia_08_MDS5q[['17_1_4']] <- 'CLP'
annotation_Sofia_08_MDS5q[['17_1_5']] <- 'CLP'
annotation_Sofia_08_MDS5q[['17_1_6']] <- 'GMP'
annotation_Sofia_08_MDS5q[['18_0'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['18_1'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['18_2'  ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['18_3'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['18_4'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['19'  ]] <- 'pro-B'
annotation_Sofia_08_MDS5q[['20'  ]] <- 'Basophil'
annotation_Sofia_08_MDS5q[['21_0'  ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['21_1'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['21_2'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['21_3'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['21_4'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['21_5'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['21_6'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['22'  ]] <- 'T'
# annotation_Sofia_08_MDS5q[['23']] <- 
annotation_Sofia_08_MDS5q[['24'  ]] <- 'Monocytes'
# annotation_Sofia_08_MDS5q[['25']] <- 
annotation_Sofia_08_MDS5q[['26'  ]] <- 'DendriticCell'
annotation_Sofia_08_MDS5q[['27'  ]] <- 'LateErythroid'
# annotation_Sofia_08_MDS5q[['28']] <- 
# annotation_Sofia_08_MDS5q[['29']] <- 



annotation_lists <- list()
annotation_lists[['5qSamples']]     <- annotation_Sofia_08_MDS5q


filter_and_annotate <- function(data, annotation){
    clusters <- setNames(as.data.frame(data$integrated_snn_res.0.8_sub), c('Cluster'))
    clusters <- clusters[paste(clusters$Cluster) %in% names(annotation),, drop=FALSE]
    clusters$Cluster_names <- apply(clusters, 1, function(x){ annotation[as.character(x[['Cluster']])][[1]] })
    Idents(data) <- 'integrated_snn_res.0.8_sub'
    data <- subset(data, idents= c(names(annotation)))
    data$Cluster_names <- clusters$Cluster_names
    return(data)
}


get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  dev.off()
  return(list(height = plot_height, width = plot_width))
}


# genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
# genes_expanded <- as.character(unique(genes_expanded$V5))


# clusters_df <- readRDS('./Data/Annotation.rds') # data.frame having the annotation per cell with columns: TimePoint = phenotypes, cell_id = cell identifier, cluster_id = cluster of the cell
plasma <- viridis(50, direction = 1, option = "C")

#### NEED THE INPUT OF THE SEURAT DATA WITH THE NAME OF THE PHENOTYPE
# seurat_data <- readRDS('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')
sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
sc_data_MDS5q_subseted <- filter_and_annotate(sc_data_MDS5q_subseted, annotation_lists[['5qSamples']])

tmp <- data.frame(Case= sc_data_MDS5q_subseted$cell_5q, Cluster_name =sc_data_MDS5q_subseted$Cluster_names)
tmp <- tmp[tmp$Case %in% c('del5q', 'normal'), ]

seurat_data <- tmp



# The phenotypes studied, be aware that must have the same order to the Simic input 
# 0 <- del5q
# 1 <- normal
phenotypes <- c('del5q', 'normal')

# binaryze the phenotypes with base 0
assignment <- as.character(seq(0,length(phenotypes)-1))
names(assignment) <- phenotypes

# Load the weigths and AUCs from SimiC
weights_file <- './Data/SimiC/Results/5qVsnon5q_1000_L10.01_L20.1_Ws_filtered_BIC.pickle'
SimiC_weights <- py_load_object(filename =weights_file)
pd <- import("pandas")
AUCs_file <-'./Data/SimiC/Results/5qVsnon5q_1000_L10.01_L20.1_AUCs.pickle'

# filter low R2 targets
unselected_targets <- list()
for (phenotype in phenotypes){
    unselected_targets[[phenotype]] <- SimiC_weights$query_targets[which(SimiC_weights$adjusted_r_squared[[assignment[[phenotype]]]] < 0.7)]

}

pdf(paste0('./Plots/SimiC/5qVsNon5q//Simic_Filtered_R2_targets_hist_5qVsNon5q.pdf'))
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




pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_TF_weigths_5qVsNon5q.pdf'), onefile = TRUE, width=20)
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
    scale_fill_manual(values=c(del5q="#084c61", normal="#db3a34")) + 
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle(drv) 
  print(p)
}
dev.off()


# print stats

# tmp <- SimiC_weights_df
# tmp <- tmp[!tmp$driver %in% TF_2_remove,]
# tmp <- tmp[!tmp$target %in% unselected_targets,]
# tmp <- split(tmp, tmp$driver)
# plotter <- data.frame(driver=NULL, targets_5q = NULL, total_targets=NULL, stringsAsFactors = F)
# for (tf in names(tmp)){
# 	stats <- tmp[[tf]]
# 	stats <- stats[stats$value != 0,]
# 	all_targets <- unique(stats$target)
# 	targets_5q <- sum(unique(stats$target) %in% genes_expanded)
# 	print(paste0(tf, ' :: ', targets_5q, '/', length(all_targets)))
# 	plotter <- rbind(plotter, data.frame(driver=tf, targets_5q = targets_5q, total_targets = length(all_targets)))
# }
# plotter <- plotter[!plotter$driver %in% TF_2_remove,]
# plotter$TF <- 'TF'
# plotter$ratio <- plotter$targets_5q/plotter$total_targets
# pdf(paste0('./Plots/SimiC/5qVsnon5q/Simic_TF_stats_5qVsNon5q.pdf'), onefile = TRUE)
# ggplot(plotter, aes(x=TF, y=ratio, color = ratio, label = driver)) + geom_point() + 
# scale_color_viridis(option="cividis", direction =-1) + theme_bw() + 
# ggrepel::geom_label_repel(data= subset(plotter, ratio > 0.04), size = 4, box.padding = 0.5, point.padding = 0.5, force = 100, segment.size  = 0.2, segment.color = "grey50") +
# theme(axis.text.x = element_text(angle = 45, hjust=1, size=4)) + ggtitle('SimiC TFs')
# dev.off()

pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_Targets_weigths_5qVsNon5q.pdf'), onefile = TRUE, width=20)
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
    scale_fill_manual(values=c(del5q="#084c61", normal="#db3a34"))+ 
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

AUCs <- list()
aucs <- py_load_object(filename =AUCs_file)
for(i in seq_along(phenotypes)){
  AUCs[[i]] <- aucs[[i]]
}

# AUCs <- list()
# for(i in seq_along(phenotypes)){
# 	AUCs[[i]] <- read.table(paste0('./Data/SimiC/Results//5qVsnon5q_1000_L10.01_L20.01_AUCs_', assignment[[i]], '_BIS.csv'), header=T, sep='\t')
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

df <- do.call('rbind', AUCs_by_state)
df$cell_id <- stringr::str_extract(rownames(df), '(?<=\\.)[A-Z0-9-_]+$')
rownames(df) <- df$cell_id

clusters_df <- setNames( seurat_data[, 'Cluster_name', drop=FALSE] , c('cluster_id'))
clusters_df$cell_id <- rownames(clusters_df)
df_w_cluster <- merge(df,clusters_df[, c('cell_id', 'cluster_id')],by="cell_id")

df_auc <- melt(df_w_cluster , id.vars = c('TimePoint','cell_id', 'cluster_id'), variable.name = 'driver', stringsAsFactors =F)



# check the population per phenotype and cluster
tmp <- with(df_w_cluster, table(cluster_id, TimePoint))
tmp <- reshape2::melt(tmp)
 
pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_Population_phenotype_5qVsNon5q.pdf'))
  ggplot(tmp, aes(x=cluster_id, y=value, fill=TimePoint)) + 
  geom_bar(stat='identity' , position='dodge',  color='black') + 
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
  scale_fill_manual(values=c(del5q="#084c61", normal="#db3a34")) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1, size=4))
dev.off()


##### keep the clusters with >10 representatnts in both populations

clusters_2_keep <- as.data.frame.matrix(with(df_auc, table(cluster_id, TimePoint)))
clusters_2_keep$cluster_id <- rownames(clusters_2_keep)
clusters_2_keep <- clusters_2_keep[rowSums(clusters_2_keep > 10) == ncol(clusters_2_keep), 'cluster_id']


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
  pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_Auc_Cluster_',clust,'_5qVsNon5q.pdf'), width = 15, onefile = TRUE)
  plot_counter <- 1
  for (tf in unique(plotter2$driver)){
    assign( paste0('p', plot_counter), 
    ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=TimePoint)) + 
    geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
    scale_fill_manual(values=c(del5q="#084c61", normal="#db3a34")) +
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


pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_Auc_all_cells_5qVsNon5q.pdf'), width = 15, onefile = TRUE)
plot_counter <- 1
plotter2 <- df_auc
for (tf in unique(plotter2$driver)){
	assign( paste0('p', plot_counter), 
	ggplot(plotter2[plotter2$driver ==tf,], aes(x=value, fill=TimePoint)) + 
	geom_density(alpha = 0.6, adjust = 1/8) + theme_classic() + 
  scale_fill_manual(values=c(del5q="#084c61", normal="#db3a34")) +
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


# pdf(paste0('./Plots/SimiC/5qVsnon5q/Simic_Auc_all_cells_perPatient_5qVsNon5q.pdf'), width = 15, onefile = TRUE)
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


# Plot the heatmap of the regulatory dissimilarity score
p <- pheatmap::pheatmap(MinMax_clust,color=plasma, fontsize=5, angle_col =45, cellwidth=40, silent=TRUE)
plot_dims <- get_plot_dims(p)
pdf(paste0('./Plots/SimiC/5qVsNon5q/Simic_HeatMaps_RegDissScore_5qVsNon5q.pdf'), height = plot_dims$height, width = plot_dims$width )
print(p)
dev.off()


# dev.off()
theme_mod <-  theme(axis.text.y=element_text(angle=45, vjust=0.5, size=8), legend.position='none')

pdf('./Plots/SimiC/5qVsNon5q/Enrich_Regulons.pdf')
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
  if( results_BP@readable == TRUE){ plot_list[['BP']] <- clusterProfiler::dotplot(results_BP) + theme_mod}else{NULL}
  if( results_MF@readable == TRUE){ plot_list[['MF']] <- clusterProfiler::dotplot(results_MF) + theme_mod}else{NULL}
  if( results_CC@readable == TRUE){ plot_list[['CC']] <- clusterProfiler::dotplot(results_CC) + theme_mod}else{NULL}
  if( results_KG@readable == TRUE){ plot_list[['KG']] <- clusterProfiler::dotplot(results_KG) + theme_mod}else{NULL}

  print(cowplot::plot_grid(title_theme, 
        cowplot::plot_grid(plotlist=plot_list,
                           ncol=2, labels=c('BP', 'MF', 'CC', 'KEGG')), 
                           nrow=2, rel_heights=c(0.1,0.9)))
}
dev.off()


# all_data <- as.data.frame( seurat_data@assays$RNA@scale.data)
# all_data <- all_data[, colnames(all_data) %in% unique(df_auc[df_auc$TimePoint %in% phenotypes, 'cell_id'])]
# tsne_out <- Rtsne(as.matrix(t(all_data)))
# plotter <-  as.data.frame(tsne_out$Y)
# rownames(plotter) <- colnames(all_data)
# plotter$cell_id <- colnames(all_data)


# genes_2_plot <- c('MYC', 'RUNX3', 'EOMES')
# pdf(paste0('./Plots/SimiC/Simic_tSNE_selection_5qVsNon5q.pdf'), height =8, width=18)
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



