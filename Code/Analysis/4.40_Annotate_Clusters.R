
library(Seurat)
library(ggplot2)

combined_sct_geneset <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/all_seurat_integrated_sct.rds')

annotation_Sofia_02 <- list()
annotation_Sofia_02[['0' ]] <- 'HSC'
annotation_Sofia_02[['1' ]] <- 'GMP/Granul/Mono'
annotation_Sofia_02[['2' ]] <- 'MEP/E-E'
annotation_Sofia_02[['3' ]] <- 'E-L'
annotation_Sofia_02[['4' ]] <- 'CLP/Pro-B'
annotation_Sofia_02[['5' ]] <- 'Unknown'
annotation_Sofia_02[['6' ]] <- 'DC/LMPP/T'
annotation_Sofia_02[['7' ]] <- 'Unknown'
annotation_Sofia_02[['8' ]] <- 'Baso'
annotation_Sofia_02[['9' ]] <- 'CFU-MK'
annotation_Sofia_02[['10']] <- 'Unknown'



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



annotation_Sofia_08 <- list()
annotation_Sofia_08[['0' ]] <- 'HSC'
annotation_Sofia_08[['1' ]] <- 'HSC'
annotation_Sofia_08[['2' ]] <- 'HSC'
annotation_Sofia_08[['3' ]] <- 'EarlyErythroid'
annotation_Sofia_08[['4' ]] <- 'Pro-B'
annotation_Sofia_08[['5' ]] <- 'LMPP'
annotation_Sofia_08[['6' ]] <- 'Mono/GMP'
annotation_Sofia_08[['7' ]] <- 'LMPP/CLP'
annotation_Sofia_08[['8' ]] <- 'GMP'
annotation_Sofia_08[['9' ]] <- 'EarlyErythroid'
annotation_Sofia_08[['10' ]] <- 'EarlyL'
annotation_Sofia_08[['11']] <- 'GRAN/GMP'
annotation_Sofia_08[['12']] <- 'EarlyErythroid'
annotation_Sofia_08[['13']] <- 'EarlyErythroid'
annotation_Sofia_08[['14']] <- 'EarlyErythroid'
annotation_Sofia_08[['15']] <- 'CLP/ProB'
annotation_Sofia_08[['16']] <- 'MEP'
annotation_Sofia_08[['17']] <- 'MEO/GMP'
annotation_Sofia_08[['18']] <- 'EarlyErythroid'
annotation_Sofia_08[['19']] <- 'ProB'
annotation_Sofia_08[['20']] <- 'BASO'
annotation_Sofia_08[['21']] <- 'LMPP'
annotation_Sofia_08[['22']] <- 'T'
annotation_Sofia_08[['23']] <- 'Unknown'
annotation_Sofia_08[['24']] <- 'MONO'
annotation_Sofia_08[['25']] <- 'Unknown'
annotation_Sofia_08[['26']] <- 'DC'
annotation_Sofia_08[['27']] <- 'EL'
annotation_Sofia_08[['28']] <- 'Unknown'
annotation_Sofia_08[['29']] <- 'Unknown'



annotation_Aintz <- list()
annotation_Aintz[['0' ]] <- 'HSC/LMPP'
annotation_Aintz[['1' ]] <- 'HSC'
annotation_Aintz[['2' ]] <- 'Late-erythroid'
annotation_Aintz[['3' ]] <- 'Early-erythroid'
annotation_Aintz[['4' ]] <- 'MEP/Early-megakaryocyte'
annotation_Aintz[['5' ]] <- 'T_cell'
annotation_Aintz[['6' ]] <- 'GMP/Early-monocyte'
annotation_Aintz[['7' ]] <- 'ProB'
annotation_Aintz[['8' ]] <- 'Unknown'
annotation_Aintz[['9' ]] <- 'CMP'
annotation_Aintz[['10']] <- 'CMP'
annotation_Aintz[['11']] <- 'GMP'
annotation_Aintz[['12']] <- 'Unknown'
annotation_Aintz[['13']] <- 'Unknown'
annotation_Aintz[['14']] <- 'CLP'
annotation_Aintz[['15']] <- 'Basophil'
annotation_Aintz[['16']] <- 'Unknown'



annotate_clusters <- function(data, annotation_list){
	cluster_col <- grep('Cluster', colnames(data))
	cluster_names <- c()
	pb <- progress::progress_bar$new(total = nrow(data), format="[:bar] :percent eta: :eta")
	for (cl in as.character(data[[cluster_col]])){
		cluster_names <- c(cluster_names , annotation_list[[as.character(cl)]])
		pb$tick()
	}
	data$ClusterName <- cluster_names
	return(data)
	
}



combined_sct_geneset <- SetIdent(combined_sct_geneset, value = "integrated_snn_res.0.4")
cell_typist <- t(as.data.frame(combined_sct_geneset@assays$RNA@counts))
# write.csv(cell_typist, file = '/home/sevastopol/data/gserranos/MDS/Data/cell_type_counts_4_celltypist.csv')

celtypist <- read.csv('./Data/predicted_labels_celtypist.csv')
rownames(celtypist) <- celtypist$X


# markers <- readRDS('./Data/all_seurat_integrated_sct_markers.rds')
# tmp <- markers[["integrated_snn_res.0.4"]]
# top_10_markers <- c()
# for (cluster in unique(tmp$cluster)){
# 	tmp2 <- tmp[tmp$cluster == cluster,]
# 	tmp2 <- tmp2[order(tmp2$avg_logFC),]
# 	top_10_markers <- c(top_10_markers, tmp2[1:10,'gene'])
# }
# pdf('./Plots/Test.pdf')
# DoHeatmap(combined_sct_geneset, features = top_10_markers) + NoLegend()
# dev.off()








get_umap_ann <- function(plotter, title){
	n_clusters <- length(unique(plotter$ClusterName))
	colors <-  colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
	p <- ggplot(plotter, aes(x= UMAP_1, y = UMAP_2, color =ClusterName)) + geom_point(alpha=0.6, size = 1) + theme_void() + 
	scale_color_manual(values= colors) + theme(legend.position='none') + ggtitle(title) + 
	guides(colour = guide_legend(override.aes = list(size=5)))
	return(p)
}

coords <- as.data.frame(combined_sct_geneset@reductions$umap@cell.embeddings)

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.2), c('Cluster'))
clustering <- annotate_clusters(clustering, annotation_Sofia_02)
tmp <- merge(coords, clustering, by = 0)
clusters_02 <- unique(tmp$ClusterName)
p02 <- get_umap_ann(tmp, 'Sofia_02') + theme(legend.position='right')

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.4), c('Cluster'))
clustering <- annotate_clusters(clustering, annotation_Sofia_04)
tmp <- merge(coords, clustering, by = 0)
clusters_04 <- unique(tmp$ClusterName)
p04 <- get_umap_ann(tmp, 'Sofia_04')+ theme(legend.position='right')

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
clustering <- annotate_clusters(clustering, annotation_Sofia_08)
tmp <- merge(coords, clustering, by = 0)
p08 <- get_umap_ann(tmp, 'Sofia_08')+ theme(legend.position='bottom')

colors <-  colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(union(clusters_04, clusters_02)))
names(colors) <-  unique(union(clusters_04, clusters_02))

pdf('./Plots/Sofia_annotation.pdf')
cowplot::plot_grid( p02 + scale_color_manual(values= colors) + theme(legend.position='none'), 
					p04 + scale_color_manual(values= colors) + theme(legend.position='bottom'), ncol=1)
dev.off()




pdf('./Plots/Integrated/Celltypist_ann.pdf')

coords <- as.data.frame(combined_sct_geneset@reductions$umap@cell.embeddings)
tmp <- merge(coords, celtypist, by.x = 0, by.y='X')
names(tmp)[names(tmp) == 'majority_voting'] <- 'Cluster_name'
print(get_umap_ann(tmp, 'Celltypist'))

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Sofia[[cl]])
}
clustering$Cluster_name <- cluster_names
tmp <- merge(coords, clustering, by = 0)
print(get_umap_ann(tmp, 'Sofia'))

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.4), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Aintz[[cl]])
}
clustering$Cluster_name <- cluster_names
tmp <- merge(coords, clustering, by = 0)
print(get_umap_ann(tmp, 'Aintzane'))

dev.off()

# harmonize annotations


harmonize_celltypist <- list()
harmonize_celltypist[['Early erythroid']]                              <-  'Early_Ery'
harmonize_celltypist[['HSC/MPP']]                                      <-  'HSC/LMPP'
harmonize_celltypist[['Early MK']]                                     <-  'MEP/EM'
harmonize_celltypist[['Megakaryocyte precursor']]                      <-  'MEP/EM'
harmonize_celltypist[['Mast cells']]                                   <-  'Gran/Baso'
harmonize_celltypist[['CMP']]                                          <-  'CMP'
harmonize_celltypist[['Mid erythroid']]                                <-  'Early_Ery'
harmonize_celltypist[['Pro-B cells']]                                  <-  'ProB'
harmonize_celltypist[['pDC']]                                          <-  'DC'
harmonize_celltypist[['Neutrophil-myeloid progenitor']]                <-  'GMP/Gran'
harmonize_celltypist[['Tcm/Naive helper T cells']]                     <-  'T'
harmonize_celltypist[['Tem/Trm cytotoxic T cells']]                    <-  'T'
harmonize_celltypist[['ELP']]                                          <-  'CLP'
harmonize_celltypist[['Classical monocytes']]                          <-  'Mono'
harmonize_celltypist[['Megakaryocyte-erythroid-mast cell progenitor']] <-  'MEP/EM'
harmonize_celltypist[['MEMP']]                                         <-  'MEP/EM'
harmonize_celltypist[['DC precursor']]                                 <-  'DC'
harmonize_celltypist[['GMP']]                                          <-  'GMP'
harmonize_celltypist[['Promyelocytes']]                                <-  'GMP/Gran'
harmonize_celltypist[['Megakaryocytes/platelets']]                     <-  'Late_Ery'
harmonize_celltypist[['Large pre-B cells']]                            <-  'ProB'
harmonize_celltypist[['Double-positive thymocytes']]                   <-  'T'




annotation_Sofia_harmonized <- list()
annotation_Sofia_harmonized[['0' ]] <- 'HSC/LMPP'
annotation_Sofia_harmonized[['1' ]] <- 'HSC/LMPP'
annotation_Sofia_harmonized[['2' ]] <- 'HSC/LMPP'
annotation_Sofia_harmonized[['3' ]] <- 'Early_Ery'
annotation_Sofia_harmonized[['4' ]] <- 'ProB'
annotation_Sofia_harmonized[['5' ]] <- 'HSC/LMPP'
annotation_Sofia_harmonized[['6' ]] <- 'GMP/Mono'
annotation_Sofia_harmonized[['7' ]] <- 'CLP'
annotation_Sofia_harmonized[['8' ]] <- 'GMP'
annotation_Sofia_harmonized[['9' ]] <- 'Early_Ery'
annotation_Sofia_harmonized[['10' ]] <- 'Early_Ery' #Late_Ery
annotation_Sofia_harmonized[['11']] <- 'GMP/Gran'
annotation_Sofia_harmonized[['12']] <- 'Early_Ery'
annotation_Sofia_harmonized[['13']] <- 'Early_Ery'
annotation_Sofia_harmonized[['14']] <- 'Early_Ery'
annotation_Sofia_harmonized[['15']] <- 'ProB'
annotation_Sofia_harmonized[['16']] <- 'MEP/EM'
annotation_Sofia_harmonized[['17']] <- 'CMP'
annotation_Sofia_harmonized[['18']] <- 'Early_Ery'
annotation_Sofia_harmonized[['19']] <- 'ProB'
annotation_Sofia_harmonized[['20']] <- 'Gran/Baso'
annotation_Sofia_harmonized[['21']] <- 'HSC/LMPP'
annotation_Sofia_harmonized[['22']] <- 'T'
annotation_Sofia_harmonized[['23']] <- 'Unknown'
annotation_Sofia_harmonized[['24']] <- 'Mono'
annotation_Sofia_harmonized[['25']] <- 'Unknown'
annotation_Sofia_harmonized[['26']] <- 'DC'
annotation_Sofia_harmonized[['27']] <- 'Late_Ery'
annotation_Sofia_harmonized[['28']] <- 'Unknown'
annotation_Sofia_harmonized[['29']] <- 'Unknown'

annotation_Aintz_harmonized <- list()
annotation_Aintz_harmonized[['0' ]] <- 'HSC/LMPP'
annotation_Aintz_harmonized[['1' ]] <- 'HSC/LMPP'
annotation_Aintz_harmonized[['2' ]] <- 'Late_Ery'
annotation_Aintz_harmonized[['3' ]] <- 'Early_Ery'
annotation_Aintz_harmonized[['4' ]] <- 'MEP/EM'
annotation_Aintz_harmonized[['5' ]] <- 'T'
annotation_Aintz_harmonized[['6' ]] <- 'GMP/Mono'
annotation_Aintz_harmonized[['7' ]] <- 'ProB'
annotation_Aintz_harmonized[['8' ]] <- 'Unknown'
annotation_Aintz_harmonized[['9' ]] <- 'CMP'
annotation_Aintz_harmonized[['10']] <- 'CMP'
annotation_Aintz_harmonized[['11']] <- 'GMP'
annotation_Aintz_harmonized[['12' ]] <- 'Unknown'
annotation_Aintz_harmonized[['13' ]] <- 'Unknown'
annotation_Aintz_harmonized[['14']] <- 'CLP'
annotation_Aintz_harmonized[['15']] <- 'Gran/Baso'
annotation_Aintz_harmonized[['16' ]] <- 'Unknown'



get_cell_dist <- function(data, colors){
	p <- ggplot(data, aes(x=ClusterName, fill=ClusterName)) + geom_bar() + 
	theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position='none') +
	scale_fill_manual(values=colors)
	return(p)
}


all_possible <- unique(c(paste(unlist(harmonize_celltypist)),  paste(unlist(annotation_Aintz_harmonized)),  paste(unlist(annotation_Sofia_harmonized))))
colors <- ggthemes::tableau_color_pal('Classic 20')(length(all_possible))
names(colors) <- all_possible

harmonized_names <- c()
for (cl in celtypist$majority_voting){
	harmonized_names <- c(harmonized_names, harmonize_celltypist[[cl]])
}
celtypist$harmonized <-harmonized_names

coords <- as.data.frame(combined_sct_geneset@reductions$umap@cell.embeddings)
tmp <- merge(coords, celtypist, by.x = 0, by.y='X')
names(tmp)[names(tmp) == 'harmonized'] <- 'ClusterName'
umap_celltypist <- get_umap_ann(tmp, 'Celltypist') + scale_color_manual(values = colors)
dist_celltypist <- get_cell_dist(tmp, colors)

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Sofia_harmonized[[cl]])
}
clustering$ClusterName <- cluster_names


tmp <- merge(coords, clustering, by = 0)
umap_sofia <-  get_umap_ann(tmp, 'Sofia') + scale_color_manual(values = colors)
dist_sofia <- get_cell_dist(tmp, colors)

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.4), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Aintz_harmonized[[cl]])
}
clustering$ClusterName <- cluster_names
tmp <- merge(coords, clustering, by = 0)
umap_aintzane <- get_umap_ann(tmp, 'Aintzane') + scale_color_manual(values = colors)
dist_aintzane <- get_cell_dist(tmp, colors)

pdf('./Plots/Integrated/All_annotations_ann_harmonized.pdf')
legend <- cowplot::get_legend(umap_sofia+theme(legend.position='bottom'))
cowplot::plot_grid(
	umap_celltypist, umap_sofia, umap_aintzane, 
	NULL, legend, NULL,
	dist_celltypist, dist_sofia, dist_aintzane,
	ncol=3, nrow=3, rel_heights=c(0.4, 0.1, 0.4)
)

dev.off()




clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
tmp <- merge(coords, clustering, by = 0)
Numeric_cluster <- ggplot(tmp, aes(x= UMAP_1, y = UMAP_2, color =Cluster)) + geom_point(alpha=0.6) + theme_void() + 
					theme(legend.position='none') + ggtitle('Numeric_cluster') + 
					scale_color_manual(values=sample(destiny::cube_helix(length(unique(tmp$Cluster)))))

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Sofia[[cl]])
}
clustering$ClusterName <- cluster_names

tmp <- merge(coords, clustering, by = 0)
colors <- ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$Cluster_name)))
Sofia_Original <- get_umap_ann(tmp, 'Sofia_Original') + scale_color_manual(values = colors)


clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.8), c('Cluster'))
cluster_names <- c()
for (cl in clustering$Cluster){
	cluster_names <- c(cluster_names , annotation_Sofia_harmonized[[cl]])
}
clustering$ClusterName <- cluster_names


tmp <- merge(coords, clustering, by = 0)
colors <- ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$Cluster_name)))
Sofia_Harmonized <- get_umap_ann(tmp, 'Sofia_Harmonized') + scale_color_manual(values = colors)


pdf('./Plots/Integrated/Annotation_sofia.pdf')
cowplot::plot_grid(
	Numeric_cluster, Sofia_Original, Sofia_Harmonized
)
dev.off()



clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.2), c('Cluster'))
tmp <- merge(coords, clustering, by = 0)
Numeric_cluster_2 <- ggplot(tmp, aes(x= UMAP_1, y = UMAP_2, color =Cluster)) + geom_point(alpha=0.6) + theme_void() + 
					theme(legend.position='right') + ggtitle('Numeric_cluster') + 
					scale_color_manual(values= ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$Cluster))))

clustering <- setNames(as.data.frame(combined_sct_geneset$integrated_snn_res.0.4), c('Cluster'))
tmp <- merge(coords, clustering, by = 0)
Numeric_cluster_4 <- ggplot(tmp, aes(x= UMAP_1, y = UMAP_2, color =Cluster)) + geom_point(alpha=0.6) + theme_void() + 
					theme(legend.position='right') + ggtitle('Numeric_cluster') + 
					scale_color_manual(values= ggthemes::tableau_color_pal('Classic 20')(length(unique(tmp$Cluster))))

pdf('./Plots/Integrated/Numeric_Clusters_02_04.pdf')
cowplot::plot_grid(
	Numeric_cluster_2, Numeric_cluster_4, ncol=1
	)
dev.off()


pdf('./Plots/Integrated/TestEnrichr.pdf')
DEenrichRPlot(object =combined_sct_geneset, ident.1 = "0",ident.2 = "1",enrich.database = "BP",max.genes = 50)
dev.off()



markers <- readRDS('./Data/all_seurat_integrated_sct_markers.rds')
max.genes = 50
for(cluster in sort(unique(combined_sct_geneset$integrated_snn_res.0.4))){
	markers_clt <- markers[[0.2]]
	all.markers <- markers_clt[markers_clt$cluster == cluster,]


	pos.markers <- all.markers[all.markers[, 2] > 0.25 & all.markers[, 1] < 0.05, , drop = FALSE]
	pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
	pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = 'GO_Biological_Process_2013')
	pos.er <- do.call(what = cbind, args = pos.er)
	pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
	pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
	pos.er <- pos.er[1:num.pathway, ]
	pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
	gene.list <- list(pos = pos.er)

}


# combined_sct_geneset$tmp <- ifelse(combined_sct_geneset$integrated_snn_res.0.4 == cluster, 1, 0)
# Idents(combined_sct_geneset) <- 'tmp'
# all.markers <- FindMarkers(
#     object = combined_sct_geneset,
#     ident.1 = 1,
#     ident.2 = 0,
#     only.pos = FALSE,
#     logfc.threshold = 0.25,
#     test.use = 'wilcox',
#     assay = DefaultAssay(combined_sct_geneset)
#   )

# pos.markers <- all.markers[all.markers[, 2] > 0.25 & all.markers[, 1] < 0.05, , drop = FALSE]

  if(nrow(pos.markers) == 0){
    message("No positive markers pass the logfc.thershold")
    pos.er <- c()
  }

  else{
#   pos.markers.list <- rownames(x = pos.markers)[1:min(max.genes, nrow(x = pos.markers))]
#   pos.er <- enrichR::enrichr(genes = pos.markers.list, databases = enrich.database)
#   pos.er <- do.call(what = cbind, args = pos.er)
#   pos.er$log10pval <- -log10(x = pos.er[, paste(enrich.database, sep = ".", "P.value")])
#   pos.er$term <- pos.er[, paste(enrich.database, sep = ".", "Term")]
#   pos.er <- pos.er[1:num.pathway, ]
#   pos.er$term <- factor(x = pos.er$term, levels = pos.er$term[order(pos.er$log10pval)])
#   gene.list <- list(pos = pos.er)
  }



p2 <- ggplot(data = neg.er, aes_string(x = "term", y = "log10pval")) +
        geom_bar(stat = "identity", fill = "indianred2") +
        coord_flip() + xlab("Pathway") +
        scale_fill_manual(values = cols, drop = FALSE) +
        ylab("-log10(pval)") +
        ggtitle(paste(enrich.database, ident.1, sep = "_", "negative markers")) +
        theme_classic() +
        geom_text(aes_string(label = "term", y = 0),
                  size = 5,
                  color = "black",
                  position = position_dodge(1),
                  hjust = 0)+
        theme(axis.title.y= element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
      p <- p2
	


 p <- ggplot(data = pos.er, aes_string(x = "term", y = "log10pval")) +
    geom_bar(stat = "identity", fill = "dodgerblue") +
    coord_flip() + xlab("Pathway") +
    scale_fill_manual(values = cols, drop = FALSE) +
    ylab("-log10(pval)") +
    ggtitle(paste(enrich.database, ident.1, sep = "_", "positive markers")) +
    theme_classic() +
    geom_text(aes_string(label = "term", y = 0),
              size = 5,
              color = "black",
              position = position_dodge(1),
              hjust = 0)+
    theme(axis.title.y= element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())




celtypist_probs <- read.csv('./Data/predicted_labels_celtypist_probs.csv')

test <- merge(coords, celtypist[, 1:2], by=0)
test <- merge(test, celtypist_probs[, 1:2], by='X')
library(ggplot2)

pdf('./Plots/Test_contour.pdf')
ggplot(test[test$predicted_labels == 'B cells', ], aes(UMAP_1, UMAP_2, z = B.cells))+  stat_contour()
dev.off()






