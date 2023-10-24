

library(Seurat)
library(ggplot2)


data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/all_seurat_Post_integrated_sct_sub.rds')




marker_list <- list()
marker_list[['hsc']]       <-c('CRHBP', 'HOPX', 'KYT', 'CD34', 'AVP', 'PRSS2','MLLT3', 'IDS', 'BTS2')  
marker_list[['lmpp']]      <-c('PTPRC', 'FLT3', 'PROM1', 'SATB1', 'ACY3','IGJ', 'LSP1', 'NEGR1', 'SPINK2')  
marker_list[['gmp']]       <-c('CSF3R', 'CTSG', 'PRTN3', 'MPO','CFD', 'CTSG', 'CSTA', 'CST7')  
marker_list[['granul']]    <-c('ELANE', 'AZU1', 'CEBPA', 'CEBPE', 'CST7','RNASE2')  
marker_list[['mono']]      <-c('LYZ', 'CSTA', 'FCER1G', 'TLR1', 'IFNGR2')  
marker_list[['dc']]        <-c('IRF8', 'IRF7', 'IL3RA', 'CLEC4', 'IRF4', 'ZEB2', 'KLF4', 'AXL')  
marker_list[['t']]         <-c('JCHAIN', 'IKZF1', 'CYTH1', 'GZMA', 'SH2D1A', 'IL32', 'MT1E', 'CXCR3')  
marker_list[['clp']]       <-c('IL7R', 'DNTT', 'CYGB', 'LTB', 'IGJ', 'DNTT', 'ADA')  
marker_list[['prob']]      <-c('VPREB1', 'EBF1', 'CD79A', 'CD79B', 'TCL1A', 'IRF4', 'CCDC42B', 'FCRLA', 'IGLL5')  
marker_list[['mep']]       <-c('NFE2', 'HFS1', 'TAL1', 'FCER1A', 'PBX1', 'PLEK', 'DAD1', 'IGSF10')  
marker_list[['mk']]        <-c('PBX1', 'MPL', 'VWF', 'FLI1', 'ITGA22B', 'GP1BA', 'CMTM5', 'PF4', 'GP9', 'CLEC1B', 'PPBP')  
marker_list[['ery_early']] <- c('CNRIP1', 'SLC40A1', 'PKIG', 'PDZD8', 'FCER1A', 'CSF2RB', 'EMP3', 'HBD')    
marker_list[['ery_late']]  <-c('APOC1', 'BLVRB', 'APOE', 'FAM178B', 'CA1', 'CNRIP1', 'AHSP', 'KCNH2', 'TFR2', 'HBB', 'KLF1', 'S100A6')  
marker_list[['baso']]      <-c('RUNX1', 'HDC', 'MS4A2', 'MS4A3', 'TPSAB1')



for (celltype in names(marker_list)){
	message(celltype)
	features <- list(c(marker_list[[celltype]]))
	tmp <- AddModuleScore(data, features=features, name = paste0(celltype, '_score'))
	coords <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.1.2',paste0(celltype, '_score1')))
	coords <- coords[order(coords[[paste0(celltype, '_score1')]]),]
	plot <- ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= get(paste0(celltype, '_score1')))) + 
			geom_point(size=1, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = celltype) + ggtitle(paste0(celltype, '_score1')) + 
			ggdark::dark_theme_grey() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank())
	pdf(paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2/UMAP_SubCluster_',celltype,'_Score.pdf'))
		print(plot)
		print(plot + 
			geom_point(data = transform(coords, integrated_snn_res.1.2 = NULL), colour = "grey85", size=0.2, alpha=0.3) +
			geom_point(size=1, alpha=0.6) + facet_wrap(~integrated_snn_res.1.2))
	dev.off()
}




get_all_markers <- function(){
	all_markers <- read.table('./Data/Annotation/E-HCAD-6-marker-genes-files/E-HCAD-6.marker_genes_inferred_cell_type_-_ontology_labels.tsv', header=TRUE, sep='\t')
	all_markers$cluster <- toupper(all_markers$cluster)
	all_selected_markers <- list()
	annotation <- setNames(read.table('./Data/Annotation/Homo_sapiens.GRCh38.Gene_name_ID.tsv', sep='\t'), c('genes', 'gene_name'))
	for (marker in unique(all_markers$cluster)){
		tmp <- all_markers[all_markers$cluster == toupper(marker) & all_markers$pvals_adj < 0.05,]
		# pdf(paste0('./Data/Annotation/',  gsub(" ", "_", cell_type), '_logfc_hist.pdf'))
		# 	hist(tmp$logfoldchanges, breaks=100, col='black', 
		# 	main=paste0(cell_type, ' Nº of markers :: ', nrow(tmp)), xlab='logFC', ylab='Frequency')
		# dev.off()
		tmp <- tmp[tmp$logfoldchanges >1,]
		tmp <- merge(tmp, annotation,  by='genes')
		all_selected_markers[[marker]] <- tmp$gene_name
	}
	return(all_selected_markers)
}



markers_other <- get_all_markers()

for (celltype in names(markers_other)){
	message(celltype)
	features <- list(c(markers_other[[celltype]]))
	celltype <- gsub(' ', '_', celltype)
	tmp <- AddModuleScore(data, features=features, name = paste0(celltype, '_score'))
	coords <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.1.2', paste0(celltype, '_score1')))
	coords <- coords[order(coords[[paste0(celltype, '_score1')]]),]
	plot <- ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= get(paste0(celltype, '_score1')))) + 
			geom_point(size=1, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = celltype) + ggtitle(paste0(celltype, '_score1')) + 
			ggdark::dark_theme_grey() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank())
	pdf(paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2/2_UMAP_SubCluster_',celltype,'_Score_markers_other.pdf'))
		print(plot)
		print(plot + 
			geom_point(data = transform(coords, integrated_snn_res.1.2 = NULL), colour = "grey85", size=0.2, alpha=0.3) +
			geom_point(size=1, alpha=0.6) + facet_wrap(~integrated_snn_res.1.2))
	dev.off()
}




results_all <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WriteXLS::WriteXLS(results_all, ExcelFileName=paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2/DE_Vs_ALL.xlsx'), SheetNames = 'DE_Vs_ALL',  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

DefaultAssay(data) <- 'SCT'
Idents(data) <- 'integrated_snn_res.1.2'
list_DE <- list()
for (cl1 in sort(unique(data$integrated_snn_res.1.2))){
	for (cl2 in sort(unique(data$integrated_snn_res.1.2))){
		if (as.numeric(cl1) >= as.numeric(cl2)){
			next
		}
		else{
			print(paste0(cl1, ' vs ', cl2))
			list_DE[[paste0(cl1, 'vs', cl2)]]<- FindMarkers(
			object = data,
			ident.1 = cl1,
			ident.2 = cl2,
			assay = "SCT",
			recorrect_umi = FALSE)
		}
	}
}
WriteXLS::WriteXLS(list_DE, ExcelFileName=paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2/DE_per_Cluster.xlsx'), SheetNames = names(list_DE),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)




# SUBCLUSTER
cluster_2_subCluster <- c(4, 7, 8, 10, 14, 18)
Idents(data) <- 'integrated_snn_res.1.2'
for (cluster in cluster_2_subCluster){
	print(cluster)
	data <- FindSubCluster(data, cluster, 'integrated_snn', subcluster.name = "integrated_snn_res.1.2_sub", resolution = 0.4)
	Idents(data) <- "integrated_snn_res.1.2_sub"
}



for (celltype in names(marker_list)){
	message(celltype)
	features <- list(c(marker_list[[celltype]]))
	tmp <- AddModuleScore(data, features=features, name = paste0(celltype, '_score'))
	coords <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.1.2_sub',paste0(celltype, '_score1')))
	coords <- coords[order(coords[[paste0(celltype, '_score1')]]),]
	plot <- ggplot(coords, aes(x=UMAP_1, y= UMAP_2, color= get(paste0(celltype, '_score1')))) + 
			geom_point(size=1, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = celltype) + ggtitle(paste0(celltype, '_score1')) + 
			ggdark::dark_theme_grey() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank())
	pdf(paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/UMAP_SubCluster_',celltype,'_Score.pdf'))
		print(plot)
		print(plot + 
			geom_point(data = transform(coords, integrated_snn_res.1.2_sub = NULL), colour = "grey85", size=0.2, alpha=0.3) +
			geom_point(size=1, alpha=0.6) + facet_wrap(~integrated_snn_res.1.2_sub))
	dev.off()
}

tmp <- data
for (celltype in names(marker_list)){
	message(celltype)
	features <- list(c(marker_list[[celltype]]))
	tmp <- AddModuleScore(tmp, features=features, name = paste0(celltype, '_score'))
}
coords <- FetchData(tmp, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.1.2_sub', paste0(names(marker_list), '_score1')))
plotter_all_Scores <- reshape2::melt(coords, id=c('UMAP_1', 'UMAP_2', 'integrated_snn_res.1.2_sub'))
plotter_all_Scores <- plotter_all_Scores[order(plotter_all_Scores$value),]
plotter_all_Scores <- ggplot(plotter_all_Scores, aes(x=UMAP_1, y= UMAP_2, color= value)) + 
			geom_point(size=0.4, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = celltype) + ggtitle(paste0(celltype, '_score1')) + 
			ggdark::dark_theme_grey() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank())

png(paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/UMAP_SubCluster_All_CT_Score.png'))
print(plotter_all_Scores + facet_wrap(~variable))
dev.off()


saveRDS(data, '/home/tereshkova/data/gserranos/MDS/Data/all_seurat_Post_integrated_sct_sub.rds')


library(future)
plan("multiprocess", workers = 32)
plan()

Idents(data) <- 'integrated_snn_res.1.2_sub'
results_all <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WriteXLS::WriteXLS(results_all, ExcelFileName=paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/DE_Vs_ALL.xlsx'), SheetNames = 'DE_Vs_ALL',  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


DefaultAssay(data) <- 'SCT'
Idents(data) <- 'integrated_snn_res.1.2_sub'
list_DE <- list()
for (cl1 in sort(unique(data$integrated_snn_res.1.2_sub))){
	for (cl2 in sort(unique(data$integrated_snn_res.1.2_sub))){
		if (as.numeric(sub('_', '.', cl1)) >= as.numeric(sub('_', '.', cl2))){
			next
		}
		else{
			print(paste0(cl1, ' vs ', cl2))
			list_DE[[paste0(cl1, 'vs', cl2)]]<- FindMarkers(
			object = data,
			ident.1 = cl1,
			ident.2 = cl2,
			assay = "SCT",
			recorrect_umi = FALSE)
		}
	}
}
WriteXLS::WriteXLS(list_DE, ExcelFileName=paste0('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/DE_per_Cluster.xlsx'), SheetNames = names(list_DE),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)




filter_and_annotate <- function(data_obj, annotation, resolution){
    clusters <- setNames(as.data.frame(data_obj[[resolution]]), c('Cluster'))
    clusters <- clusters[paste(clusters$Cluster) %in% names(annotation),, drop=FALSE]
    clusters$Cluster_names <- apply(clusters, 1, function(x){ annotation[as.character(x[['Cluster']])][[1]] })
    Idents(data_obj) <- resolution
    data_obj <- subset(data_obj, idents= c(names(annotation)))
    data_obj$Cluster_names <- clusters$Cluster_names
    return(data_obj)
}



annotation_Nerea_1.2_sub_POST <- list()
annotation_Nerea_1.2_sub_POST[["0"]]    <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["1"]]    <- 'EarlyErythroid'
annotation_Nerea_1.2_sub_POST[["2"]]    <- 'LMPP'
annotation_Nerea_1.2_sub_POST[["3"]]    <- 'CLP'
annotation_Nerea_1.2_sub_POST[["4_0"]]  <- 'MEP'
annotation_Nerea_1.2_sub_POST[["4_1"]]  <- 'MEP'
annotation_Nerea_1.2_sub_POST[["4_2"]]  <- 'MK_Prog'
annotation_Nerea_1.2_sub_POST[["4_3"]]  <- 'MEP'
annotation_Nerea_1.2_sub_POST[["5"]]    <- 'HSC'
annotation_Nerea_1.2_sub_POST[["6"]]    <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["7_0"]]  <- 'Monocytes'
annotation_Nerea_1.2_sub_POST[["7_1"]]  <- 'Granulocyte'
annotation_Nerea_1.2_sub_POST[["7_2"]]  <- 'Monocytes'
annotation_Nerea_1.2_sub_POST[["8_0"]]  <- 'EarlyErythroid'
annotation_Nerea_1.2_sub_POST[["8_1"]]  <- 'LateErythroid'
# annotation_Nerea_1.2_sub_POST[["8_2"]]  <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["8_3"]]  <- 'LateErythroid'
annotation_Nerea_1.2_sub_POST[["9"]]    <- 'GMP'
annotation_Nerea_1.2_sub_POST[["10_0"]] <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["10_1"]] <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["10_2"]] <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["10_3"]] <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["10_4"]] <- 'LMPP'
annotation_Nerea_1.2_sub_POST[["11"]]   <- 'DendriticCell'
annotation_Nerea_1.2_sub_POST[["12"]]   <- 'CLP'
annotation_Nerea_1.2_sub_POST[["13"]]   <- 'LateErythroid'
annotation_Nerea_1.2_sub_POST[["14_0"]] <- 'LMPP'
annotation_Nerea_1.2_sub_POST[["14_1"]] <- 'CLP'
annotation_Nerea_1.2_sub_POST[["14_2"]] <- 'CLP'
annotation_Nerea_1.2_sub_POST[["14_3"]] <- 'CLP'
# annotation_Nerea_1.2_sub_POST[["14_4"]] <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["15"]]   <- 'pro-B'
annotation_Nerea_1.2_sub_POST[["16"]]   <- 'EarlyErythroid'
annotation_Nerea_1.2_sub_POST[["17"]]   <- 'EarlyErythroid'
# annotation_Nerea_1.2_sub_POST[["18_0"]] <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["18_1"]] <- 'Basophil'
# annotation_Nerea_1.2_sub_POST[["18_2"]] <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["19"]]   <- 'LateErythroid'
# annotation_Nerea_1.2_sub_POST[["20"]]   <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["21"]]   <- 'Monocytes'
# annotation_Nerea_1.2_sub_POST[["22"]]   <- 'ADIOS'
annotation_Nerea_1.2_sub_POST[["23"]]   <- 'DendriticCell'
annotation_Nerea_1.2_sub_POST[["24"]]   <- 'LMPP'






Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#006d77')
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog', '???', 'Platelets')

data_ann <- filter_and_annotate(data, annotation_Nerea_1.2_sub_POST, 'integrated_snn_res.1.2_sub')

plotter <- FetchData(data_ann, vars = c('UMAP_1', 'UMAP_2','integrated_snn_res.1.2_sub', 'Cluster_names', 'Sample'))
pdf('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/Annotated_UMAP.pdf')
p <- ggplot(plotter, aes(x=UMAP_1, y= UMAP_2, color= Cluster_names)) + 
			geom_point(size=0.4, alpha=0.9) +
			scale_color_manual(values=Cluster_colors)  + 
			# ggdark::dark_theme_grey() + 
			theme_classic() +
			theme(legend.position="bottom", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank()) +
			guides(color = guide_legend(nrow=3, byrow=TRUE, override.aes = list(size=3, alpha=0.9)))
print(p)
print(p + facet_wrap(~integrated_snn_res.1.2_sub))
print(p + facet_wrap(~Sample))
dev.off()


table(data_ann$Cluster_names, data_ann$integrated_snn_res.1.2_sub)



get_plot_4_cluster <- function(dataplotter, cluster, resolution='integrated_snn_res.1.2_sub'){
	dataplotter <- dataplotter[dataplotter[[resolution]] == cluster,]
	p <- ggplot(dataplotter, aes(x=UMAP_1, y= UMAP_2, color= value)) + 
			geom_point(size=0.4, alpha=0.6) +
			scale_color_distiller(palette = "Spectral", direction = -1, name = 'Score')  + ggtitle(paste0('Cluster ', cluster)) +
			ggdark::dark_theme_grey() + theme(legend.position="right", panel.background = element_blank(), 
			axis.ticks=element_blank(), axis.text=element_blank()) + facet_wrap(~variable)
	return(p)
}

pdf('/mnt/md0/gserranos/MDS/Plots/Integrated/Post/Markers/integrated_snn_res.1.2_sub/Check_ann_cluster.pdf')
print(get_plot_4_cluster(plotter_all_Scores, '2'))
print(get_plot_4_cluster(plotter_all_Scores, '0'))
print(get_plot_4_cluster(plotter_all_Scores, '15'))
print(get_plot_4_cluster(plotter_all_Scores, '18_2'))
print(get_plot_4_cluster(plotter_all_Scores, '3'))
print(get_plot_4_cluster(plotter_all_Scores, '8_2'))
print(get_plot_4_cluster(plotter_all_Scores, '20'))
dev.off()

saveRDS(data, '/home/tereshkova/data/gserranos/MDS/Data/all_seurat_Post_integrated_sct_sub.rds')
saveRDS(data_ann, '/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')