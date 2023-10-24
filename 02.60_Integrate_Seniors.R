
library(Seurat)
library(argparse)
library(ggplot2)



how_many_PCs <- function(obj, pdf_name, print_pdf = TRUE) {
  # determine the correct PC number
  # Determine percent of variation associated with each PC
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  # last point where change of % of variation is more than 0.1%.
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  plot_df <- data.frame(pct = pct, 
           cumu = cumu, 
           rank = 1:length(pct))
  # Elbow plot to visualize 
  if(print_pdf){
	pdf(paste0(PLOT_PATH ,pdf_name,'.pdf'))
		print(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
		geom_text() + 
		geom_vline(xintercept = 90, color = "grey") + 
		geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
		theme_bw())
	dev.off()
	}
  return(pcs)
}

get_cells_per_cluster <- function(obj, assay){
  n_cells <- FetchData(obj, vars = c(paste0(assay, seq(0.2,1.2, by=0.2)))) 
  n_cells <- apply(n_cells, 2, table)
  results <- data.frame(Res=NULL, Cluster=NULL, Ncells=NULL)
  for (name in names(n_cells)){
    tmp <- setNames(as.data.frame(t(n_cells[[name]])), c('Res', 'Cluster', 'Ncells'))
    tmp$Res <- name
    results <- rbind(results, tmp)
  }
  return(results)
}


get_cell_distribution <- function(dataset, cols){
	p <- ggplot(dataset, aes(x=Cluster, y=Ncells, fill=Cluster)) + geom_bar(alpha=0.8, stat="identity") + 
	ggprism::theme_prism()  + theme(legend.position='none') + scale_fill_manual(values=cols)
	return(p)
}


plot_resolution_umap <- function(coords, resolution){
	resolution <- setNames(as.data.frame(all_seurat_integrated_sct[[resolution]]), 'Cluster')
	plotter <- merge(coords, resolution, by=0)
	p <- ggplot(plotter, aes(x=UMAP_1, y=UMAP_2, fill = Cluster)) + 
		geom_point(data = coords[sample(nrow(coords), 1e3),], fill='grey', color='grey', alpha=0.2, size=1) + 
		geom_point(size=1, alpha=1, pch=21, color='black') + theme_classic() + 
		theme(legend.position = "none", axis.ticks=element_blank(), axis.text=element_blank()) +
		scale_fill_manual(values = c(colors)) + facet_wrap(~Cluster) #+ guides(fill=guide_legend(nrow=2))
	return(p)
}


marker_list <- list()
marker_list[['hsc']]    <-c("CRHBP", "HOPX", "KYT", "CD34")
marker_list[['lmpp']]   <-c("PTPRC", "FLT3", "PROM1", "SATB1")
marker_list[['cc']]     <-c("CDC20", "TOP2A")
marker_list[['gmp']]    <-c("CSF3R", "CTSG", "PRTN3", "MPO")
marker_list[['granul']] <-c("ELANE", "AZU1", "CEBPA", "CEBPE", "CST7")
marker_list[['mono']]   <-c("LYZ", "CSTA")
marker_list[['dc']]     <-c("IRF8", "IRF7", "IL3RA", "CLEC4")
marker_list[['t']]      <-c("JCHAIN", "IKZF1", "CYTH1")
marker_list[['clp']]    <-c("IL7R", "DNTT")
marker_list[['prob']]   <-c("VPREB1", "EBF1", "CD79A", "CD79B")
marker_list[['mep']]    <-c("NFE2", "HFS1", "TAL1")
marker_list[['mk']]     <-c("PBX1", "MPL", "VWF", "FLI1", "ITGA22B", "GP1BA")
marker_list[['ery']]    <-c("GATA1",  "KLF1", "HBD", "HBB", "CA1", "AHSP")
marker_list[['baso']]   <-c("RUNX1", "HDC", "MS4A2", "MS4A3", "TPSAB1")


# SAMPLE_NAMES_5Q <- c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303')
#    SMD34459 SMD35109 SMD37209 SMD35303

# parser <- ArgumentParser(description='Process sc Data, no regression')
# parser$add_argument('-S', '--SampleNames', default="", type="character", nargs='+', action='append',
#                     help='Samples to integrate')
# args <- parser$parse_args()

# samples <- args$SampleNames
# samples <- samples[2]




SAMPLE_NAMES <- c('GSM5460411', 'GSM5460412', 'GSM5460413')


PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/Elder')
dir.create(PLOT_PATH, showWarnings = FALSE)



message('Reading Seurat samples: ')
all_seurat <- list()
colors <- ggthemes::tableau_color_pal('Classic 10')(length(SAMPLE_NAMES))
for (SAMPLE in SAMPLE_NAMES){
	message(SAMPLE)
	all_seurat[[SAMPLE]] <- readRDS(paste0(getwd(), '/Data/Normal_Data/', SAMPLE, '_seurat_obj_norm.rds'))
}

merged_seurat <- merge(x = all_seurat[['GSM5460411']], y = c(all_seurat[['GSM5460412']], all_seurat[['GSM5460413']]), 
					add.cell.id = c('GSM5460411', 'GSM5460412', 'GSM5460413'),  project='merged')

VariableFeatures(merged_seurat[["SCT"]]) <- rownames(merged_seurat[["SCT"]]@scale.data)
merged_seurat <-  RunPCA(merged_seurat, verbose = FALSE)
PCA_dims <- how_many_PCs(merged_seurat, 'merged_PCA_dim')
merged_seurat <- RunUMAP(merged_seurat, reduction = "pca", dims = 1:PCA_dims)


pdf(paste0(PLOT_PATH, 'Umap_merged_samples.pdf'), width=12, height=10)
	DimPlot(merged_seurat, reduction = "umap", group.by = "Sample", cols=colors)
dev.off()

message('Running Seurat integration') 
# seurat integration
features <- SelectIntegrationFeatures(object.list = all_seurat, nfeatures = 3000)
all_seurat <- PrepSCTIntegration(object.list = all_seurat, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = all_seurat, normalization.method = "SCT", anchor.features = features)
saveRDS(anchors, 'Elder_Integration_anchors.rds')
all_seurat_integrated_sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Elder_integrated_sct.rds'))

message('Running clusterings and plotting')
all_seurat_integrated_sct <-  RunPCA(all_seurat_integrated_sct, verbose = FALSE)
PCA_dims <- how_many_PCs(all_seurat_integrated_sct, 'integrated_PCA_dim', FALSE)
all_seurat_integrated_sct <- FindNeighbors(object = all_seurat_integrated_sct, dims = 1:PCA_dims)
all_seurat_integrated_sct <- FindClusters(object = all_seurat_integrated_sct, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))
all_seurat_integrated_sct <- RunUMAP(all_seurat_integrated_sct, reduction = "pca", dims = 1:PCA_dims)

saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Elder_integrated_sct.rds'))
# all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))

### WARNING the umap plots may differ, I din't save the original coords :c

pdf(paste0(PLOT_PATH, '/Umap_Elder_integrated_samples.pdf'), width=12, height=10)
	DimPlot(all_seurat_integrated_sct, reduction = "umap", group.by = "Sample", cols=colors)
dev.off()


library(clustree)
pdf(paste0(PLOT_PATH , '/Umap_Stats.pdf'), height = 8)
	print(clustree::clustree(all_seurat_integrated_sct@meta.data, prefix = "SCT_snn_res."))
	print(DimPlot(all_seurat_integrated_sct,  reduction = "umap", group.by = "orig.ident",  cols = c(ggthemes::tableau_color_pal('Classic 10 Medium')(4))) ) 
	features <-  c('nGene', 'nUMI', 'mitoRatio', 'RPSRatio', 'CD34')
	print(FeaturePlot(all_seurat_integrated_sct, features = features, order = TRUE) & NoLegend())
dev.off()

pdf(paste0(PLOT_PATH, '/allResUmaps.pdf'), height = 8)
populations <- get_cells_per_cluster(all_seurat_integrated_sct, 'integrated_snn_res.')
for (res in (paste0('integrated_snn_res.', seq(0.2,1.2, by=0.2)))){
	message(res)
	Idents(object = all_seurat_integrated_sct) <- res
	n_clusters <- length(populations[populations$Res == res,'Cluster'])
	if( n_clusters < 20){
		colors <- ggthemes::tableau_color_pal('Classic 20')(n_clusters)
	}else{
		colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
	}
	names(colors) <-  populations[populations$Res == res,'Cluster']
	print(
		cowplot::plot_grid(
		DimPlot(all_seurat_integrated_sct,
			reduction = "umap",
			label = TRUE,
			label.size = 6,
			cols= colors) + ggtitle(paste0('Resolution::',res)),
		get_cell_distribution(populations[populations$Res == res,], colors),
		nrow=2, rel_heights = c(3,1)
		)
	)
}
dev.off()

DefaultAssay(all_seurat_integrated_sct) <- "integrated"
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
for (res in (paste0('integrated_snn_res.', seq(0.2,1.2, by=0.2)))){
	message(res)
	pdf(paste0(PLOT_PATH, '/Umaps_per_cluster',res,'.pdf'))
		print(plot_resolution_umap(coords, res))
	dev.off()
}


library(future)
plan("multicore", workers = 64)


DefaultAssay(all_seurat_integrated_sct) <- "RNA"
all_markers <- list()
for (res in (paste0('integrated_snn_res.', 0.8))){
	message(res)
	Idents(object = all_seurat_integrated_sct) <- res
	markers <- FindAllMarkers(object = all_seurat_integrated_sct, 
								only.pos = TRUE,
								logfc.threshold = 0.25)
	all_markers[[res]] <- markers
	Markers_sep <-  split(markers, f=markers$cluster)
	WriteXLS::WriteXLS(Markers_sep, ExcelFileName=paste0(PLOT_PATH, 'Markers_', res, '.xlsx'), SheetNames = names(Markers_sep),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
}

saveRDS(all_markers, paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_markers.rds'))


# for (res in (paste0('integrated_snn_res.', seq(0.2,1.2, by=0.2)))){
#   message(res)
#   Idents(object = all_seurat_integrated_sct) <- res
#   Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster.1=NULL, cluster.2=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
#   for (cluster1 in  levels(all_seurat_integrated_sct[[res]][[1]])){
#     for (cluster2 in  levels(all_seurat_integrated_sct[[res]][[1]])){
#     	if(cluster1==cluster2){
#         	next
#     	}else{
# 			tmp <- FindMarkers(all_seurat_integrated_sct,
# 							ident.1 = cluster1,
# 							ident.2 = cluster2)
# 			tmp <- tmp[order(tmp$p_val_adj),]   
# 			tmp$gene_name <- rownames(tmp)
# 			tmp$cluster.1 <- as.character(cluster1)
# 			tmp$cluster.2 <- as.character(cluster2)
# 			tmp <- tmp[, c('gene_name', 'pct.1', 'pct.2', 'cluster.1', 'cluster.2',  'p_val', 'avg_logFC', 'p_val_adj')]
# 			tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_logFC) >2,]
# 		}
#       }
#       Results <- rbind(Results, tmp)
#     }
#   }
# 	Results_sp <- split(Results, Results$cluster.1)
#     WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'Markers_1Vs1_', res, '_parallel.xlsx'), SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
# }



get_markers <- function(cell_type){
	# cell_type <- 'common myeloid progenitor'
	all_markers <- read.table('./Data/Annotation/E-HCAD-6-marker-genes-files/E-HCAD-6.marker_genes_inferred_cell_type_-_ontology_labels.tsv', header=TRUE, sep='\t')
	all_markers$cluster <- toupper(all_markers$cluster)
	all_markers <- all_markers[all_markers$cluster == toupper(cell_type) & all_markers$pvals_adj < 0.05,]
	pdf(paste0('./Data/Annotation/',  gsub(" ", "_", cell_type), '_logfc_hist.pdf'))
		hist(all_markers$logfoldchanges, breaks=100, col='black', 
		main=paste0(cell_type, ' Nº of markers :: ', nrow(all_markers)), xlab='logFC', ylab='Frequency')
	dev.off()
	annotation <- setNames(read.table('./Data/Annotation/Homo_sapiens.GRCh38.Gene_name_ID.tsv', sep='\t'), c('genes', 'gene_name'))
	all_markers <- merge(all_markers, annotation,  by='genes')
	return(all_markers$gene_name)
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



get_violin_marker <- function(marker, clusters, data=norm_data){
	data = data[, marker, drop=FALSE]
	plotter <- setNames(merge(clusters, data, by=0), c('cell_id','cluster', 'value'))
	n_clusters <- length(unique(plotter$cluster))
	colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(n_clusters)
	plot <- ggplot(plotter, aes(x=cluster, y=value, fill=cluster)) + geom_violin(scale = "width") + theme_classic() +
	theme(legend.position='none') + scale_fill_manual(values=colors) + ggtitle(marker)
	return(plot)
}





# all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_Elder_integrated_sct.rds'))

# DefaultAssay(all_seurat_integrated_sct) <- 'integrated'
DefaultAssay(all_seurat_integrated_sct) <- 'RNA'

# options(future.globals.maxSize = 8000 * 1024^2)
# all_seurat_integrated_sct <- NormalizeData(all_seurat_integrated_sct)
# all_seurat_integrated_sct <- ScaleData(all_seurat_integrated_sct)
# norm_data <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$RNA@data)))
norm_data <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$SCT@data)))



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
# marker_list <- get_all_markers()

RESOLUTION <- '0.8'
res <- paste0('integrated_snn_res.', RESOLUTION)
clusters <- all_seurat_integrated_sct[[res]]

print(res)



# Only dotplots 
PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/Elder')

pdf(paste0(PLOT_PATH,'/Markers_per_Cluster/Markers_',RESOLUTION,'_SCT.pdf'))
for (marker in names(marker_list)){
	message(marker)
	features <- unique(marker_list[[marker]])
	print(DotPlot(all_seurat_integrated_sct, features = features, 
	group.by = res) + coord_flip() + ggtitle(marker) + 
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted('red'), mid ='white') +
	theme(legend.position='bottom', axis.text.x= element_text(angle = 45, vjust = 0.5, hjust=1),
			axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin()))
	
}
dev.off()

Idents(object = all_seurat_integrated_sct) <- res
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster.1=NULL, cluster.2=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster1 in  levels(all_seurat_integrated_sct[[res]][[1]])){
	for (cluster2 in  levels(all_seurat_integrated_sct[[res]][[1]])){
		if(cluster1==cluster2){
			next
		}else{
			message(paste0(cluster1, 'Vs', cluster2))
			tmp <- FindMarkers(all_seurat_integrated_sct,
							ident.1 = cluster1,
							ident.2 = cluster2)
			tmp <- tmp[order(tmp$p_val_adj),]
			tmp$gene_name <- rownames(tmp)
			tmp$cluster.1 <- as.character(cluster1)
			tmp$cluster.2 <- as.character(cluster2)
			tmp <- tmp[, c('gene_name', 'pct.1', 'pct.2', 'cluster.1', 'cluster.2',  'p_val', 'avg_log2FC', 'p_val_adj')]
			tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_log2FC) >0.2,]
		}
		Results <- rbind(Results, tmp)
	}
}

Results_sp <- split(Results, Results$cluster.1)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'Markers_1Vs1_', res, '_parallel.xlsx'), SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


for (marker in names(marker_list)){
	message(marker)
	features <- unique(marker_list[[marker]])
	violins <- list()
	for (feature in features){
		if(feature %in% colnames(norm_data)){
			violins[[feature]] <- get_violin_marker(feature, clusters)
		}
	}
	pdf(paste0(PLOT_PATH,'Markers_per_Cluster/Markers_',marker,'_',RESOLUTION,'_SCT.pdf'))

	print(DotPlot(all_seurat_integrated_sct, features = features, 
	group.by = res) + coord_flip() + 
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted('red'), mid ='white') +
	theme(legend.position='bottom', axis.text.x= element_text(angle = 45, vjust = 0.5, hjust=1),
			axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin()))
	
	full_pages <- length(violins)%/%2 +1
	for (page in seq(0, full_pages, by=2)){
		plot1 <- names(violins)[page+1]
		plot2 <- names(violins)[page+2]
		print(paste0(plot1, '---' ,plot2))
		print(cowplot::plot_grid(
			violins[[plot1]], violins[[plot2]],
		nrow=2))
	}
	if (length(violins)%%2 != 0){
		plot3 <- names(violins)[length(violins)]
		print(plot3)
		print(cowplot::plot_grid(
			violins[[plot3]], NULL,
		nrow=2))
	}
	
	dev.off()
}




# SubCluster Time!!

all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_Elder_integrated_sct.rds'))
DefaultAssay(all_seurat_integrated_sct) <- 'RNA'

# options(future.globals.maxSize = 8000 * 1024^2)
norm_data <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$SCT@data)))
Idents(all_seurat_integrated_sct) <- 'integrated_snn_res.0.8'



clusters_2_subcluster <- c('5', '6', '8', '12', '13', '14', '18', '20', '21')

for (cluster in clusters_2_subcluster){
	print(cluster)
	all_seurat_integrated_sct <- FindSubCluster(all_seurat_integrated_sct, cluster, 'integrated_snn', subcluster.name = paste0("integrated_snn_res.0.8_sub"), resolution = 0.1)
	Idents(all_seurat_integrated_sct) <- "integrated_snn_res.0.8_sub"
}

saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_subcluster.rds'))
# dotplots

PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/Elder')



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



pdf(paste0(PLOT_PATH,'/Markers_per_Cluster/Markers_0.8_SCT_Subcluster.pdf'))
for (marker in names(marker_list)){
	message(marker)
	features <- unique(marker_list[[marker]])
	print(DotPlot(all_seurat_integrated_sct, features = features, 
	group.by =  "integrated_snn_res.0.8_sub") + coord_flip() + ggtitle(marker) + 
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted('red'), mid ='white') +
	theme(legend.position='bottom', axis.text.x= element_text(angle = 90, vjust = 0.5, hjust=1),
			axis.title.x=element_blank(),legend.box="vertical", legend.margin=margin()))
	
}
dev.off()



# markers 1VsAll


DefaultAssay([all_seurat_integrated_sct]) <- "RNA"
all_markers <- list()
for (res in (paste0('integrated_snn_res.', '0.8_sub'))){
	message(res)
	Idents(object = all_seurat_integrated_sct) <- res
	markers <- FindAllMarkers(object = all_seurat_integrated_sct, 
								only.pos = TRUE,
								logfc.threshold = 0.25)
	all_markers[[res]] <- markers
	Markers_sep <-  split(markers, f=markers$cluster)
	WriteXLS::WriteXLS(Markers_sep, ExcelFileName=paste0(PLOT_PATH, 'Markers_', res, '_Subcluster.xlsx'), SheetNames = names(Markers_sep),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
}

saveRDS(all_markers, paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_markers_subcluster.rds'))


# markers 1Vs1

Idents(object = all_seurat_integrated_sct) <- "integrated_snn_res.0.8_sub"
res <- 'integrated_snn_res.0.8_sub'
Results <- data.frame(gene_name=NULL, pct.1=NULL, pct.2=NULL, cluster.1=NULL, cluster.2=NULL, p_val=NULL, avg_logFC=NULL, p_val_adj=NULL )
for (cluster1 in  unique(all_seurat_integrated_sct[[res]][[1]])[stringr::str_which(unique(all_seurat_integrated_sct[[res]][[1]]), '_')]){
	for (cluster2 in  unique(all_seurat_integrated_sct[[res]][[1]])){
		if(cluster1==cluster2){
			next
		}else{
			message(paste0(cluster1, 'Vs', cluster2))
			tmp <- FindMarkers(all_seurat_integrated_sct,
							ident.1 = cluster1,
							ident.2 = cluster2)
			tmp <- tmp[order(tmp$p_val_adj),]
			tmp$gene_name <- rownames(tmp)
			tmp$cluster.1 <- as.character(cluster1)
			tmp$cluster.2 <- as.character(cluster2)
			tmp <- tmp[, c('gene_name', 'pct.1', 'pct.2', 'cluster.1', 'cluster.2',  'p_val', 'avg_log2FC', 'p_val_adj')]
			tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_log2FC) >0.2,]
		}
		Results <- rbind(Results, tmp)
	}
}

Results_sp <- split(Results, Results$cluster.1)
WriteXLS::WriteXLS(Results_sp, ExcelFileName=paste0(PLOT_PATH,'Markers_1Vs1_', res, '_Subcluster.xlsx'), SheetNames = names(Results_sp),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


DefaultAssay(all_seurat_integrated_sct) <- "integrated"
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)

coords <- merge(coords, all_seurat_integrated_sct$integrated_snn_res.0.8_sub, by=0)
coords <- merge(coords, all_seurat_integrated_sct$RPSRatio, by.x='Row.names', by.y=0)
coords <- setNames(coords, c('Row.names', 'UMAP1', 'UMAP2',  'Cluster', 'RiboRatio'))

pdf('./Plots/tmp.pdf')
ggplot()+ 
geom_point(data= coords[!coords$Cluster %in% c("6_1", "6_0"), ],aes(x=UMAP1, y=UMAP2), color='grey', alpha=0.3) +
geom_point(data= coords[ coords$Cluster %in% c("6_1", "6_0"), ],aes(x=UMAP1, y=UMAP2, color=RiboRatio), alpha=0.8) + theme_classic() + scale_color_viridis_c(option='turbo')
ggplot()+ 
geom_point(data= coords[!coords$Cluster %in% c("6_1", "6_0"), ],aes(x=UMAP1, y=UMAP2), color='grey', alpha=0.3) + 
geom_point(data= coords[ coords$Cluster %in% c("6_1", "6_0"), ],aes(x=UMAP1, y=UMAP2, color=Cluster), alpha=0.8) + scale_color_manual(values=ccolss[c(1,4)]) + theme_classic()
dev.off()


ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f",
"#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f",
"#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")

colors <-  c((ggthemes::tableau_color_pal('Classic 20')(20)), ccolss)
coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
res <- "integrated_snn_res.0.8_sub"
	pdf(paste0(PLOT_PATH, '/Umaps_per_cluster',res,'.pdf'))
		print(plot_resolution_umap(coords, res))
	dev.off()


