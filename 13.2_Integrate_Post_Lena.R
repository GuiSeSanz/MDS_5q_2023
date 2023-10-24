

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





SAMPLE_NAMES <- c('FS-0406-post', 'FS-0634-post')


PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/Post/')
dir.create(PLOT_PATH, showWarnings = FALSE)



message('Reading Seurat samples: ')
all_seurat <- list()
colors <- ggthemes::tableau_color_pal('Classic 10')(length(SAMPLE_NAMES))
for (SAMPLE in SAMPLE_NAMES){
	message(SAMPLE)
	all_seurat[[SAMPLE]] <- readRDS(paste0(getwd(), '/Data/', SAMPLE, '/', SAMPLE, '_seurat_obj_norm.rds'))
}

merged_seurat <- merge(x = all_seurat[['FS-0406-post']], y = all_seurat[['FS-0634-post']], 
					add.cell.id = c('FS-0406-post', 'FS-0634-post'),  project='merged_Post')

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
saveRDS(anchors, 'Post_Integration_anchors.rds')
all_seurat_integrated_sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Post_integrated_sct.rds'))

message('Running clusterings and plotting')
all_seurat_integrated_sct <-  RunPCA(all_seurat_integrated_sct, verbose = FALSE)
PCA_dims <- how_many_PCs(all_seurat_integrated_sct, 'integrated_PCA_dim', FALSE)
all_seurat_integrated_sct <- FindNeighbors(object = all_seurat_integrated_sct, dims = 1:PCA_dims)
all_seurat_integrated_sct <- FindClusters(object = all_seurat_integrated_sct, resolution = c(0.2,0.4, 0.6, 0.8, 1.0, 1.2))
all_seurat_integrated_sct <- RunUMAP(all_seurat_integrated_sct, reduction = "pca", dims = 1:PCA_dims)

saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Post_integrated_sct.rds'))
# all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))

### WARNING the umap plots may differ, I din't save the original coords :c

pdf(paste0(PLOT_PATH, '/Umap_Post_integrated_samples.pdf'), width=12, height=10)
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

all_seurat_integrated_sct <- PrepSCTFindMarkers(object = all_seurat_integrated_sct)
saveRDS(all_seurat_integrated_sct, paste0(getwd(), '/Data/','all_seurat_Post_integrated_sct.rds'))
