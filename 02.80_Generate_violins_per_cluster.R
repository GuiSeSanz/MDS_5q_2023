library(Seurat)
library(argparse)
library(ggplot2)
library(future)
plan("multicore", workers = 64)


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

parser <- ArgumentParser(description='Get the DE genes between clusters, 1Vs1 on the resolution selected')
parser$add_argument('-R', '--Resolution', default="0.2", type="double",
                    help='Resolution of the clustering')
args <- parser$parse_args()
RESOLUTION <- args$Resolution


PLOT_PATH <- paste0(getwd(), '/Plots/Integrated/')
dir.create(PLOT_PATH, showWarnings = FALSE)
all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/all_seurat_integrated_sct.rds'))

DefaultAssay(all_seurat_integrated_sct) <- 'integrated'
# DefaultAssay(all_seurat_integrated_sct) <- 'RNA'

# options(future.globals.maxSize = 8000 * 1024^2)
# all_seurat_integrated_sct <- NormalizeData(all_seurat_integrated_sct)
# all_seurat_integrated_sct <- ScaleData(all_seurat_integrated_sct)
# norm_data <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$RNA@data)))
norm_data <- as.data.frame(t(as.data.frame(all_seurat_integrated_sct@assays$SCT@data)))

res <- paste0('integrated_snn_res.', RESOLUTION)
clusters <- all_seurat_integrated_sct[[res]]

print(res)




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
marker_list <- get_all_markers()


for (marker in names(marker_list)){
	message(marker)
	features <- unique(marker_list[[marker]])
	violins <- list()
	for (feature in features){
		if(feature %in% colnames(norm_data)){
			violins[[feature]] <- get_violin_marker(feature, clusters)
		}
	}
	pdf(paste0('./Plots/Integrated/Markers_per_Cluster/Markers_',marker,'_',RESOLUTION,'_SCT.pdf'))

	print(DotPlot(all_seurat_integrated_sct, features = features, 
	group.by = res) + coord_flip() + 
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	viridis::scale_colour_viridis(option="plasma") + 
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




marker_list <- get_all_markers()

for (marker in names(marker_list)){
	message(marker)
	features <- unique(marker_list[[marker]])
	violins <- list()
	for (feature in features){
		if(feature %in% colnames(norm_data)){
			violins[[feature]] <- get_violin_marker(feature, clusters)
		}
	}
	pdf(paste0('./Plots/Integrated/Markers_per_Cluster/Custom_Markers_',marker,'_',RESOLUTION,'_SCT.pdf'))

	print(DotPlot(all_seurat_integrated_sct, features = features, 
	group.by = res) + coord_flip() + 
	geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
	viridis::scale_colour_viridis(option="plasma") + 
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