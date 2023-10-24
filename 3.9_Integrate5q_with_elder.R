
library(Seurat)
library(argparse)
library(ggplot2)



how_many_PCs <- function(obj){
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
  return(pcs)
}

get_umap_signature <- function(df, title){
	plot <- ggplot(df, aes(x=UMAP_1, y=UMAP_2, color= Signature)) + 
	geom_point(alpha=0.9, size = 0.8) + 
	# scale_color_gradientn(low="grey90", high ="blue", name = 'Expression') + 
	# guides(color = guide_colourbar(barwidth = 0.5, barheight = 2, label = TRUE)) + 
	ggprism::theme_prism() + ggtitle(title)+  labs(x = "UMAP 1", y = 'UMAP 2') +
	theme(legend.position='none', plot.title = element_text(hjust = 0.5, size = 12),  
	axis.text=element_blank(), axis.ticks=element_blank())#, legend.key.height = unit(2, 'mm'), legend.key.width = unit(1, 'mm')) + 
	if(is.numeric(df$Signature)){
		  plot <- plot + viridis::scale_color_viridis()
	}else{
	  plot <- plot + scale_color_manual(values=c('High'='#30A3CC', 'Low'='#d3d3d3'))
	}
	return(plot)
}

get_expression_signature <- function(gene_list, df, coords = coords, upplim=NULL, split_HL=FALSE, noPlot=FALSE){
		check_genes <- function(gene_list, df){
			not_found_genes <- gene_list[which(!gene_list %in% colnames(df))]
			if(length(not_found_genes)!=0){
				message(paste0('Left behind ', length(not_found_genes), ': ', paste0(not_found_genes, collapse='; ')))
				gene_list <- gene_list[which(gene_list %in% colnames(df))]
			}else{
				gene_list <- gene_list
			}

			return(gene_list)
		}
	gene_list <- check_genes(gene_list, df)
	tmp <- df[, colnames(df) %in% gene_list]
	tmp$Signature <- rowMeans(tmp)
	# tmp$cell_id <- sub('-', '\\.',rownames(tmp))
	tmp <- merge(tmp, coords, by=0)
	if(!is.null(upplim)){
		upp_outliers <- which(tmp$Signature>upplim) 
		tmp[upp_outliers , 'Signature'] <- upplim
	}
	plot <- get_umap_signature(tmp, 'title')
	if(split_HL){
		plot <- plot + facet_wrap(~BinScore)
	}
	if(noPlot){
		return(tmp)
	}else{
		return(plot)
	}
}



parser <- ArgumentParser(description='Perform the DE analysis with the normals')
parser$add_argument('-T', '--Type', default="ANCHOR", type="character",
					help='ANCHOR or 5Q')
args <- parser$parse_args()
TYPE <- args$Type


data_path <- '/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/RawData/'
elder_list <- c('GSM5460411_elderly1_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460412_elderly2_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460413_elderly3_filtered_feature_bc_matrix_counts.tsv')


all_normal_data <- list()
for (elder_file in elder_list){
	elder <- stringr::str_extract(elder_file,'(?<=_)[a-z]+[\\d]{1}')
	SAMPLE <- stringr::str_extract(elder_file, '[A-Z0-9]+(?=_young|_elderly)')
	tmp <- readRDS(paste0('/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/', SAMPLE, '_seurat_obj_norm.rds'))
	all_normal_data[[SAMPLE]] <- tmp
}

# merge all the samples on a single object
all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))
DefaultAssay(all_seurat_integrated_sct) <- "SCT"
all_seurat_integrated_sct[['integrated']] <- NULL
mds_5q_samples <- SplitObject(all_seurat_integrated_sct, split.by='Sample')


all_samples <- mds_5q_samples
for (sample in names(all_normal_data)){
	all_samples[[sample]] <- all_normal_data[[sample]]
}
# integrate them all

# features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)ç

genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
genes_expanded <- as.character(unique(genes_expanded$V5))

comon_genes <- c()
for (sample in names(all_samples)){
	tmp <- rownames(all_samples[[sample]])
	if (sample == names(all_samples)[1]){
		comon_genes <- tmp
	}else{
		comon_genes <- intersect(comon_genes, tmp)
	}
}

common_genes_5q <- genes_expanded[which(genes_expanded %in% comon_genes)]

# delete <- c('CD180', 'SERF1B', 'ZNF366', 'SV2C', 'BHMT2', 'ADGRV1', 'LUCAT1', 
#'LINC00992', 'PRR16', 'ACSL6', 'IL4', 'GDF9', 'TIFAB', 'SLC23A1', 'PCDHB1', 
#'PCDHGA7', 'PLAC8L1', 'PPP2R2B', 'ADRB2', 'ITK', 'PCDHA10', 'NUDT12', 'CATSPER3')

all_samples <- lapply(X = all_samples, FUN = function(x){
	 x <- SCTransform(x, variable.features.n = nrow(x))
})

### geneset + anchor
if (TYPE == "ANCHOR"){
	delete <- c('CD180','SERF1B','ZNF366','SV2C','BHMT2','ADGRV1','LUCAT1','LINC00992',
	'PRR16','ACSL6','IL4','GDF9','TIFAB','SLC23A1','PCDHB1','PCDHGA7','PLAC8L1','PPP2R2B',
	'ADRB2','ITK', 'PCDHA10', 'NUDT12', 'CATSPER3')
	features <- SelectIntegrationFeatures(object.list = all_samples, nfeatures = 2000)
	features <- unique(c(features, common_genes_5q))
	features <- features[!(features %in% delete)]
	all_samples <- PrepSCTIntegration(object.list = all_samples,  anchor.features = features)
	anchors      <- FindIntegrationAnchors(object.list = all_samples,
													normalization.method = "SCT", 
													anchor.features = features)
	combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
	combined.sct <- RunPCA(combined.sct, verbose = FALSE)
	combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)

	saveRDS(combined.sct, '/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_AnchorsAndChr5_sct.rds')
	}




#### Only geneset
if (TYPE == "5Q"){
	features <- common_genes_5q
	delete <- c('CD180','SERF1B','ZNF366','SV2C','BHMT2','ADGRV1','LUCAT1','LINC00992','PRR16','ACSL6','IL4',
	'GDF9','TIFAB','SLC23A1','PCDHB1','PCDHGA7','PLAC8L1','PPP2R2B','ADRB2','ITK', 'PCDHA10', 'NUDT12', 'CATSPER3')
	features <- features[!(features %in% delete)]
	all_samples         <- PrepSCTIntegration(object.list = all_samples, 
											anchor.features = features)
	anchors      <- FindIntegrationAnchors(object.list = all_samples,
													normalization.method = "SCT", 
													anchor.features = features)
	combined.sct_geneset  <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
	combined.sct_geneset  <- RunPCA(combined.sct_geneset, verbose = FALSE)
	PCAs <- how_many_PCs(combined.sct_geneset)
	combined.sct_geneset  <- RunTSNE(combined.sct_geneset, reduction = "pca", dims = 1:PCAs)
	combined.sct_geneset  <- RunUMAP(combined.sct_geneset, reduction = "pca", dims = 1:PCAs)
	saveRDS(combined.sct_geneset, '/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
}


elder_list <- c('GSM5460411_elderly1_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460412_elderly2_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460413_elderly3_filtered_feature_bc_matrix_counts.tsv')


genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
genes_expanded <- as.character(unique(genes_expanded$V5))


elder_names <- stringr::str_extract(elder_list, '^[A-Z0-9]+')
combined.sct_geneset <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
combined.sct_geneset_markers <- FindAllMarkers(combined.sct_geneset, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)


coords <- as.data.frame(combined.sct_geneset@reductions$umap@cell.embeddings)
coords$Sample <- stringr::str_extract(rownames(coords), "(^[A-Z0-9]+)")


coords$Control <- ifelse(coords$Sample %in% elder_names, "Elder", "5q")

pdf('./Plots/Integrated/Integration_5q_elder_Only5q.pdf')
ggplot(coords, aes(x=UMAP_1, y = UMAP_2, col=Sample)) + 
geom_point(alpha=0.8, size = 0.2) + theme_classic() + ggthemes::scale_colour_tableau()
ggplot(coords, aes(x=UMAP_1, y = UMAP_2, col=Control)) + 
geom_point(alpha=0.8, size = 0.2) + theme_classic() + ggthemes::scale_colour_tableau()
dev.off()


test_comp <- as.data.frame(table(combined.sct_geneset$seurat_clusters, combined.sct_geneset$Sample))
combined.sct_geneset <- FindNeighbors(combined.sct_geneset, dims = 1:PCAs)
combined.sct_geneset <- FindClusters(combined.sct_geneset)
combined.sct_geneset$Control <-  ifelse(combined.sct_geneset$Sample %in% elder_names, "Elder", "5q")
composition <-  as.data.frame(table(combined.sct_geneset$seurat_clusters, combined.sct_geneset$Control))
composition$Var1 <- as.numeric(composition$Var1)
test_comp$Var1 <- as.numeric(test_comp$Var1)
pdf('./Plots/Integrated/Integrated_elder_anchors_Cluster_comp_5q.pdf')
	ggplot(composition, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat='identity', position='fill') + 
	theme_classic() + ggthemes::scale_fill_tableau()
	ggplot(test_comp, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat='identity', position='fill') + 
	theme_classic() + ggthemes::scale_fill_tableau() + facet_wrap(~Var2)
dev.off()





coords <- as.data.frame(combined.sct_geneset@reductions$umap@cell.embeddings)
coords$Sample <- stringr::str_extract(rownames(coords), "(^[A-Z0-9]+)")
coords$Control <- ifelse(coords$Sample %in% elder_names, "Elder", "5q")

normData <- t(as.data.frame(combined.sct_geneset@assays$SCT@data))
tmp <- setNames(as.data.frame(rowMeans(normData[,colnames(normData) %in% genes_expanded])), c('Signature'))
coords <- merge(coords, tmp, by=0)
coords <- merge( coords, setNames(as.data.frame(combined.sct_geneset$seurat_clusters), 'Cluster'), by.x= 'Row.names', by.y=0)

cmp <- as.data.frame(table(combined.sct_geneset$seurat_clusters))
cmp$Var1 <- as.numeric(cmp$Var1)

cmp2 <- as.data.frame(table(combined.sct_geneset$seurat_clusters, combined.sct_geneset$Sample))

pdf('./Plots/Integrated/Integrated_elder_Only_5q_cluster_comp_sample.pdf')
	ggplot(cmp2, aes(x=Var1, y=Freq, fill=Var2)) + 
	geom_bar(stat='identity', position='fill') + 
	theme_classic() + ggthemes::scale_fill_tableau() 
dev.off()

clusters_2_plot <- c(28,29, 30, 31, 32, 33, 34)
pdf('./Plots/Integrated/Integrated_elder_Only_5q_Signature_hist.pdf')
	ggplot(coords, aes(x = UMAP_1, y = UMAP_2, color= Signature)) + geom_point(size=0.2) + 
	theme_classic() + viridis::scale_color_viridis() + facet_wrap(~Control)
	ggplot(coords, aes(x= Signature, fill=Sample)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.8)
	ggplot(coords, aes(x= Signature, fill=Cluster)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.8)
dev.off()

a <- aggregate(coords[, 'Signature'], list(coords$Cluster), mean)
colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(35)
pdf('./Plots/Integrated/Integrated_elder_Only_5q_Signature_per_Cluster.pdf')
ggplot(coords , aes(x=Cluster, y=Signature, fill=Cluster)) + geom_boxplot() + 
	theme_classic() + scale_fill_manual(values=colors) + theme(legend.position='none')
dev.off()


##### Plot the nUMIS 
tmp <- setNames(as.data.frame(combined.sct_geneset$nUMI), 'nUMI_all_genes')
coords <- merge(coords, tmp, by.x='Row.names', by.y=0)
tmp <- as.data.frame(combined.sct_geneset@assays$RNA@counts)
tmp_5q <- setNames(as.data.frame(colSums(tmp[rownames(tmp) %in% genes_expanded,])), 	'nUMI_5q_genes')
tmp_no_5q <- setNames(as.data.frame(colSums(tmp[!rownames(tmp) %in% genes_expanded,])), 	'nUMI_no_5q_genes')
tmp <- merge(tmp_no_5q, tmp_5q, by=0)
coords <- merge(coords, tmp, by.x='Row.names')

# a <- aggregate(coords[, 'Signature'], list(coords$Cluster), mean)
colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(35)
pdf('./Plots/Integrated/Integrated_elder_Only_5q_UMIS_per_Cluster.pdf')
ggplot(coords , aes(x=Cluster, y=nUMI_all_genes, fill=Cluster)) + geom_boxplot() + 
	theme_classic() + scale_fill_manual(values=colors) + theme(legend.position='none')
ggplot(coords , aes(x=Cluster, y=nUMI_no_5q_genes, fill=Cluster)) + geom_boxplot() + 
	theme_classic() + scale_fill_manual(values=colors) + theme(legend.position='none')
ggplot(coords , aes(x=Cluster, y=nUMI_5q_genes, fill=Cluster)) + geom_boxplot() + 
	theme_classic() + scale_fill_manual(values=colors) + theme(legend.position='none')
dev.off()





genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
genes_expanded <- as.character(unique(genes_expanded$V5))

normData <- as.data.frame(combined.sct_geneset@assays$SCT@data)
normData_t <- t(normData)

tmp <- setNames(as.data.frame(rowMeans(normData_t[,colnames(normData_t) %in% genes_expanded])), c('Signature'))

coords <- merge(coords, tmp, by=0)
pdf('./Plots/Integrated/Integrated_elder_5q_Signature_hist.pdf')
	ggplot(coords, aes(x= Signature, fill=Control)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.6) + facet_wrap(~Sample)
	ggplot(coords, aes(x= Signature, fill=Control)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.6)
	ggplot(coords, aes(x= Signature, fill=Sample)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.6)
dev.off()

pdf('./Plots/Integrated/Integrated_elder_5q_Signature_UMAP.pdf')
	ggplot(coords, aes(x = UMAP_1, y = UMAP_2, color= Signature)) + geom_point() + 
	theme_classic() + viridis::scale_color_viridis() + facet_wrap(~Control)
dev.off()

pdf('./Plots/Integrated/Integrated_elder_5q_clusters.pdf')
DimPlot(combined.sct_geneset, reduction = "umap")
dev.off()

combined.sct_geneset_markers <- FindAllMarkers(combined.sct_geneset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

RPS14
CD74







################################################################################

combined.sct <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_AnchorsAndChr5_sct.rds')
PCAs <- how_many_PCs(combined.sct)
combined.sct <- FindNeighbors(combined.sct, dims = 1:PCAs)
combined.sct <-  FindClusters(combined.sct, resolution = 0.5)
combined.sct_markers <- FindAllMarkers(combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


features_2_plot <- c()
for( cluster in unique(combined.sct_markers$cluster)){
	tmp <- combined.sct_markers[combined.sct_markers$cluster == cluster,]
	tmp  <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
	features_2_plot <- c(features_2_plot, tmp[1:10, 'gene'])
}



pdf('./Plots/Integrated/Integrated_elder_anchors.pdf')
DefaultAssay(combined.sct) <- "integrated"
DimPlot(combined.sct, reduction = "umap")
DefaultAssay(combined.sct) <- "SCT"
DoHeatmap(combined.sct, features = features_2_plot) + NoLegend()
dev.off()

combined.sct$Control <-  ifelse(combined.sct$Sample %in% elder_names, "Elder", "5q")
pdf('./Plots/Integrated/Integrated_elder_anchors_check.pdf')
DimPlot(combined.sct, reduction = "umap", group.by = "Control")
dev.off()

features_2_plot <- rownames(combined.sct)[rownames(combined.sct) %in% genes_expanded]
pdf('./Plots/Integrated/Integrated_elder_anchors_2.pdf')
DoHeatmap(combined.sct, features = features_2_plot) + NoLegend()
dev.off()


features_2_plot <- c('ZNF608','RPS23','KLHL3','PIK3R1','SPARC','ANXA6','S100Z',
'NREP','ADAM19','DTWD2','HBEGF','SSBP2','RPS14','CAST','CYFIP2','MEF2C','SMIM3',
'DHFR','GLRX','P4HA2','EGR1','IRF1','CRHBP','CCNB1','TGFBI','MZB1','CD74','EBF1')
d1 <- split(features_2_plot, ceiling(seq_along(features_2_plot)/4))

for (feature in features_2_plot){
pdf(paste0('./Plots/Integrated/Integrated_elder_anchors_violins_', feature,'.pdf'))
	print(VlnPlot(combined.sct, features = feature,pt.size = FALSE))
dev.off()
}


composition <-  as.data.frame(table(combined.sct$seurat_clusters, combined.sct$Control))
pdf('./Plots/Integrated/Integrated_elder_anchors_Cluster_comp.pdf')
	ggplot(composition, aes(x=Var1, y=Freq, fill=Var2)) + geom_bar(stat='identity', position='fill') + 
	theme_classic() + scale_fill_viridis_d(alpha=0.9)
dev.off()


coords <- as.data.frame(combined.sct@reductions$umap@cell.embeddings)
coords$Sample <- stringr::str_extract(rownames(coords), "(^[A-Z0-9]+)")
coords$Control <- ifelse(coords$Sample %in% elder_names, "Elder", "5q")

normData <- t(as.data.frame(combined.sct@assays$SCT@data))

tmp <- setNames(as.data.frame(rowMeans(normData[,colnames(normData) %in% genes_expanded])), c('Signature'))
coords <- merge(coords, tmp, by=0)
pdf('./Plots/Integrated/Integrated_elder_5q_Signature_hist.pdf')
	ggplot(coords, aes(x= Signature, fill=Control)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.6) + facet_wrap(~Sample)
	ggplot(coords, aes(x= Signature, fill=Control)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.6)
	ggplot(coords, aes(x= Signature, fill=Sample)) + geom_histogram(bins=200) + 
	theme_classic() + scale_fill_viridis_d(alpha=0.8)
dev.off()

pdf('./Plots/Integrated/Integrated_elder_5q_Signature_UMAP.pdf')
	ggplot(coords, aes(x = UMAP_1, y = UMAP_2, color= Signature)) + geom_point(size=0.2) + 
	theme_classic() + viridis::scale_color_viridis() + facet_wrap(~Control)
dev.off()



