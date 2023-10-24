
library(Seurat)
library(ggplot2)
options(stringsAsFactors=FALSE)



# less restrictive

get_genes <- function(sample_name_filter){
	# sample_name_filter <- elder or non5q
	genes <- data.frame(X=NULL, log2FoldChange=NULL)
	for (fl in de_files){
		if (stringr::str_detect(fl, stringr::fixed(sample_name_filter, ignore_case=TRUE))){
			tmp <- read.table(paste(data_folder, fl, sep = '/'), sep=',', header=TRUE)
			tmp <- tmp[!is.na(tmp$padj),]
			tmp <- tmp[tmp$padj<0.05, c('X', 'log2FoldChange')]
			genes <- rbind(genes, tmp)
		}
	}
	print(dim(genes))
	genes <- setNames(genes, c('gene_name', 'log2FoldChange'))
	return(genes)
}

HL_signature <- function(data, signature_df){
	data <- data[rownames(data) %in% signature_df$gene_name,]
    logpluspointfive <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(data, 2, logpluspointfive))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))
	zscore$gene_name <- rownames(zscore)
    ann_markers <- merge(zscore, signature_df, by='gene_name')
	ann_markers <- ann_markers[ !duplicated(ann_markers$gene_name),]
	res <- ann_markers[, !colnames(ann_markers) %in% c('gene_name', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
	rownames(res) <- ann_markers$gene_name
	res <- setNames(as.data.frame(colSums(res[, !colnames(res) %in% c('gene_name', 'log2FoldChange', 'padj'), drop=FALSE], )), c('High_pondered'))

} 

get_histogram <- function(data, mode='CASE'){
	# plotter <- merge(data, coords, by=0)
	plotter <- data
	plotter$Sample <- stringr::str_extract(plotter$Row.names, '^[A-Z0-9]+')
	plotter$Case <- ifelse( stringr::str_detect(plotter$Sample, '^GSM'), 'Elder', '5q')
	if (mode == 'SAMPLE'){
		p <- ggplot(plotter, aes(x= High_pondered, fill=Sample)) + geom_histogram(bins=200, alpha = 0.7, position = "identity") + 
			theme_classic() + scale_fill_viridis_d(option='turbo')
	} else {
	p <- ggplot(plotter, aes(x= High_pondered, fill=Case)) + geom_histogram(bins=200, alpha = 0.7, position = "identity") + 
		theme_classic() + scale_fill_viridis_d(option='turbo')
	}
	return(p)
}

get_histogram_all <- function(data, mode='CASE'){
	# plotter <- merge(data, coords, by=0)
	plotter <- data
	plotter$Sample <- stringr::str_extract(plotter$Row.names, '^[A-Z0-9]+')
	plotter$Case <- ifelse( plotter$Sample %in%  c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303'), '5q',
					ifelse( plotter$Sample %in%  c('GSM5460411','GSM5460412','GSM5460413'), 'Elder', 'MDS_Non_5q'))
	if (mode == 'SAMPLE'){
		p <- ggplot(plotter, aes(x= High_pondered, fill=Sample)) + geom_histogram(bins=200, alpha = 0.7, position = "identity") + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)
	} else {
	p <- ggplot(plotter, aes(x= High_pondered, fill=Case)) + geom_histogram(bins=200, alpha = 0.7, position = "identity") + 
		theme_classic() + scale_fill_viridis_d(option='turbo',alpha=0.7)
	}
	return(p)
}


get_ridges <- function(data, color_by = 'Sample'){
	# plotter <- merge(data, coords, by=0)
	plotter <- data
	plotter$Sample <- stringr::str_extract(plotter$Row.names, '^[A-Z0-9]+')
	plotter$Case <- ifelse( plotter$Sample %in%  c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303'), '5q',
					ifelse( plotter$Sample %in%  c('GSM5460411','GSM5460412','GSM5460413'), 'Elder', 'MDS_Non_5q'))
	if ( color_by == 'Case'){
		p<-  ggplot(plotter, aes(x= High_pondered, y = Sample, fill=Case)) + ggridges::geom_density_ridges2() + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)
	}else{
		p <- ggplot(plotter, aes(x= High_pondered, y = Sample, fill=Sample)) + ggridges::geom_density_ridges2() + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)

	}
	return(p)
}

get_umap <- function(data){
	plotter <- merge(data, coords, by=0)
	plotter$Sample <- stringr::str_extract(plotter$Row.names, '^[A-Z0-9]+')
	plotter$Case <- ifelse( stringr::str_detect(plotter$Sample, '^GSM'), 'Elder', '5q')
	ggplot(plotter, aes(x= UMAP_1, y= UMAP_2, color=High_pondered)) + geom_point(alpha = 0.8) + 
	theme_classic() + scale_colour_gradientn(colours = shades::saturation(viridis::turbo(250), shades::scalefac(0.6)))
}

get_boxplot <- function(data){
	colors <- colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(35)
	plotter <- merge(data, coords, by=0)
	plotter <- merge( plotter, setNames(as.data.frame(combined.sct_geneset$seurat_clusters), 'Cluster'), by.x= 'Row.names', by.y=0)
	plotter$Sample <- stringr::str_extract(plotter$Row.names, '^[A-Z0-9]+')
	plotter$Case <- ifelse( stringr::str_detect(plotter$Sample, '^GSM'), 'Elder', '5q')
	p <- ggplot(plotter , aes(x=Cluster, y=High_pondered, fill=Cluster)) + geom_boxplot() + 
		theme_classic() + scale_fill_manual(values=colors) + theme(legend.position='none')
	return(p)
}

combined.sct_geneset <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
coords <- as.data.frame(combined.sct_geneset@reductions$umap@cell.embeddings)

de_files <- c(  'HSC5qvsElderly.tsv',
				'HSCs5qvsnon5q.tsv',
				'MEPs5qvsElderly.csv',
				'MEPs5qvsnon5q.tsv',
				'CMPs5qvsElderly.tsv',
				'CMPs5qvsnon5q.tsv', 
				'GMPs5qvsElderly.tsv',
				'GMPs5qvsnon5q.tsv')

data_folder <- '/home/sevastopol/data/gserranos/MDS/Data/Annotation'


data <- as.data.frame(combined.sct_geneset@assays$SCT@data)
plotter_non_5q <- HL_signature(data, get_genes('non5q'))
plotter_elder  <- HL_signature(data, get_genes('elder'))


pdf('./Plots/Integrated/Integrated_HL_Signature_elder.pdf')
get_histogram(plotter_elder)
get_histogram(plotter_elder, 'SAMPLE')
get_umap(plotter_elder) + facet_wrap(~Case)
get_boxplot(plotter_elder)
dev.off()


pdf('./Plots/Integrated/Integrated_HL_Signature_non5q.pdf')
get_histogram(plotter_non_5q)
get_histogram(plotter_non_5q, 'SAMPLE')
get_umap(plotter_non_5q) + facet_wrap(~Case)
get_boxplot(plotter_non_5q)
dev.off()

pdf('./Plots/Integrated/Integrated_HL_Signature_non5q_Lite.pdf')
get_histogram(plotter_non_5q)
get_histogram(plotter_non_5q, 'SAMPLE')
get_histogram(plotter_non_5q, 'SAMPLE') + facet_wrap(~Sample, nrow=2) + 
theme(legend.position='none')
get_ridges(plotter_non_5q)
dev.off()



genes <- get_genes('non5q')

all_samples <-  as.data.frame(combined.sct_geneset@assays$SCT@data)
# all_samples_bck <- all_samples
# all_samples <- all_samples_bck
all_samples <- all_samples[rownames(all_samples) %in% genes$gene_name,]
normal_mds <- list.files('./Data/', pattern= '^SMD')
normal_mds <- setdiff(normal_mds , unique(combined.sct_geneset$Sample))
for (smp in normal_mds){
	print(smp)
	tmp <- readRDS(paste0('./Data/', smp, '/', smp, '_seurat_obj_norm.rds'))
	tmp <- as.data.frame(tmp@assays$SCT@data)
	tmp <- tmp[rownames(tmp) %in% genes$gene_name,]
	all_samples <- merge(all_samples, tmp, by=0)
	rownames(all_samples) <- all_samples$Row.names
	all_samples$Row.names <- NULL
	print(dim(all_samples))
}






all_samples_non_5q <- HL_signature(all_samples, get_genes('non5q'))
all_samples_non_5q$Row.names <- rownames(all_samples_non_5q)
all_samples_non_5q$Sample <- stringr::str_extract(all_samples_non_5q$Row.names, '^[A-Z0-9]+')

all_samples_non_5q$Case <- ifelse(
									all_samples_non_5q$Sample %in% c('GSM5460411', 'GSM5460412', 'GSM5460413') , 'Elder', 
							ifelse( all_samples_non_5q$Sample %in% c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303') , '5q', 
									'OtherMDS'))

# all_samples_non_5q$Case <- factor(all_samples_non_5q$Case, levels=c('Elder', 'OtherMDS', '5q'))

pdf('./Plots/Integrated/All_samples_HL_Signature_non5q_Lite.pdf')
get_histogram_all(all_samples_non_5q)
ggplot(all_samples_non_5q, aes(x= High_pondered, y = Case, fill=Case)) + ggridges::geom_density_ridges2() + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)
# get_histogram_all(all_samples_non_5q, 'SAMPLE')
get_ridges(all_samples_non_5q, 'Case')
# get_ridges(all_samples_non_5q)
ggplot(all_samples_non_5q, aes(x=Case, y=High_pondered, fill=Case)) + geom_boxplot(alpha=0.7) +
scale_fill_manual( values = destiny::cube_helix(length(unique(all_samples_non_5q$Case))) ) +
ggprism::theme_prism() + theme(legend.position='none') + 
ggsignif::geom_signif(comparisons = list(c('5q', 'Elder'), c('5q', 'OtherMDS'), c('OtherMDS', 'Elder')), map_signif_level=TRUE)
dev.off()






HL_signature_relu <- function(data, signature_df){
	data <- data[rownames(data) %in% signature_df$gene_name,]
    logpluspointfive <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(data, 2, logpluspointfive))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))
	zscore[zscore < 0] <- 0
	zscore$gene_name <- rownames(zscore)
    ann_markers <- merge(zscore, signature_df, by='gene_name')
	ann_markers <- ann_markers[ !duplicated(ann_markers$gene_name),]
	res <- ann_markers[, !colnames(ann_markers) %in% c('gene_name', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
	rownames(res) <- ann_markers$gene_name
	res <- setNames(as.data.frame(colSums(res[, !colnames(res) %in% c('gene_name', 'log2FoldChange', 'padj'), drop=FALSE], )), c('High_pondered'))

} 

all_samples_non_5q <- HL_signature_relu(all_samples, get_genes('non5q'))
all_samples_non_5q$Row.names <- rownames(all_samples_non_5q)
all_samples_non_5q$Sample <- stringr::str_extract(all_samples_non_5q$Row.names, '^[A-Z0-9]+')

all_samples_non_5q$Case <- ifelse(
									all_samples_non_5q$Sample %in% c('GSM5460411', 'GSM5460412', 'GSM5460413') , 'Elder', 
							ifelse( all_samples_non_5q$Sample %in% c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303') , '5q', 
									'OtherMDS'))

# all_samples_non_5q$Case <- factor(all_samples_non_5q$Case, levels=c('Elder', 'OtherMDS', '5q'))

pdf('./Plots/Integrated/All_samples_HL_Signature_non5q_Lite_RELU.pdf')
get_histogram_all(all_samples_non_5q)
ggplot(all_samples_non_5q, aes(x= High_pondered, y = Case, fill=Case)) + ggridges::geom_density_ridges2() + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)
# get_histogram_all(all_samples_non_5q, 'SAMPLE')
get_ridges(all_samples_non_5q, 'Case')
# get_ridges(all_samples_non_5q)
ggplot(all_samples_non_5q, aes(x=Case, y=High_pondered, fill=Case)) + geom_boxplot(alpha=0.7) +
scale_fill_manual( values = destiny::cube_helix(length(unique(all_samples_non_5q$Case))) ) +
ggprism::theme_prism() + theme(legend.position='none') + 
ggsignif::geom_signif(comparisons = list( c('OtherMDS', 'Elder')), map_signif_level=TRUE)
dev.off()


wilcox.test(all_samples_non_5q[all_samples_non_5q$Case == 'OtherMDS','High_pondered'],all_samples_non_5q[all_samples_non_5q$Case == 'Elder','High_pondered'])





signature_genes <- get_genes('non5q')

correlation_signature <- function(data, lfc){
	order <- sort(rownames(data))
	lfc <- lfc[lfc$gene_name %in% order,]
	data <- data[match(order, rownames(data)),]
	lfc <- lfc[match(order, lfc$gene_name),]
	signature <- data.frame(Signature = apply(data, 2, function(x) {cor(x, lfc$log2FoldChange)}))
	return(signature)
}


correlation_signature <- correlation_signature(all_samples, get_genes('non5q'))

pdf('./Plots/Integrated/All_samples_correlation_signature.pdf')
correlation_signature$Sample <- stringr::str_extract(rownames(correlation_signature), '^[A-Z0-9]+')
correlation_signature$Case <- ifelse(
									correlation_signature$Sample %in% c('GSM5460411', 'GSM5460412', 'GSM5460413') , 'Elder', 
							ifelse( correlation_signature$Sample %in% c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303') , '5q', 
									'OtherMDS'))

ggplot(correlation_signature, aes(x= Signature, y=Sample ,fill=Case)) + ggridges::geom_density_ridges2() + 
			theme_classic() + scale_fill_viridis_d(option='turbo', alpha=0.7)
dev.off()