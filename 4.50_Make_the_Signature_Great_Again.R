

library(ggplot2)
library(Seurat)
options(stringsAsFactors=FALSE)

libsize_nomrmalization <- function(counts){
	libSizes <- colSums(counts)
	sizeFactors <- libSizes/mean(libSizes)
	normCounts <- counts/sizeFactors
	log_norm_counts <- log2(normCounts+1)
	return(log_norm_counts)
}

get_genes <- function(sample_name_filter){
	# sample_name_filter <- elder or non5q
	de_files <- c(  'HSC5qvsElderly.tsv',
				'HSCs5qvsnon5q.tsv',
				'MEPs5qvsElderly.csv',
				'MEPs5qvsnon5q.tsv',
				'CMPs5qvsElderly.tsv',
				'CMPs5qvsnon5q.tsv', 
				'GMPs5qvsElderly.tsv',
				'GMPs5qvsnon5q.tsv')
	data_folder <- '/home/tereshkova/data/gserranos/MDS/Data/Annotation'
	genes <- data.frame(X=NULL, log2FoldChange=NULL)
	for (fl in de_files){
		if (stringr::str_detect(fl, stringr::fixed(sample_name_filter, ignore_case=TRUE))){
			tmp <- read.table(paste(data_folder, fl, sep = '/'), sep=',', header=TRUE)
			tmp <- tmp[!is.na(tmp$padj),]
			tmp <- tmp[tmp$padj<0.05, c('X', 'log2FoldChange')]
			genes <- rbind(genes, tmp)
		}
	}
	gene_2_reannotate <- which(stringr::str_detect(genes$X, '^ENSG+'))
	tmp <- setNames(genes[gene_2_reannotate,], c('gene_id', 'log2FoldChange'))
	tmp <- merge(tmp, gene_annotation, by='gene_id')
	tmp <- tmp[, c('gene_name', 'log2FoldChange')]
	genes <- setNames(genes, c('gene_name', 'log2FoldChange'))
	genes <- genes[!genes$gene_name %in%  genes$gene_name[gene_2_reannotate],]
	genes <- rbind(genes,tmp )
	return(genes)
}

HL_signature <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
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

HL_signature_relu <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
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

correlation_signature <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
	data <- data[rownames(data) %in% signature_df$gene_name, ]
	signature_df <- signature_df[signature_df$gene_name %in% rownames(data),]
	signature_df <- signature_df[!duplicated(signature_df$gene_name),]
	data <- data[order(rownames(data)),]
	signature_df <- signature_df[order(signature_df$gene_name),]
	signature <- data.frame(Signature = apply(data, 2, function(x) {cor(x, signature_df$log2FoldChange)}))
	return(signature)
}


pearson_rank <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
	gene_names <- rownames(data)
	data <- as.data.frame(apply(data, 2, scale))
	ranked_data <- as.data.frame(apply(data, 2, rank))
	rownames(ranked_data) <- gene_names
	ranked_data <- ranked_data[rownames(ranked_data) %in% signature_df$gene_name,]
	signature_df <- signature_df[signature_df$gene_name %in% rownames(ranked_data),]
	signature_df <- signature_df[!duplicated(signature_df$gene_name),]
	ranked_data <- ranked_data[order(rownames(ranked_data)),]
	signature_df <- signature_df[order(signature_df$gene_name),]
	signature_df_rank <- rank(signature_df$log2FoldChange)
	res <- as.data.frame(apply( ranked_data, 2, function(x){cor(signature_df_rank, x, method='pearson')}))
	return(res)
}


pearson_rank_relu <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
	signature_df <- signature_df[signature_df$log2FoldChange > 0,]
	gene_names <- rownames(data)
	data <- as.data.frame(apply(data, 2, scale))
	ranked_data <- as.data.frame(apply(data, 2, rank))
	rownames(ranked_data) <- gene_names
	ranked_data <- ranked_data[rownames(ranked_data) %in% signature_df$gene_name,]
	signature_df <- signature_df[signature_df$gene_name %in% rownames(ranked_data),]
	signature_df <- signature_df[!duplicated(signature_df$gene_name),]
	ranked_data <- ranked_data[order(rownames(ranked_data)),]
	signature_df <- signature_df[order(signature_df$gene_name),]
	signature_df_rank <- rank(signature_df$log2FoldChange)
	res <- as.data.frame(apply( ranked_data, 2, function(x){cor(signature_df_rank, x, method='pearson')}))
	return(res)
}


HL_signature_zscore_filter <- function(data, signature_df=NULL, refenrence = 'non5q'){
	if( is.null(signature_df)){
		signature_df <-  get_genes(refenrence)
	}
	# data <- data[rownames(data) %in% signature_df$gene_name,]
    logpluspointfive <- function(x) {log(x + 0.5)}
    l <- as.data.frame(apply(data, 2, logpluspointfive))
    # head(l[,1:5])
    zscore <- as.data.frame(scale(l))
	zscore$gene_name <- rownames(zscore)
	zscore <- zscore[rownames(zscore) %in% signature_df$gene_name,]
	signature_df <- signature_df[signature_df$gene_name %in% rownames(zscore),]
	signature_df <- signature_df[!duplicated(signature_df$gene_name),]
	ann_markers <- merge(zscore, signature_df, by='gene_name')
	ann_markers$log2FoldChange <- scale(ann_markers$log2FoldChange)
	ann_markers <- ann_markers[ !duplicated(ann_markers$gene_name),]
	res <- ann_markers[, !colnames(ann_markers) %in% c('gene_name', 'log2FoldChange', 'padj')] * ann_markers[, 'log2FoldChange']
	rownames(res) <- ann_markers$gene_name
	res <- setNames(as.data.frame(colSums(res[, !colnames(res) %in% c('gene_name', 'log2FoldChange', 'padj'), drop=FALSE], )), c('High_pondered'))

} 





normal_K_non5qMDS<-c('A112','A114','A116','A119','A120','A121','A122','A123','A124','A125','A126','A127','A130','A131','A132','A134','A136','A142','A143','A147','A148','A151','A153','A154','A155','A158','A161','A164','A165','A170','A173','A176','A178','A180','A182','A183','A184','A186','A187','A188','A189','A190','A191','A193','A195','A196','A197','A198','A204','A208')
special_K_non5qMDS<-c('A128','A138','A139','A141','A149','A150','A156','A166','A167','A168','A171','A172','A175','A177','A179','A181','A185','A192','A200','A205','A206','A207','A203')
healthy<-c('A115','A133','A135','A137','A146','A159','A162','A169')
MDS_5q<-c('A117','A118','A129','A140','A163','A194','A209','A145')



public_data <- read.table('./Data/Normal_Data/GSE114922_Count_table.txt', header=TRUE, row.names=1)
colnames(public_data) <- stringr::str_extract(colnames(public_data), '(A[0-9]{3})')
public_data <- libsize_nomrmalization(public_data)
gene_annotation <- setNames(read.table('./Data/Annotation/Homo_sapiens.GRCh38.Gene_name_ID.tsv', header=FALSE, sep='\t'), c('gene_id' , 'gene_name'))
public_data <- merge(public_data, gene_annotation, by.x=0, by.y='gene_id')

# remove bullshit names
# genes_2_keep <-  names(which(table(public_data$gene_name ) == 1))
# public_data <- public_data[public_data$gene_name %in% genes_2_keep,]

# signature_genes <- get_genes('non5q')
# signature_genes$gene_name <- as.character(signature_genes$gene_name)

interesting_samples <- c(normal_K_non5qMDS, special_K_non5qMDS, healthy, MDS_5q)
sum(!interesting_samples %in% colnames(public_data)) == 0

# FILTER THE GENES
# public_data <- public_data[public_data$gene_name %in% signature_genes$gene_name,]
duplicated_gene_names <- as.character(public_data$gene_name[ which(duplicated(public_data$gene_name))])
public_data <- public_data[!public_data$gene_name %in% duplicated_gene_names,]
rownames(public_data) <- public_data$gene_name
public_data <- public_data[, !colnames(public_data) %in% c('Row.names', 'gene_name')]


Signature         <- setNames(HL_signature(public_data, refenrence = 'elder'), c('Signature'))
SignatureRelu     <- setNames(HL_signature_relu(public_data, refenrence = 'elder'), c('SignatureRelu'))
Signature_cor     <- setNames(correlation_signature(public_data, refenrence = 'elder'), c('Signature_cor'))
Pearson_Rank      <- setNames(pearson_rank(public_data, refenrence = 'elder'), c('Pearson_Rank'))
Zscore_LogFC      <- setNames(HL_signature_zscore_filter(public_data, refenrence = 'elder'), c('Zscore_LogFC'))
Pearson_Rank_Relu <- setNames(pearson_rank_relu(public_data, refenrence = 'elder'), c('pearson_rank_relu'))



Signatures <- data.frame(   Signature = Signature, 
							SignatureRelu = SignatureRelu, 
							Signature_cor = Signature_cor, 
							Pearson_Rank = Pearson_Rank, 
							Pearson_Rank_Relu = Pearson_Rank_Relu,
							Zscore_LogFC = Zscore_LogFC
							)

Signatures$Sample <- rownames(Signatures)
Signatures$Case <- ifelse(Signatures$Sample %in% normal_K_non5qMDS, 'Normal_MDS', ifelse(Signatures$Sample %in% special_K_non5qMDS, 'Special_Karyotype', ifelse(Signatures$Sample %in% healthy, 'Healthy', 'MDS_5q')))

saveRDS(Signatures, './Data/Signatures_BULK.rds')

pdf('./Plots/Test_Bulk_elder.pdf',)
ggplot( reshape2::melt(Signatures), aes(x=Case, y=value, fill =Case)) + geom_boxplot(alpha=0.7) + 
theme_classic() + scale_fill_viridis_d(option='turbo') + 
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom', 
axis.title.x=element_blank(), axis.title.y=element_blank()) + 
facet_wrap(~variable, nrow=2, scale='free')
dev.off()








de_files <- c(  'HSC5qvsElderly.tsv',
				'HSCs5qvsnon5q.tsv',
				'MEPs5qvsElderly.csv',
				'MEPs5qvsnon5q.tsv',
				'CMPs5qvsElderly.tsv',
				'CMPs5qvsnon5q.tsv', 
				'GMPs5qvsElderly.tsv',
				'GMPs5qvsnon5q.tsv')

data_folder <- '/home/tereshkova/data/gserranos/MDS/Data/Annotation'

combined.sct_geneset <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
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
	if (smp == normal_mds[1]){
		normals <- tmp
	}
	else{
		tmp <- tmp[rownames(tmp) %in% rownames(normals),]
		normals <- normals[rownames(normals) %in% rownames(tmp),]

		normals <- cbind(normals, tmp)
	}

}

genes_2_keep <- intersect(rownames(normals), rownames(all_samples))
all_samples_mds_5q_elder <- cbind(all_samples[rownames(all_samples) %in% genes_2_keep,], normals[rownames(normals) %in% genes_2_keep,])



Pearson_Rank      <- setNames(pearson_rank(all_samples_mds_5q_elder , refenrence='elder'),               c('Pearson_Rank'))
Zscore_LogFC      <- setNames(HL_signature_zscore_filter(all_samples_mds_5q_elder , refenrence='elder'), c('Zscore_LogFC'))
Signature         <- setNames(HL_signature(all_samples_mds_5q_elder , refenrence='elder'),               c('Signature'))
SignatureRelu     <- setNames(HL_signature_relu(all_samples_mds_5q_elder , refenrence='elder'),          c('SignatureRelu'))
Signature_cor     <- setNames(correlation_signature(all_samples_mds_5q_elder , refenrence='elder'),      c('Signature_cor'))
Pearson_Rank_Relu <- setNames(pearson_rank_relu(all_samples_mds_5q_elder, refenrence = 'elder'),         c('pearson_rank_relu'))

# saveRDS(Pearson_Rank, './Data/Pearson_Rank_signature.rds')
# saveRDS(Zscore_LogFC, './Data/Zscore_LogFC_signature.rds')
# saveRDS(Signature, './Data/Signature_signature.rds')
# saveRDS(SignatureRelu, './Data/SignatureRelu_signature.rds')
# saveRDS(Signature_cor, './Data/Signature_cor_signature.rds')




Signatures <- data.frame(	Signature = Signature, 
							SignatureRelu = SignatureRelu, 
							Signature_cor = Signature_cor, 
							Pearson_Rank = Pearson_Rank, 
							Pearson_Rank_Relu= Pearson_Rank_Relu,
							Zscore_LogFC = Zscore_LogFC)

Signatures$Sample <- stringr::str_extract(rownames(Signatures), '^[A-Z0-9]+')
Signatures$Case <- ifelse(
						Signatures$Sample %in% c('GSM5460411', 'GSM5460412', 'GSM5460413') , 'Elder', 
				ifelse( Signatures$Sample %in% c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303') , '5q', 
						'OtherMDS'))

pdf('./Plots/Test_Sc_elder.pdf.pdf',)
ggplot( reshape2::melt(Signatures), aes(x=Case, y=value, fill =Case)) + geom_boxplot(alpha=0.7) + 
theme_classic() + scale_fill_viridis_d(option='turbo') + 
theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom', 
axis.title.x=element_blank(), axis.title.y=element_blank()) + 
facet_wrap(~variable, nrow=2, scale='free')
dev.off()






get_annotation_sofia_04 <- function(){
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
	return(annotation_Sofia_04)
}


all_seurat_integrated_sct <- readRDS('./Data/all_seurat_integrated_sct.rds')

coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
# annotation resolution is 0.4
coords <- merge(coords, setNames(as.data.frame(all_seurat_integrated_sct$integrated_snn_res.0.4), c('Cluster_04')), by=0)
ann_04 <- get_annotation_sofia_04()
coords$ClusterName <- apply(tmp, 1, FUN=function(x) ann_04[[as.character(x[['Cluster_04']])]])

cluster_names <- c()
for (cl in as.character(coords$Cluster_04)){
	cluster_names <- c(cluster_names , ann_04[[as.character(cl)]])

}
coords$ClusterName <- cluster_names


ridges_patients <- ggplot(Signatures, aes(x=Pearson_Rank, y= Sample, fill = Case)) +  ggridges::geom_density_ridges2(alpha = 0.8) + theme_classic() +
scale_fill_manual(values =destiny::cube_helix(length(unique(Signatures$Case)))) + 
theme(legend.position='left', axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

tmp <- merge(coords, Signatures, by.x='Row.names', by.y=0)

ridges_cell_types <- ggplot(tmp, aes(x=Pearson_Rank, y= ClusterName, fill = ClusterName)) +  ggridges::geom_density_ridges2(alpha = 0.8) + theme_classic() +
scale_fill_viridis_d(option='turbo') + 
theme(legend.position='none', axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))

pdf('./Plots/signature_pearson.pdf')
print(ridges_patients)
print(ridges_cell_types)
dev.off()