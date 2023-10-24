
suppressPackageStartupMessages({
	library(argparse)
	library(Seurat)
	library(DESeq2)
	library(limma)
})
# read all the counts
young_list <- c('GSM5460406_young1_filtered_feature_bc_matrix_counts.tsv',
				# 'GSM5460407_young2_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460408_young3_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460409_young4_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460410_young5_filtered_feature_bc_matrix_counts.tsv')
elder_list <- c('GSM5460411_elderly1_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460412_elderly2_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460413_elderly3_filtered_feature_bc_matrix_counts.tsv')


parser <- ArgumentParser(description='Perform the DE analysis with the normals')
parser$add_argument('-N', '--Normals', default="YOUNG", type="character",
                    help='YOUNG or ELDER')
args <- parser$parse_args()
NORMAL <- args$Normals
print(NORMAL)
# set normal to young or elder
if (NORMAL == 'YOUNG') {
	NORMAL_list <- young_list
} else if (NORMAL == 'ELDER') {
	NORMAL_list <- elder_list
} else {
	stop("NORMAL must be 'YOUNG' or 'ELDER'")
}




all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))
MDS_5q_data <- as.data.frame(all_seurat_integrated_sct@assays$SCT@data)

genes_2_keep <- rownames(MDS_5q_data) 
for (sample in NORMAL_list){
	SAMPLE <- stringr::str_extract(sample, '[A-Z0-9]+(?=_young|_elderly)')
	print(SAMPLE)
	tmp <- readRDS(paste0('/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/', SAMPLE, '_seurat_obj_norm.rds'))
	tmp <- as.data.frame(tmp@assays$SCT@data)
	if (sample == NORMAL_list[1]) {
		norm_data <- tmp
	} else {
		norm_data <- base::merge(norm_data, tmp, by=0, all=TRUE)
		rownames(norm_data) <- norm_data$Row.names
		norm_data <- norm_data[,-1]
	}
	genes_2_keep <- intersect(genes_2_keep, rownames(tmp))
}

MDS_5q_data <- MDS_5q_data[genes_2_keep,]
norm_data <- norm_data[genes_2_keep,]

samples_normal <- unique(stringr::str_extract(colnames(norm_data), '^[A-Z0-9^_]+'))
samples_5q     <- unique(stringr::str_extract(colnames(MDS_5q_data), '^[A-Z0-9^_]+'))

all_data <- base::merge(MDS_5q_data, norm_data, by=0)
rownames(all_data) <- all_data$Row.names
all_data <- all_data[,-1]

perform_DE <- function(data, group){
  ### DE
  # create the linear model
  fit_tmp <- lmFit(data, group)
  # model correction
  fit_tmp <- eBayes(fit_tmp)
  # results <- topTable(fit_tmp, n=Inf)
  x <- paste0('MDS_5q', '-', 'Normal')
  contrast_mat_tmp <- makeContrasts(contrasts=x, levels= c('MDS_5q', 'Normal'))
  fit2_tmp <- contrasts.fit(fit_tmp, contrast_mat_tmp)
  fit2_tmp <- eBayes(fit2_tmp)
  tmp   <- topTable(fit2_tmp, adjust="fdr", n=Inf)
  tmp$gene_name <- rownames(tmp)
  
#   tmp[tmp$P.Value == 0, 'P.Value'] <- 1.445749e-281
#   tmp[tmp$adj.P.Val == 0, 'adj.P.Val'] <- 1.445749e-281
  
  return(tmp)
}

# remove the genes with sd = 0
stand_dev <- apply(all_data,1, sd)
all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]


DESIGN <- data.frame(cells = colnames(all_data))
DESIGN$labels <- ifelse(stringr::str_extract(DESIGN$cells, '^[A-Z0-9]+') %in% samples_normal, 'Normal', 'MDS_5q')

design_tmp <- as.matrix(DESIGN[, 'labels'])
design_tmp[design_tmp == 'Normal'] <-1
design_tmp[design_tmp == 'MDS_5q'] <-0
design_tmp <- model.matrix(~0+as.factor(design_tmp[,1]))
colnames(design_tmp) <- c('MDS_5q', 'Normal')
rownames(design_tmp) <- colnames(all_data)

# test <- limma::removeBatchEffect(all_data, DESIGN$labels, desing = design_tmp)


### DE
tmp <- perform_DE(all_data, design_tmp )

tmp2 <- tmp[tmp$adj.P.Val < 0.05,]
tmp2 <- tmp2[order(-abs(tmp2$logFC)),]


genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
genes_expanded <- as.character(unique(genes_expanded$V5))
pathway <- list('MDS_5q'=genes_expanded)

stats_5qVsnormal <- tmp2$logFC
names(stats_5qVsnormal) <- tmp2$gene_name
fgseaRes <- fgsea::fgsea(pathways = pathway, 
                  stats    = stats_5qVsnormal ,
                  minSize  = 15,
				  nperm=1000,
                  maxSize  = 500)

pdf(paste0('./Plots/fgsea_',NORMAL,'.pdf'))
fgsea::plotEnrichment(pathway[['MDS_5q']],stats_5qVsnormal)
dev.off()

pdf(paste0('./Plots/Volcano_', NORMAL,'.pdf'))
EnhancedVolcano::EnhancedVolcano(tmp,
    # lab = reasults_CD8$GeneID,
	lab=NA,
    x = 'logFC',
    y = 'adj.P.Val',
    title = '',
	subtitle = "5q versus Young",
	legendPosition='right',
    pCutoff = 0.05,
    FCcutoff = 0.1,
    pointSize = 1.5,
	colAlpha = 0.8,
	col=c('grey30', 'grey30', 'grey30', 'red2'),
	caption='P-value < 0.05; log2FC > 1',
    labSize = 6.0)
dev.off()
# all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))


