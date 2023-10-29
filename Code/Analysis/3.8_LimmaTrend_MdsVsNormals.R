# https://www.nature.com/articles/nmeth.4612
suppressPackageStartupMessages({
	library(argparse)
	library(Seurat)
	library(DESeq2)
	library(limma)
	library(edgeR)
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
parser$add_argument('-N', '--Normals', default="ELDER", type="character",
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
MDS_5q_data <- as.data.frame(all_seurat_integrated_sct@assays$RNA@counts)

genes_2_keep <- rownames(MDS_5q_data) 
for (sample in NORMAL_list){
	SAMPLE <- stringr::str_extract(sample, '[A-Z0-9]+(?=_young|_elderly)')
	print(SAMPLE)
	tmp <- readRDS(paste0('/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/', SAMPLE, '_seurat_obj_norm.rds'))
	tmp <- as.data.frame(tmp@assays$RNA@counts)
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




run_limmatrend <- function(counts, design_vector) {
	message("limmatrend")
	dge <- DGEList(counts, group = design_vector)
	dge <- calcNormFactors(dge)
	design <- model.matrix(~0+design_vector)
	y <- new("EList")
	y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
	fit <- lmFit(y, design = design)
	fit <- eBayes(fit, trend = TRUE, robust = TRUE)
	tT <- topTable(fit, n = Inf, adjust.method = "BH")
	tt <- toptable(fit, n = Inf, adjust.method = "BH")
	return(list(tt = tt, fit=fit, dge = dge))
}

# remove the genes with sd = 0
stand_dev <- apply(all_data,1, sd)
all_data <- all_data[!rownames(all_data) %in% names(which(stand_dev ==0)), ]


DESIGN <- data.frame(cells = colnames(all_data))
DESIGN$labels <- ifelse(stringr::str_extract(DESIGN$cells, '^[A-Z0-9]+') %in% samples_normal, 'Normal', 'MDS_5q')

### DE
my_design <- DESIGN$labels
tmp <- run_limmatrend(all_data, DESIGN$labels )


# pdf(paste0('./Plots/Test_',NORMAL,'.pdf'))
#   hist(tmp$tt$P.Value, 50)
#   hist(tmp$tt$adj.P.Val, 50)
#   limma::plotMDS(tmp$dge, col = as.numeric(as.factor(DESIGN$labels)), pch = 19)
#   limma::plotMD(tmp$fit)
# dev.off()



tmp2 <- tmp$tt[tmp$tt$adj.P.Val < 0.05,]
tmp2 <- tmp2[order(-abs(tmp2$logFC)),]


genes_expanded <- read.table( './Data/Annotation/5q13-33_TheGoodOne.txt', fill=TRUE, header=FALSE)
genes_expanded <- as.character(unique(genes_expanded$V5))
pathway <- list('MDS_5q'=genes_expanded)

stats_5qVsnormal <- tmp2$logFC
names(stats_5qVsnormal) <- rownames(tmp2)
fgseaRes <- fgsea::fgsea(pathways = pathway, 
				  stats    = stats_5qVsnormal ,
				  minSize  = 15,
				  nperm=1000,
				  maxSize  = 500)

saveRDS(fgseaRes, paste0('./Plots/fgsea_',NORMAL,'_trend.rds'))
pdf(paste0('./Plots/fgsea_',NORMAL,'_trend.pdf'))
fgsea::plotEnrichment(pathway[['MDS_5q']],stats_5qVsnormal)
dev.off()

pdf(paste0('./Plots/Volcano_', NORMAL,'_trend.pdf'))
EnhancedVolcano::EnhancedVolcano(tmp$tt,
	lab = genes_expanded,
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
	caption='P-value < 0.05; log2FC > 0.1',
	labSize = 6.0)
dev.off()
# all_seurat_integrated_sct <- readRDS(paste0(getwd(), '/Data/','all_seurat_integrated_sct.rds'))


