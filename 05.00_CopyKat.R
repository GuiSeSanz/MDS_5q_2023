library(Seurat)
library(copykat)
library(argparse)

parser <- ArgumentParser(description='Launch CopyKat with fixed mixed controls')
parser$add_argument('-S', '--SampleName', 
					default="SMD34459",
					type="character",
					help='Sample to process from: SMD34459, SMD35109, SMD34459, SMD37209')

args <- parser$parse_args()

sample_name <- args$SampleName


# It might take a while to run a dataset with more than 10,000 single cells. 
# It is suggested to run one sample at a time. Combining different sample would pick up batch effects between samples.
# Final final note, CopyKAT had difficulty in predicting tumor and normal cells in the cases of pediatric and liquid 
# tumors that have a few CNAs. CopyKAT provides two ways to bypass this to give certain output instead of being dead staright: 
# 1) input a vector of cell names of known normal cells from the same dataset 2) or try to search for T cells.


# GSM5460411 GSM5460412 GSM5460413   SMD34459   SMD35109   SMD35303   SMD37209 
#       8360      19112      17070       5352      14759      12811      13850


print(paste0('Processing: ', sample_name))

combined.sct_geneset <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/integrated_5q_and_elder_OnlyChr5_sct.rds')
raw_data <- as.data.frame(combined.sct_geneset@assays$RNA@counts)

control_cells <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Shuffled_3_333_cells_per_elder_for_CASPER.tsv')
control_cells <- as.character(control_cells$V1)

raw_data_lite <- cbind( raw_data[, stringr::str_detect(colnames(raw_data), paste0('^', sample_name))] , 
						raw_data[, control_cells] )


working_dir = file.path('/home/tereshkova/data/gserranos/MDS/Data/CopyKat', sample_name)

if (file.exists(working_dir)){
    setwd(working_dir)
} else {
    dir.create(working_dir)
    setwd(working_dir)
}

print(paste0('Working directory: ', working_dir))
print(paste0('Number of cells: ', ncol(raw_data_lite)))
print('Starting CopyKat analysis')

copykat.obj <- copykat(rawmat=raw_data_lite, 
						id.type="S", 
						ngene.chr=5, 
						win.size=25, 
						KS.cut=0.1, 
						sam.name=sample_name, 
						distance="euclidean", 
						norm.cell.names= control_cells, 
						output.seg="TRUE", 
						plot.genes="TRUE", 
						genome="hg20",
						n.cores=60)


# saveRDS(copykat.obj, paste0(sample_name, '_copykat_obj.rds'))



# pred.test <- data.frame(copykat.test$prediction)
# pred.test <- pred.test[-which(pred.test$copykat.pred=="not.defined"),]  ##remove undefined cells
# CNA.test <- data.frame(copykat.test$CNAmat)

# annotation <- data.frame(Sample = colnames(raw_data_lite))
# annotation$barcode <- stringr::str_extract(annotation$Sample, '(?<=[-_])([ACTG]+)(?=[-_])')
# annotation$Type <- ifelse(stringr::str_detect(annotation$Sample, '^GSM'), 'Elder', '5q')

# predictions <- read.table('/home/sevastopol/data/gserranos/MDS/test_copykat_prediction.txt', header=TRUE)
# predictions_annotated <- merge(annotation, predictions, by.x="Sample", by.y="cell.names")


# table(predictions_annotated$Type, predictions_annotated$copykat.pred)





# all_seurat_integrated_sct <- readRDS('./Data/all_seurat_integrated_sct.rds')
# coords <- as.data.frame(all_seurat_integrated_sct@reductions$umap@cell.embeddings)
# casper_results_smp <- casper_results(samples[1])


# casper_results <- function(samples){
# 	all_results <- data.frame(X1=NULL, X2=NULL, value=NULL, value2=NULL, phenotype=NULL)
# 	for (sample_name in samples){
# 		print(paste0('Processing the',sample_name))
# 		coords_smp <- coords[stringr::str_detect(coords$Row.names, sample_name), ]
# 		chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'.rds'))
# 		casper_results <- reshape::melt(chrMat)
# 		casper_results$value2 <- "neutral"
# 		casper_results$value2[casper_results$value > 0] <- "amplification"
# 		casper_results$value2[casper_results$value < 0] <- "deletion"
# 		casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
# 			"deletion", "neutral"))
# 		casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
# 		casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[A-Z0-9^]+')
# 		all_results <- rbind(all_results, casper_results)
# 	}
# 	return(all_results)
# }

# get_venn <- function(data, percent = TRUE){
# 	p <- ggvenn::ggvenn(data,
# 			fill_color = destiny::cube_helix(length(data)),
# 			stroke_size = 0.4,
# 			show_percentage = percent,
# 			fill_alpha = 0.4,
# 			stroke_color = 'white',
# 			stroke_alpha = 1,
# 			stroke_linetype = 'solid',
# 			text_color = 'black',
# 			set_name_size = 4,
# 			text_size = 3.5,
# 			label_sep = ','
# 			)+ theme(plot.title = element_text(hjust = 0.5))
# 	return(p)
# }


# casper_results_smp <- casper_results_smp[casper_results_smp$phenotype == 'SMD34459' & casper_results_smp$value2 != 'neutral',]

# aneuploid_samples <- predictions_annotated[predictions_annotated$Type == '5q' & predictions_annotated$copykat.pred == 'aneuploid', ]

# data <- list()
# data[['anueploid_copykat']] <- aneuploid_samples$Sample
# data[['casper_no_neutral']] <- casper_results_smp$X1

# pdf('./Plots/CasperVscopykat.pdf')
# get_venn(data)
# dev.off()