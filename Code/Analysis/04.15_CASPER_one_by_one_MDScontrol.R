
library(argparse)
library(Seurat)
library(CaSpER)
library(ggplot2)
data("hg38_cytoband")

# TODO: launch them with this new samples
# drwxrwxr-x  2 sevastopol sevastopol 4.0K May  5 11:35 senior35210/
# drwxrwxr-x  2 sevastopol sevastopol 4.0K May  5 11:35 SMD29402/

# SAMPLE_NAMES_5Q <- c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303')
# SAMPLE_NAMES <- c("SMD21949","SMD29402","SMD34363","SMD35302","SMD35344","SMD36730","SMD38312","SMD41197","SMD_MFF" )


# SAMPLES TO RUN THIS SCRIPT  c('SMD34459', 'SMD35109', 'SMD37209', 'SMD35303', 'SMD29402')

parser <- ArgumentParser(description='Launc CASPER with control')
parser$add_argument('-S', '--SampleName', 
					default="SMD34459",
					type="character",
					help='Sample to process')
parser$add_argument('-C', '--ControlName', 
					default="SMD35344",
					type="character",
					help='Control to process compare: also can be SMD29402')
args <- parser$parse_args()
sample_name <- args$SampleName
control_name <-	args$ControlName


# controls <- SMD35344 or 
control <- readRDS(paste0('/home/sevastopol/data/gserranos/MDS/Data/',control_name,'/', control_name ,'_seurat_obj_norm.rds'))
control_counst <- as.matrix(control@assays$RNA@counts)



print(paste0('Processing: ', sample_name, ' With control: ', control_name))
if (sample_name == 'SMD29402'){
	data <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/SMD29402/SMD29402_seurat_obj_norm.rds')
	log.ge  <- as.matrix(data@assays$RNA@counts)
}else{
	if(!exists('all_seurat_integrated_sct')){
		all_seurat_integrated_sct <- readRDS('/home/sevastopol/data/gserranos/MDS/Data/all_seurat_integrated_sct.rds')
	}
	Idents(all_seurat_integrated_sct) <- 'Sample'
	log.ge <- subset(all_seurat_integrated_sct, idents = sample_name)
	log.ge  <- as.matrix(log.ge@assays$RNA@counts)
}

log.ge <- merge(log.ge, control_counst, by =0)
rownames(log.ge) <- log.ge$Row.names
log.ge <- log.ge[,-1]

# SMD35344 <- this is the control
control <- colnames(log.ge)[stringr::str_detect(colnames(log.ge), paste0('^', control_name))]

#https://github.com/akdess/CaSpER/issues/13

samples <- unique(stringr::str_extract(colnames(log.ge), '^[A-Z0-9^]+'))

genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
# log.ge <- log2(log.ge +1)
loh <- list()
path <- paste0('/home/sevastopol/data/gserranos/MDS/Data/BAM_files/', sample_name, '/')
baf_tmp <- readBAFExtractOutput(path=path, sequencing.type="bulk", suffix="baf")
loh[[sample_name]] <- baf_tmp[[paste0(sample_name, '.baf')]]
# path <- paste0('/home/sevastopol/data/gserranos/MDS/Data/BAM_files/', control_name, '/')
# baf_tmp <- readBAFExtractOutput(path=path, sequencing.type="bulk", suffix="baf")
# loh[[control_name]] <- baf_tmp[[paste0(control_name, '.baf')]]

loh.name.mapping <- data.frame(
					loh.name= factor(stringr::str_extract(colnames(log.ge), '^[A-Z0-9^]+'), levels = c(sample_name, control_name)),
					sample.name= colnames(log.ge))

object <- CreateCasperObject(raw.data=log.ge, 
							loh.name.mapping=loh.name.mapping, 
							sequencing.type="single-cell", 
							cnv.scale=3, loh.scale=3, 
							expr.cutoff=0.1, filter="median", matrix.type="raw",
							annotation=annotation, method="iterative", loh=loh, 
							control.sample.ids=control, cytoband=cytoband)

print(paste0('Casper Object created for: ', sample_name))

saveRDS(object, paste0('CASPER_OBJECT_',sample_name,'Vs',control_name,'.rds'))

## runCaSpER
final.objects <- runCaSpER(object, removeCentromere=T, cytoband=cytoband, method="iterative")
saveRDS(final.objects, paste0('./Data/CASPER/FINAL_CASPER_OBJECT_',sample_name,'Vs',control_name,'.rds'))
print(paste0('CNVs computed for: ', sample_name))
## summarize large scale events 

# title extractLargeScaleEvents()
# description generates coherent set of large scale CNV events using the pairwise comparison of all scales from BAF and expression signals
# param final.objects casper object
# param thr gamma threshold determining the least number of scales required to support 
# return final large scale event summary reported as a matrix 
finalChrMat <- extractLargeScaleEvents (final.objects, thr=0.75) 
saveRDS(finalChrMat, paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat',sample_name,'Vs',control_name,'.rds'))
print(paste0('Results ready for: ', sample_name))

plotHeatmap(object, fileName=paste0('./Plots/CASPER/Heatmap_',sample_name,'Vs',control_name,'.pdf'), cnv.scale= 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T)
plotLargeScaleEvent2(finalChrMat, fileName=paste0('./Plots/CASPER/LargeScaleEventsSummarized_',sample_name,'Vs',control_name,'.pdf')) 
plotBAFAllSamples(loh = final.objects[[1]]@loh.median.filtered.data,  fileName=paste0('./Plots/CASPER/LOH_',sample_name,'Vs',control_name,'.pdf'))
plotBAFOneSample (object, fileName=paste0('./Plots/CASPER/LOH_allScales_',sample_name,'Vs',control_name,'.pdf')) 
results <- extractMUAndCooccurence (finalChrMat, loh, loh.name.mapping)
plotMUAndCooccurence(results)

print(paste0('Succesfully processed: ', sample_name))

