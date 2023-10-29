
library(Seurat)
library(CaSpER)
library(ggplot2)
data("hg38_cytoband")
library(argparse)



parser <- ArgumentParser(description='Launch CASPER with fixed  control')
parser$add_argument('-S', '--SampleName', 
					default="FS-0634-post",
					type="character",
					help='Sample to process from: FS-0634-post // FS-0406-post')

args <- parser$parse_args()

sample_name <- args$SampleName

# sample_name <- 'FS-0634-post'
# sample_name <- 'FS-0406-post'

control_name <- 'SMD29402'
control <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/',control_name,'/', control_name ,'_seurat_obj_norm.rds'))
control_counst <- as.matrix(control@assays$RNA@counts)



seurat_data <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', sample_name,'/',sample_name,'_seurat_obj_norm.rds' ))
log.ge  <- as.matrix(seurat_data@assays$RNA@counts)


log.ge <- merge(log.ge, control_counst, by =0)
rownames(log.ge) <- log.ge$Row.names
log.ge <- log.ge[,-1]

# SMD29402 <- this is the control
control <- colnames(log.ge)[stringr::str_detect(colnames(log.ge), paste0('^', control_name))]

#https://github.com/akdess/CaSpER/issues/13

samples <- unique(stringr::str_extract(colnames(log.ge), '^[\\w\\d-]+(?=_)'))

genes <- rownames(log.ge)
annotation <- generateAnnotation(id_type="hgnc_symbol", genes=genes, centromere=centromere, ishg19 = F)
log.ge <- log.ge[match( annotation$Gene,rownames(log.ge)) , ]
rownames(log.ge) <- annotation$Gene
# log.ge <- log2(log.ge +1)

loh <- list()
path <- paste0('/home/tereshkova/data/gserranos/MDS/Data/BAM_files/', sample_name, '/')
baf_tmp <- readBAFExtractOutput(path=path, sequencing.type="bulk", suffix="baf")
loh[[sample_name]] <- baf_tmp[[paste0(sample_name, '.baf')]]
path <- paste0('/home/tereshkova/data/gserranos/MDS/Data/BAM_files/', control_name, '/')
baf_tmp <- readBAFExtractOutput(path=path, sequencing.type="bulk", suffix="baf")
loh[[control_name]] <- baf_tmp[[paste0(control_name, '.baf')]]

loh.name.mapping <- data.frame(
					loh.name= factor(stringr::str_extract(colnames(log.ge), '^[\\w\\d-]+(?=_)'), levels = c(sample_name, control_name)),
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

# plotBAFOneSample (object, fileName="LohPlotsAllScales.pdf") 