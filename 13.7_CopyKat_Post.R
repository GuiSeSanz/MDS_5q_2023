library(Seurat)
library(copykat)
library(argparse)



# sample_name <- 'FS-0634-post'
sample_name <- 'FS-0406-post'
print(paste0('Processing: ', sample_name))

data_Sample <- readRDS(paste0('/home/tereshkova/data/gserranos/MDS/Data/', sample_name,'/',sample_name,'_seurat_obj_norm.rds' )) 

control_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

control_cells <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Shuffled_3_333_cells_per_elder_for_CASPER.tsv')
control_cells <- as.character(control_cells$V1)
control_data <- subset(control_data, cells=control_cells)


raw_data <- merge(data_Sample, control_data)

control_cells <- control_cells[control_cells %in% colnames(raw_data)]

working_dir = file.path('/home/tereshkova/data/gserranos/MDS/Data/CopyKat', sample_name)

if (file.exists(working_dir)){
    setwd(working_dir)
} else {
    dir.create(working_dir)
    setwd(working_dir)
}

print(paste0('Working directory: ', working_dir))
print(paste0('Number of cells: ', ncol(raw_data)))
print('Starting CopyKat analysis')

raw_data <- as.matrix(raw_data@assays$RNA@counts)
copykat.obj <- copykat(rawmat=raw_data, 
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

saveRDS(copykat.obj, paste0(working_dir, '/', sample_name, '_copykat_obj.rds'))