library(Seurat)
library(SeuratData)
library(SeuratDisk)




data_integrated <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

DefaultAssay(data_integrated) <- 'RNA'

SaveH5Seurat(data_integrated, filename = "/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_COUNTS.h5Seurat")
Convert('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_COUNTS.h5Seurat', dest = "h5ad")


write.csv(data_integrated@meta.data, sep='\t', quote=FALSE, row.names=TRUE,
file='/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_meta.csv')

write.table(data_integrated@reductions$umap@cell.embeddings, sep='\t', quote=FALSE, row.names=TRUE,
file='/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_umap.csv')



write.csv(as.matrix(data_integrated@assays$RNA@counts), sep='\t', quote=FALSE, row.names=TRUE,
file='/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_counts.csv')




data_integrated <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
# data_integrated@assays$integrated <- NULL
DefaultAssay(data_integrated)<-"RNA"
VariableFeatures(data_integrated) <-rownames(data_integrated)
SaveH5Seurat(data_integrated, filename = "/home/tereshkova/data/gserranos/MDS/Data/Elder_data_integrated_COUNTS.h5Seurat", overwrite =TRUE)
Convert('/home/tereshkova/data/gserranos/MDS/Data/Elder_data_integrated_COUNTS.h5Seurat', dest = "h5ad", overwrite =TRUE)





library(anndata)
my_anndata_object <- as.h5ad(data_integrated)
