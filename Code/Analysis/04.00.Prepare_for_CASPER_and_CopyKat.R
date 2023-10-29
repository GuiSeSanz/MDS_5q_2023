library(Seurat)



# We need to select 333 samples from each elder (3 * 3333 = 9999 cells) 

data_path <- '/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/RawData/'
elder_list <- c('GSM5460411_elderly1_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460412_elderly2_filtered_feature_bc_matrix_counts.tsv',
				'GSM5460413_elderly3_filtered_feature_bc_matrix_counts.tsv')


elder_samples <- list()
for (elder_file in elder_list){
	print(elder_file)
	elder <- stringr::str_extract(elder_file,'(?<=_)[a-z]+[\\d]{1}')
	SAMPLE <- stringr::str_extract(elder_file, '[A-Z0-9]+(?=_young|_elderly)')
	tmp <- readRDS(paste0('/home/sevastopol/data/gserranos/MDS/Data/Normal_Data/', SAMPLE, '_seurat_obj_norm.rds'))
	elder_samples[[SAMPLE]] <- colnames(tmp)
}



mean(unlist(lapply(elder_samples, length))) # more or less 15K cells, so we will subsample 5k per sample

selected_samples <- lapply(elder_samples, function(x) x[sample(length(x), 3333)])
write.table(unlist(selected_samples), 'Shuffled_3_333_cells_per_elder_for_CASPER.tsv', sep='\t', row.names=FALSE, col.names=FALSE, quote =FALSE)

