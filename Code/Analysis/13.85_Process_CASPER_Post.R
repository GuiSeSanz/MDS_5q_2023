

library(Seurat)
library(CaSpER)
library(ggplot2)
data("hg38_cytoband")

data_post <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/all_seurat_Post_integrated_sct.rds')

deletion_percentages <- data.frame(Samples = c('FS-0634-post',	'FS-0406-post'), 
								Karyotype = c(0,0), 
								CASPER= c(0,0))

coords <- FetchData(data_post, vars = c("UMAP_1", "UMAP_2", "Sample"))

all_5q_depleted_cells <- c()
CONTROL = 'SMD29402'
for (sample_name in c('FS-0634-post', 'FS-0406-post')){	
	# sample_name <- 'SMD34459'
	print(sample_name)
	coords_smp <- coords[coords$Sample==sample_name, ]

	chrMat <- readRDS(paste0('./Data/CASPER/FINAL_CASPER_OBJECT_finalChrMat', sample_name,'Vs',CONTROL,'.rds'))
	casper_results <- melt(chrMat)
	casper_results$value2 <- "neutral"
	casper_results$value2[casper_results$value > 0] <- "amplification"
	casper_results$value2[casper_results$value < 0] <- "deletion"
	casper_results$value2 <- factor(casper_results$value2, levels = c("amplification", 
		"deletion", "neutral"))
	casper_results$X2 <- factor(casper_results$X2, levels = colnames(chrMat))
	casper_results$phenotype <- stringr::str_extract(casper_results$X1, '^[-A-Z0-9^]+')
	deletion_percentages[deletion_percentages$Samples == sample_name,'CASPER'] <- (prop.table(table(casper_results[casper_results$phenotype != 'SMD29402' & casper_results$ X2 == '5q', 'value2']))*100)[['deletion']]
	all_5q_depleted_cells <- c(all_5q_depleted_cells, as.character(casper_results[casper_results$X2 %in% c('5q') & casper_results$value2 == 'deletion','X1']))
	pdf(paste0('./Plots/CASPER/Results_CASPER_',sample_name,'Vs',CONTROL,'_umap.pdf'), height=5)
		coords_smp$deletion <- ifelse(rownames(coords_smp) %in% all_5q_depleted_cells, 'del5q', 'normal') 
		print(ggplot(coords_smp, aes(x=UMAP_1, y=UMAP_2, color=deletion)) + geom_point() + theme_classic())
	dev.off()

}


saveRDS(all_5q_depleted_cells, './Data/CASPER/all_5q_depleted_cells_POST.rds')



# all samples

get_gvenn <- function(data, scale_Y=TRUE){
	p <-  ggvenn::ggvenn(data,
	fill_color = destiny::cube_helix(length(data)),
	stroke_size = 0.4,
	show_percentage = TRUE,
	fill_alpha = 0.4,
	stroke_color = 'white',
	stroke_alpha = 1,
	stroke_linetype = 'solid',
	text_color = 'black',
	set_name_size = 4,
	text_size = 2,
	label_sep = ','
	)+ theme(plot.title = element_text(hjust = 0.5, family = "Helvetica", size = 7),
		plot.margin = unit(c(0, 0, 0, 0), "cm"),
		axis.text.y=element_blank(),
		axis.title.y=element_blank(),
		axis.ticks.y=element_blank())
	if(scale_Y){
	p <- p + scale_y_continuous(limits = c(-1, 1.5))
	}
	return(p)
}

all_5q_depleted_cells_COPYKAT <- readRDS('./Data/CopyKat/all_5q_CopyKat_depleted_cells_POST.rds')
all_5q_depleted_cells_CASPER <- readRDS('./Data/CASPER/all_5q_depleted_cells_POST.rds')
plot_list <- list()

for(Sample in  c('FS-0634-post', 'FS-0406-post')){
plot_list[[Sample]] <- get_gvenn(list(
						CopyKat = gsub('\\.', '-', all_5q_depleted_cells_COPYKAT[all_5q_depleted_cells_COPYKAT$Sample == Sample, 'Cell_id']),
						CASPER = all_5q_depleted_cells_CASPER[stringr::str_detect(all_5q_depleted_cells_CASPER, Sample)]))
}

pdf('./Plots/CopyKat_CASPER_5q_depleted_cells_POST.pdf', height=5)
cowplot::plot_grid(plotlist = plot_list, ncol = 2, labels=names(plot_list))
dev.off()



all_selected_cells <- intersect( gsub('\\.', '-',all_5q_depleted_cells_COPYKAT$Cell_id), all_5q_depleted_cells_CASPER)
table(stringr::str_extract(all_selected_cells, '^[\\w\\d-]+(?=_)'))
saveRDS(all_selected_cells, './Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')