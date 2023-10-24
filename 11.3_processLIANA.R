library(ggplot2)
 
 
data <- read.table('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/ResultsPerSample.csv', sep='\t', header=TRUE)
data_flt <-  data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]


# cehck all common between them
data_elder <- data_flt[grepl('^GSM', data_flt$Comparison),  ]
data_elder$Unique_Interaction <- paste(data_elder$source, data_elder$target, data_elder$ligand_complex, data_elder$receptor_complex, sep='_')

data_MDS <- data_flt[grepl('^SMD', data_flt$Comparison),  ]
data_MDS$Unique_Interaction <- paste(data_MDS$source, data_MDS$target, data_MDS$ligand_complex, data_MDS$receptor_complex, sep='_')


tmp_elder <- data_elder[, c('Comparison', 'Unique_Interaction')]
tmp_elder <- split(tmp_elder,tmp_elder$Comparison)

tmp_MDS <- data_MDS[, c('Comparison', 'Unique_Interaction')]
tmp_MDS <- split(tmp_MDS,tmp_MDS$Comparison)

split_by <- function(data){
	for (i in 1:length(data)){
		data[[i]] <- data[[i]][,2]
	}
return(data)
}
tmp_elder <- split_by(tmp_elder)
tmp_MDS <- split_by(tmp_MDS)


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
text_size = 2)
return(p)
}


get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)),
	order.by = "freq", empty.intersections = "on")
	# plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
	# 		nrow=2, align='hv', rel_heights = c(3,1),
	# 		rel_widths = c(2,3))
	return(plot)
}


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_venn.pdf')
print(get_gvenn(tmp_elder))
print(get_upset(tmp_MDS))
dev.off()


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_SCores.pdf')
tmp_scores <- data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]
tmp_scores <- tmp_scores[, c('Comparison', 'specificity_rank', 'magnitude_rank')]
ggplot(tmp_scores, aes(x=log(specificity_rank),  color=Comparison)) + geom_density()
ggplot(tmp_scores, aes(x=log(magnitude_rank),  color=Comparison)) + geom_density()
dev.off()









data <- read.table('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/ResultsPerSample_CT_and_Condition.csv', sep='\t', header=TRUE, row.names=1)

data_flt <-  data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]


# cehck all common between them
data_elder <- data_flt[grepl('^GSM', data_flt$Sample),  ]
data_elder$Unique_Interaction <- paste(data_elder$source, data_elder$target, data_elder$ligand_complex, data_elder$receptor_complex, sep='_')

data_MDS <- data_flt[grepl('^SMD', data_flt$Sample),  ]
data_MDS$Unique_Interaction <- paste(data_MDS$source, data_MDS$target, data_MDS$ligand_complex, data_MDS$receptor_complex, sep='_')


tmp_elder <- data_elder[, c('Sample', 'Unique_Interaction')]
tmp_elder <- split(tmp_elder,tmp_elder$Sample)

tmp_MDS <- data_MDS[, c('Sample', 'Unique_Interaction')]
tmp_MDS <- split(tmp_MDS,tmp_MDS$Sample)

split_by <- function(data){
	for (i in 1:length(data)){
		data[[i]] <- data[[i]][,2]
	}
return(data)
}
tmp_elder <- split_by(tmp_elder)
tmp_MDS <- split_by(tmp_MDS)


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
text_size = 2)
return(p)
}


get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)),
	order.by = "freq", empty.intersections = "on")
	# plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
	# 		nrow=2, align='hv', rel_heights = c(3,1),
	# 		rel_widths = c(2,3))
	return(plot)
}


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_venn.pdf')
print(get_gvenn(tmp_elder))
print(get_upset(tmp_MDS))
dev.off()


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_SCores.pdf')
tmp_scores <- data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]
tmp_scores <- tmp_scores[, c('Comparison', 'specificity_rank', 'magnitude_rank')]
ggplot(tmp_scores, aes(x=log(specificity_rank),  color=Comparison)) + geom_density()
ggplot(tmp_scores, aes(x=log(magnitude_rank),  color=Comparison)) + geom_density()
dev.off()



pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_SCores.pdf')
tmp <- data_elder
tmp$from <- stringr::str_extract(tmp$source, '[\\w]+(?=&)')
tmp$to <- stringr::str_extract(tmp$target, '[\\w]+(?=&)')
tmp$Interactive_cells <- paste(tmp$source, tmp$target, sep='_')
ggplot(tmp, aes(x=-log10(specificity_rank), y=-log10(magnitude_rank), color=from)) + 
geom_point(alpha=0.8) + theme_bw() + theme(legend.position='bottom')+
scale_color_manual(values=colorRampPalette(ggthemes::tableau_color_pal('Classic 20')(20))(length(unique(tmp$from))))

dev.off()

library(circlize)


mycolor <- viridis::viridis(length(unique(data_elder$source)), alpha = 1, begin = 0, end = 1, option = "D")

pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/Test_Chords.pdf')
chordDiagram(table(data_elder$source , data_elder$target), transparency = 0.5,   grid.col = mycolor)
dev.off()






data <- read.table('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/ResultsPerSample_CT_and_Condition.csv', sep='\t', header=TRUE, row.names=1)

data_flt <-  data[data$cellchat_pvals < 0.05 & data$cellphone_pvals < 0.05, ]


get_upset <- function(data_list){
	tmp <- UpSetR::fromList(data_list)
	plot <- UpSetR::upset(tmp, sets = c(names(data_list)),
	order.by = "freq", empty.intersections = "on")
	# plot <- cowplot::plot_grid(NULL, plot$Main_bar, plot$Sizes, plot$Matrix,
	# 		nrow=2, align='hv', rel_heights = c(3,1),
	# 		rel_widths = c(2,3))
	return(plot)
}


get_HM <- function(data, title){
	pheatmap::pheatmap(as.data.frame.matrix(table(data$source, data$target)), 
	cluster_rows=FALSE, cluster_cols=FALSE, main = title, color=RColorBrewer::brewer.pal(9,"Blues"))
}


# cehck all common between them
data_elder <- data_flt[grepl('^GSM', data_flt$Sample),  ]
data_elder$Unique_Interaction <- paste(stringr::str_extract(data_elder$source, '[\\w]+(?=&)'), stringr::str_extract(data_elder$target, '[\\w]+(?=&)'), data_elder$ligand_complex, data_elder$receptor_complex, sep='_')

data_MDS <- data_flt[grepl('^SMD', data_flt$Sample),  ]
data_MDS$Unique_Interaction <- paste(stringr::str_extract(data_MDS$source, '[\\w]+(?=&)'), stringr::str_extract(data_MDS$target, '[\\w]+(?=&)'), data_MDS$ligand_complex, data_MDS$receptor_complex, sep='_')

pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/UpsetNOfiltered.pdf')
tmp <- data_elder[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
tmp <- data_MDS[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
tmp <- rbind(data_elder[, c('Sample', 'Unique_Interaction')], data_MDS[, c('Sample', 'Unique_Interaction')])
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
dev.off()

pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/HistMagnitudes.pdf')
	hist(-log10(data_elder$magnitude_rank), breaks=50)
	abline(v=5, col='red', lwd=3, lty='dashed')
	hist(-log10(data_MDS$magnitude_rank), breaks=50)
	abline(v=5, col='red', lwd=3, lty='dashed')
dev.off()

data_MDS_flt   <- data_MDS[-log10(data_MDS$magnitude_rank) > 5, ]
data_elder_flt <- data_elder[-log10(data_elder$magnitude_rank) > 5, ]


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/UpsetFiltered.pdf')
tmp <- data_elder_flt[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
tmp <- data_MDS_flt[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
tmp <- rbind(data_elder_flt[, c('Sample', 'Unique_Interaction')], data_MDS_flt[, c('Sample', 'Unique_Interaction')])
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,2])
get_upset(tmp)
dev.off()

ground_healthy <- data_elder_flt$Unique_Interaction

WriteXLS::WriteXLS(data_elder_flt, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/LIANA_data_elder_filtered.xlsx'),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

# Remove the T lymphocytes as they are not present in the MDS


pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Heatmaps_MDS_with_healthy.pdf')
get_HM(data_MDS_flt, 'All Samples MDS with healthy')
get_HM(data_MDS_flt[data_MDS_flt$Sample == 'SMD34459', ], 'SMD34459 with healthy')
get_HM(data_MDS_flt[data_MDS_flt$Sample == 'SMD35109', ], 'SMD35109 with healthy')
get_HM(data_MDS_flt[data_MDS_flt$Sample == 'SMD35303', ], 'SMD35303 with healthy')
get_HM(data_MDS_flt[data_MDS_flt$Sample == 'SMD37209', ], 'SMD37209 with healthy')
dev.off()

data_MDS_flt_NoT <- data_MDS_flt[!data_MDS_flt$target %in% c('T&normal', 'T&del5q') & !data_MDS_flt$source %in% c('T&normal', 'T&del5q'), ]

tmp <- data_MDS_flt_NoT[, c('Sample', 'Unique_Interaction')]
tmp <- split(tmp, tmp$Sample)
tmp <- lapply(tmp, function(x) x[,'Unique_Interaction'])
all_shared_comm <- Reduce(intersect,tmp)

comms_core_MDS <- data_MDS_flt_NoT[data_MDS_flt_NoT$Unique_Interaction %in% all_shared_comm , ]
WriteXLS::WriteXLS(comms_core_MDS, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/comms_core_MDS.xlsx'),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


comms_core_MDS_no_Healthy <- comms_core_MDS[!comms_core_MDS$Unique_Interaction %in% ground_healthy, ]
WriteXLS::WriteXLS(comms_core_MDS_no_Healthy, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/comms_core_MDS_no_Healthy.xlsx'),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


tmp <- data_MDS[data_MDS$Sample == 'SMD34459', ]
table(tmp$source, tmp$target)
table(data_elder$source, data_elder$target)

pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Heatmaps_MDS.pdf')
get_HM(comms_core_MDS_no_Healthy, 'All Samples MDS no healthy')
get_HM(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Sample == 'SMD34459', ], 'SMD34459 no healthy')
get_HM(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Sample == 'SMD35109', ], 'SMD35109 no healthy')
get_HM(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Sample == 'SMD35303', ], 'SMD35303 no healthy')
get_HM(comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Sample == 'SMD37209', ], 'SMD37209 no healthy')
dev.off()

pdf('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Heatmaps_Elder.pdf')
get_HM(data_elder_flt, 'All Samples healthy')
get_HM(data_elder_flt[data_elder_flt$Sample == 'GSM5460411', ], 'GSM5460411 healthy')
get_HM(data_elder_flt[data_elder_flt$Sample == 'GSM5460412', ], 'GSM5460412 healthy')
get_HM(data_elder_flt[data_elder_flt$Sample == 'GSM5460413', ], 'GSM5460413 healthy')
dev.off()








# Novel interactions detected on MDS 5q- samples
table(unique(comms_core_MDS_no_Healthy$Unique_Interaction) %in% unique(data_elder$Unique_Interaction))
novel_mds5q_interactions <- comms_core_MDS_no_Healthy[!comms_core_MDS_no_Healthy$Unique_Interaction %in% data_elder$Unique_Interaction, ]
WriteXLS::WriteXLS(novel_mds5q_interactions, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/Novel_Interactions_in_MDS5q.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

# interactions detected on MDS 5q- samples but with low strength on healthy samples
intense_mds5q_interactions <- comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Unique_Interaction %in% data_elder$Unique_Interaction, ]
intense_mds5q_interactions$log10magnitude_rank <- -log10(intense_mds5q_interactions$magnitude_rank)
elder_lower_intesity <- data_elder[data_elder$Unique_Interaction %in% intense_mds5q_interactions$Unique_Interaction, ]
elder_lower_intesity$log10magnitude_rank <- -log10(elder_lower_intesity$magnitude_rank)

WriteXLS::WriteXLS(intense_mds5q_interactions, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/mds5q_intense_interactions.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(elder_lower_intesity, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/elder_lower_interactions.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)




# table(unique(intense_mds5q_interactions$Unique_Interaction) %in% unique(elder_lower_intesity$Unique_Interaction))

table(unique(data_elder_flt$Unique_Interaction) %in% unique(data_MDS$Unique_Interaction))
loss_mds5q_interactions <- data_elder_flt[!data_elder_flt$Unique_Interaction %in% data_MDS$Unique_Interaction, ]
loss_mds5q_interactions$log10magnitude_rank <- -log10(loss_mds5q_interactions$magnitude_rank)
WriteXLS::WriteXLS(loss_mds5q_interactions, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/loss_mds5q_interactions.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)



data_MDS$log10magnitude_rank <- -log10(data_MDS$magnitude_rank)
weak_mds <-data_MDS[data_MDS$log10magnitude_rank < 5, ]

elder_higher_intesity <- data_elder_flt[data_elder_flt$Unique_Interaction %in% weak_mds$Unique_Interaction, ]
elder_higher_intesity$log10magnitude_rank <- -log10(elder_higher_intesity$magnitude_rank)

mds_lower_intesity  <- weak_mds[weak_mds$Unique_Interaction %in% data_elder_flt$Unique_Interaction, ]
mds_lower_intesity$log10magnitude_rank <- -log10(mds_lower_intesity$magnitude_rank)
WriteXLS::WriteXLS(elder_higher_intesity, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/elder_higher_intesity.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)
WriteXLS::WriteXLS(mds_lower_intesity, 
			ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Plots/Liana/MDS_and_Elder/mds_lower_intesity.xlsx'),  
			col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)



# interseting scheme 

hard_interactions <- NULL
less_hard_interactions <- NULL
for (interaction in unique(novel_mds5q_interactions$Unique_Interaction)){
	tmp <- comms_core_MDS_no_Healthy[comms_core_MDS_no_Healthy$Unique_Interaction == interaction, ]
	tmp <- tmp[length(unique(tmp$Sample)) >=3,]
	if(nrow(tmp[grepl('normal', tmp$source) & grepl('normal', tmp$target),] !=0)){
		# some normal-normal interactions has been detected
		less_hard_interactions <- rbind(less_hard_interactions, tmp)
	}else{
		# no normal-normal interactions detected
		hard_interactions <- rbind(hard_interactions, tmp)
	}
}