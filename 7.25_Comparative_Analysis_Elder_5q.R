

library(Seurat)
library(ggplot2)



########################
# rename Prog
########################

annotation_Sofia_08_Elder <- list()
annotation_Sofia_08_Elder[["0_0"   ]] <- 'LMPP'
annotation_Sofia_08_Elder[["0_1"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["0_2"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["0_3"   ]] <- 'CLP'
annotation_Sofia_08_Elder[["0_4"   ]] <- 'GMP'
annotation_Sofia_08_Elder[["1_0"   ]] <- 'LMPP'
annotation_Sofia_08_Elder[["1_1"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["1_2"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["1_3"   ]] <- 'CLP'
annotation_Sofia_08_Elder[["1_4"   ]] <- 'GMP'
annotation_Sofia_08_Elder[["2_0"   ]] <- 'LMPP'
annotation_Sofia_08_Elder[["2_1"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["2_2"   ]] <- 'HSC'
annotation_Sofia_08_Elder[["2_3"   ]] <- 'CLP'
annotation_Sofia_08_Elder[["2_4"   ]] <- 'GMP'
annotation_Sofia_08_Elder[["3"   ]] <- 'Granulocyte'
annotation_Sofia_08_Elder[["4"   ]] <- 'EarlyErythroid'
# annotation_Sofia_08_Elder[["5_0" ]] <- 
# annotation_Sofia_08_Elder[["6_0" ]] <- 
# annotation_Sofia_08_Elder[["6_1" ]] <- 
annotation_Sofia_08_Elder[["7"   ]] <- 'MEP'
annotation_Sofia_08_Elder[["8_0_0" ]] <- 'LMPP'
annotation_Sofia_08_Elder[["8_0_1" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_0_2" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_0_3" ]] <- 'CLP'
annotation_Sofia_08_Elder[["8_0_4" ]] <- 'GMP'
annotation_Sofia_08_Elder[["8_1_0" ]] <- 'LMPP'
annotation_Sofia_08_Elder[["8_1_1" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_1_2" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_1_3" ]] <- 'CLP'
annotation_Sofia_08_Elder[["8_1_4" ]] <- 'GMP'
annotation_Sofia_08_Elder[["8_2_0" ]] <- 'LMPP'
annotation_Sofia_08_Elder[["8_2_1" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_2_2" ]] <- 'HSC'
annotation_Sofia_08_Elder[["8_2_3" ]] <- 'CLP'
annotation_Sofia_08_Elder[["8_2_4" ]] <- 'GMP'
annotation_Sofia_08_Elder[["9_0"   ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["9_1"   ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["9_2"   ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["9_3"   ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["9_4"   ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["10_0"  ]] <- 'LMPP'
annotation_Sofia_08_Elder[["10_1"  ]] <- 'HSC'
annotation_Sofia_08_Elder[["10_2"  ]] <- 'HSC'
annotation_Sofia_08_Elder[["10_3"  ]] <- 'CLP'
annotation_Sofia_08_Elder[["10_4"  ]] <- 'GMP'
annotation_Sofia_08_Elder[["11"  ]] <- 'MEP'
# annotation_Sofia_08_Elder[["12_0"]] <- 
annotation_Sofia_08_Elder[["12_1"]] <- 'Granulocyte'
annotation_Sofia_08_Elder[["12_2"]] <- 'MK_Prog'
# annotation_Sofia_08_Elder[["13_0"]] <- 
# annotation_Sofia_08_Elder[["13_1"]] <- 
annotation_Sofia_08_Elder[["14_0"]] <- 'Basophil'
# annotation_Sofia_08_Elder[["14_1"]] <- 
annotation_Sofia_08_Elder[["14_2"]] <- 'Basophil'
annotation_Sofia_08_Elder[["15"  ]] <- 'MK_Prog'
annotation_Sofia_08_Elder[["16_0"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["16_1"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["16_2"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["16_3"  ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["16_4"  ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["17_0"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["17_1"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["17_2"  ]] <- 'LateErythroid'
annotation_Sofia_08_Elder[["17_3"  ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["17_4"  ]] <- 'EarlyErythroid'
annotation_Sofia_08_Elder[["18_0"]] <- 'DendriticCell'
annotation_Sofia_08_Elder[["18_1"]] <- 'DendriticCell'
annotation_Sofia_08_Elder[["19"  ]] <- 'Monocytes'
# annotation_Sofia_08_Elder[["20_0"]] <- 
annotation_Sofia_08_Elder[["21_0"]] <- 'pro-B'
# annotation_Sofia_08_Elder[["22"  ]] <- 
annotation_Sofia_08_Elder[["23"  ]] <- 'LateErythroid'


# PrePost
annotation_Nerea_08_PrePost <- list()
annotation_Nerea_08_PrePost[["0"   ]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["1_0"   ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["1_1"   ]] <- 'GMP'
annotation_Nerea_08_PrePost[["1_2"   ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["1_3"   ]] <- 'CLP'
annotation_Nerea_08_PrePost[["1_4"   ]] <- 'HSC'
annotation_Nerea_08_PrePost[["2"   ]] <- 'LateErythroid'
annotation_Nerea_08_PrePost[["3"   ]] <- 'Monocytes'
annotation_Nerea_08_PrePost[["4"   ]] <- 'LateErythroid'
annotation_Nerea_08_PrePost[["5"   ]] <- 'MK_Prog'
annotation_Nerea_08_PrePost[["6"   ]] <- 'MEP'
annotation_Nerea_08_PrePost[["7"   ]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["8"   ]] <- 'Granulocyte'
annotation_Nerea_08_PrePost[["9_0_0" ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["9_0_1" ]] <- 'GMP'
annotation_Nerea_08_PrePost[["9_0_2" ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["9_0_3" ]] <- 'CLP'
annotation_Nerea_08_PrePost[["9_0_4" ]] <- 'HSC'
annotation_Nerea_08_PrePost[["9_1_0" ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["9_1_1" ]] <- 'GMP'
annotation_Nerea_08_PrePost[["9_1_2" ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["9_1_3" ]] <- 'CLP'
annotation_Nerea_08_PrePost[["9_1_4" ]] <- 'HSC'
annotation_Nerea_08_PrePost[["10_0"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["10_1"  ]] <- 'GMP'
annotation_Nerea_08_PrePost[["10_2"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["10_3"  ]] <- 'CLP'
annotation_Nerea_08_PrePost[["10_4"  ]] <- 'HSC'
annotation_Nerea_08_PrePost[["11_0"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["11_1"  ]] <- 'GMP'
annotation_Nerea_08_PrePost[["11_2"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["11_3"  ]] <- 'CLP'
annotation_Nerea_08_PrePost[["11_4"  ]] <- 'HSC'
annotation_Nerea_08_PrePost[["12"  ]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["13_0"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["13_1"  ]] <- 'GMP'
annotation_Nerea_08_PrePost[["13_2"  ]] <- 'LMPP'
annotation_Nerea_08_PrePost[["13_3"  ]] <- 'CLP'
annotation_Nerea_08_PrePost[["13_4"  ]] <- 'HSC'
annotation_Nerea_08_PrePost[["14"  ]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["16"  ]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["17_1"]] <- 'EarlyErythroid'
annotation_Nerea_08_PrePost[["17_2"]] <- 'MEP'
annotation_Nerea_08_PrePost[["18"  ]] <- 'Basophil'
# annotation_Nerea_08_PrePost[["19_0"]] <- '????'
# annotation_Nerea_08_PrePost[["19_1"]] <- '????'
annotation_Nerea_08_PrePost[["21"  ]] <- 'MK_Prog'




# MDS 5q
annotation_Sofia_08_MDS5q <- list()
annotation_Sofia_08_MDS5q[['0_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['0_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['0_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['0_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['0_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['0_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['0_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['1_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['1_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['1_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['1_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['1_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['1_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['1_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['2_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['2_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['2_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['2_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['2_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['2_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['2_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['3_0'   ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['3_1'   ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['3_2'   ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['3_3'   ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['3_4'   ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['4'   ]] <- 'pro-B'
annotation_Sofia_08_MDS5q[['5_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['5_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['5_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['5_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['5_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['5_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['5_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['6'   ]] <- 'Monocytes'
annotation_Sofia_08_MDS5q[['7_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['7_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['7_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['7_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['7_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['7_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['7_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['8_0'   ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['8_1'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['8_2'   ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['8_3'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['8_4'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['8_5'   ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['8_6'   ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['9'   ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['10'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['11'  ]] <- 'Granulocyte'
annotation_Sofia_08_MDS5q[['12_0']] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['12_1']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['12_2']] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['13'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['14'  ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['15_0'  ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['15_1'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['15_2'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['15_3'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['15_4'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['15_5'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['15_6'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['16'  ]] <- 'MK_Prog'
# annotation_Sofia_08_MDS5q[['17_0']] <- '??'
annotation_Sofia_08_MDS5q[['17_1_0']] <- 'LMPP'
annotation_Sofia_08_MDS5q[['17_1_1']] <- 'HSC'
annotation_Sofia_08_MDS5q[['17_1_2']] <- 'HSC'
annotation_Sofia_08_MDS5q[['17_1_3']] <- 'GMP'
annotation_Sofia_08_MDS5q[['17_1_4']] <- 'CLP'
annotation_Sofia_08_MDS5q[['17_1_5']] <- 'CLP'
annotation_Sofia_08_MDS5q[['17_1_6']] <- 'GMP'
annotation_Sofia_08_MDS5q[['18_0'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['18_1'  ]] <- 'LateErythroid'
annotation_Sofia_08_MDS5q[['18_2'  ]] <- 'EarlyErythroid'
annotation_Sofia_08_MDS5q[['18_3'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['18_4'  ]] <- 'MEP'
annotation_Sofia_08_MDS5q[['19'  ]] <- 'pro-B'
annotation_Sofia_08_MDS5q[['20'  ]] <- 'Basophil'
annotation_Sofia_08_MDS5q[['21_0'  ]] <- 'LMPP'
annotation_Sofia_08_MDS5q[['21_1'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['21_2'  ]] <- 'HSC'
annotation_Sofia_08_MDS5q[['21_3'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['21_4'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['21_5'  ]] <- 'CLP'
annotation_Sofia_08_MDS5q[['21_6'  ]] <- 'GMP'
annotation_Sofia_08_MDS5q[['22'  ]] <- 'T'
# annotation_Sofia_08_MDS5q[['23']] <- 
annotation_Sofia_08_MDS5q[['24'  ]] <- 'Monocytes'
# annotation_Sofia_08_MDS5q[['25']] <- 
annotation_Sofia_08_MDS5q[['26'  ]] <- 'DendriticCell'
annotation_Sofia_08_MDS5q[['27'  ]] <- 'LateErythroid'
# annotation_Sofia_08_MDS5q[['28']] <- 
# annotation_Sofia_08_MDS5q[['29']] <- 

annotation_lists <- list()
annotation_lists[['5qSamples']]     <- annotation_Sofia_08_MDS5q
annotation_lists[['pre_post_5q']]   <- annotation_Nerea_08_PrePost
annotation_lists[['elder_Samples']] <- annotation_Sofia_08_Elder


get_volcano <- function(results, CType, sub= "5q Vs Elder"){
	p <-  EnhancedVolcano::EnhancedVolcano(results,
				lab = results$gene_name,
				x = 'avg_logFC', captionLabSize=7,
				y = 'p_val_adj',col=c('black', 'black', 'black', 'red3'),
				subtitle =sub, axisLabSize = 5, title= CType,
				pCutoff = 0.05, labSize=2,
				FCcutoff = 0.2) +
				theme(legend.position='none', text = element_text(family = "Helvetica", size = 2),
				axis.title= element_text(family = "Helvetica", size = 7))
	return(p)
}

get_tornado <- function(results, title){
	p <- ggplot(results, aes(x = NES, y = reorder(pathway, NES), fill = is_pos)) + 
			geom_bar(stat = "identity", position = "identity") + theme_classic() +
			ylab('Pathways in C2')+
			theme(legend.position="none", text = element_text(family = "Helvetica", size = 7), axis.text.y= element_text(size = 5),) + 
			scale_fill_manual(values=c(Neg='#006e90', Pos='#cc2936')) + ggtitle(title)
	return(p)
}

get_heatmap <- function(matrix){
 		anno_col <- data.frame(Sample= stringr::str_extract(colnames(matrix), '(^[A-Z0-9]+)(?=:)'), Case= stringr::str_extract(colnames(mat), '(?<=:)[\\w]+'))
        rownames(anno_col) <- colnames(matrix)
        colors_ann <- list(
			Sample = setNames(c("#f94144","#f3722c","#f8961e","#f9844a","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590","#277da1"),unique(metadata$replicate)),
			Case = setNames(c('#2a9d8f', '#f4a261'), unique(metadata$label))
		)
			pheatmap::pheatmap(matrix, 
			scale = "row",
            treeheight_row=0, treeheight_col=0,
            cluster_col = FALSE,cluster_row = TRUE,
            # color = viridis::viridis(50),
        #    gaps_row = gap_indexes,
            # cutree_rows = 2, 
            # cutree_cols = 3,
            fontsize_row = 4,
            legend=TRUE,
        #    annotation_row= anno_row,
            annotation_col=anno_col,
            annotation_colors=colors_ann,
            # annotation_legend=FALSE,
            annotation_names_row = FALSE,
            annotation_names_col = FALSE,
        #    annotation_row = anno_row,
            show_colnames = FALSE,
            show_rownames = TRUE,
            silent=TRUE,
            border_color='NA')
}

filter_and_annotate <- function(data, annotation){
    clusters <- setNames(as.data.frame(data$integrated_snn_res.0.8_sub), c('Cluster'))
    clusters <- clusters[paste(clusters$Cluster) %in% names(annotation),, drop=FALSE]
    clusters$Cluster_names <- apply(clusters, 1, function(x){ annotation[as.character(x[['Cluster']])][[1]] })
    Idents(data) <- 'integrated_snn_res.0.8_sub'
    data <- subset(data, idents= c(names(annotation)))
    data$Cluster_names <- clusters$Cluster_names
    return(data)
}


pathways_c2 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c2.cp.v7.5.1.symbols.gmt')




# if (file.exists('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')){
# 	sc_data_MDS5q_subseted <- readRDS('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')
# }else{

# 	sc_data_MDS5q <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_integrated_sct_subcluster.rds')
# 	annotation_Sofia_08_MDS5q <- list()
# 	annotation_Sofia_08_MDS5q[['0' ]] <- 'HSC'
# 	annotation_Sofia_08_MDS5q[['1' ]] <- 'HSC'
# 	annotation_Sofia_08_MDS5q[['2' ]] <- 'HSC'
# 	annotation_Sofia_08_MDS5q[['3' ]] <- 'EarlyErythroid'
# 	annotation_Sofia_08_MDS5q[['4' ]] <- 'pro-B'
# 	annotation_Sofia_08_MDS5q[['5' ]] <- 'LMPP'
# 	annotation_Sofia_08_MDS5q[['6' ]] <- 'Monocytes'
# 	annotation_Sofia_08_MDS5q[['7' ]] <- 'LMPP'
# 	annotation_Sofia_08_MDS5q[['8' ]] <- 'GMP'
# 	annotation_Sofia_08_MDS5q[['9' ]] <- 'EarlyErythroid'
# 	annotation_Sofia_08_MDS5q[['10']] <- 'LateErythroid'
# 	annotation_Sofia_08_MDS5q[['11']] <- 'Granulocyte'
# 	annotation_Sofia_08_MDS5q[['12_0']] <- 'LateErythroid'
# 	annotation_Sofia_08_MDS5q[['12_1']] <-  'LateErythroid'
# 	annotation_Sofia_08_MDS5q[['12_2']] <- 'LateErythroid'
# 	annotation_Sofia_08_MDS5q[['13']] <- 'EarlyErythroid'
# 	annotation_Sofia_08_MDS5q[['14']] <- 'EarlyErythroid'
# 	annotation_Sofia_08_MDS5q[['15']] <- 'CLP'
# 	annotation_Sofia_08_MDS5q[['16']] <- 'MEP'
# 	annotation_Sofia_08_MDS5q[['17_1']] <- 'GMP'
# 	annotation_Sofia_08_MDS5q[['18']] <- 'EarlyErythroid'
# 	annotation_Sofia_08_MDS5q[['19']] <- 'pro-B'
# 	annotation_Sofia_08_MDS5q[['20']] <- 'Basophil'
# 	annotation_Sofia_08_MDS5q[['21']] <- 'LMPP'
# 	annotation_Sofia_08_MDS5q[['22']] <- 'T'
# 	annotation_Sofia_08_MDS5q[['24']] <- 'Monocytes'
# 	annotation_Sofia_08_MDS5q[['26']] <- 'DendriticCell'
# 	annotation_Sofia_08_MDS5q[['27']] <- 'LateErythroid'

# 	Idents(sc_data_MDS5q) <- 'integrated_snn_res.0.8_sub'
# 	clusters <- setNames(as.data.frame(sc_data_MDS5q$integrated_snn_res.0.8_sub), c('Cluster'))
# 	coords   <- as.data.frame(sc_data_MDS5q@reductions$umap@cell.embeddings)
# 	Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78')
# 	names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell")
# 	MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
# 	clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_08_MDS5q),, drop=FALSE]
# 	coords_ann <- merge(clusters, coords, by=0)
# 	coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
# 	coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_08_MDS5q[as.character(x[['Cluster']])][[1]] })
# 	sc_data_MDS5q_subseted <- subset(sc_data_MDS5q, idents = c(names(annotation_Sofia_08_MDS5q)))
# 	rownames(coords_ann) <- coords_ann$Row.names
# 	coords_ann <- coords_ann[coords_ann$Row.names %in% colnames(sc_data_MDS5q_subseted),]
# 	sc_data_MDS5q_subseted$Cluster_names <- coords_ann[, 'Cluster_names', drop=FALSE]


# 	all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') 
# 	cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')
# 	real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)

# 	sc_data_MDS5q_subseted$cell_5q <- ifelse(colnames(sc_data_MDS5q_subseted) %in% real_5q_cells, 'del5q', 
# 										ifelse(colnames(sc_data_MDS5q_subseted) %in% 
# 										base::union(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id), 'maybe5q',
# 										 'normal'))

# # saveRDS(sc_data_MDS5q_subseted, './Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')
# }



# if(file.exists('/home/tereshkova/data/gserranos/MDS_std/all_seurat_Elder_integrated_sct_subcluster.rds')){
# 	elder_integrated_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_Elder_integrated_sct_subcluster.rds')
# }else{
# 	elder_integrated <-  readRDS(paste0(getwd(), '/Data/','all_seurat_Elder__integrated_sct_subcluster.rds'))
# 	annotation_Sofia_08 <- list()
# 	annotation_Sofia_08[["0"   ]] <- 'HSC'
# 	annotation_Sofia_08[["1"   ]] <- 'HSC'
# 	annotation_Sofia_08[["2"   ]] <- 'HSC'
# 	annotation_Sofia_08[["3"   ]] <- 'Granulocyte'
# 	annotation_Sofia_08[["4"   ]] <- 'EarlyErythroid'
# 	annotation_Sofia_08[["7"   ]] <- 'EarlyErythroid'
# 	annotation_Sofia_08[["8_0" ]] <- 'CLP'
# 	annotation_Sofia_08[["8_1" ]] <- 'LMPP'
# 	annotation_Sofia_08[["8_2" ]] <- 'LMPP'
# 	annotation_Sofia_08[["9"   ]] <- 'LateErythroid'
# 	annotation_Sofia_08[["10"  ]] <- 'GMP'
# 	annotation_Sofia_08[["11"  ]] <- 'MEP'
# 	annotation_Sofia_08[["12_1"]] <- 'Granulocyte'
# 	annotation_Sofia_08[["12_2"]] <- 'MEP'
# 	annotation_Sofia_08[["14_0"]] <- 'Basophil'
# 	annotation_Sofia_08[["14_2"]] <- 'Basophil'
# 	annotation_Sofia_08[["15"  ]] <- 'MEP'
# 	annotation_Sofia_08[["16"  ]] <- 'EarlyErythroid'
# 	annotation_Sofia_08[["17"  ]] <- 'LateErythroid'
# 	annotation_Sofia_08[["18_0"]] <- 'DendriticCell'
# 	annotation_Sofia_08[["18_1"]] <- 'DendriticCell'
# 	annotation_Sofia_08[["19"  ]] <- 'Monocytes'
# 	annotation_Sofia_08[["21_0"]] <- 'pro-B'
# 	annotation_Sofia_08[["23"  ]] <- 'LateErythroid'

# 	coords   <- as.data.frame(elder_integrated@reductions$umap@cell.embeddings)

# 	clusters <- setNames(as.data.frame(elder_integrated$integrated_snn_res.0.8_sub), c('Cluster'))
# 	clusters <- clusters[paste(clusters$Cluster) %in% names(annotation_Sofia_08),, drop=FALSE]
# 	coords_ann <- merge(clusters, coords, by=0)
# 	coords_ann$Sample <- stringr::str_extract(coords_ann$Row.names, '^[A-Z0-9]+')
# 	coords_ann$Cluster_names <- apply(coords_ann, 1, function(x){ annotation_Sofia_08[as.character(x[['Cluster']])][[1]] })
# 	elder_integrated_subseted <- subset(elder_integrated, idents= c(names(annotation_Sofia_08)))
# 	elder_integrated_subseted$Cluster_names <- coords_ann$Cluster_names
# 	ElderSample_colors <- c( '#264653', '#2a9d8f', '#e9c46a')
# 	names(ElderSample_colors) <- unique(coords_ann$Sample)

# 	# saveRDS(elder_integrated_subseted, '/home/tereshkova/data/gserranos/MDS_std/all_seurat_Elder_integrated_sct_subcluster.rds')
# }


sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
sc_data_MDS5q_subseted <- filter_and_annotate(sc_data_MDS5q_subseted, annotation_lists[['5qSamples']])


elder_integrated_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
elder_integrated_subseted <- filter_and_annotate(elder_integrated_subseted, annotation_lists[['elder_Samples']])



Cluster_colors <- c('#A55E34', '#C6B2D3', '#D0342B', '#8E221A', '#2E6B34', '#BBDE93', '#AECDE1', '#3C77AF', '#ED9E9B', '#DA913D', '#821851', '#643F95', '#DBBE78', '#3e5722', '#03071e', '#003049')
names(Cluster_colors) <- c("HSC","EarlyErythroid","pro-B","LMPP","Monocytes","GMP","LateErythroid","Granulocyte","CLP","MEP","Basophil","T","DendriticCell", 'MK_Prog', '????', 'HSC/LMPP')
MDS_sample_colors <- setNames(c('#85a6b2', '#495867', '#577399', '#bdd5ea'), c("SMD34459", "SMD35109", "SMD35303", "SMD37209"))
ElderSample_colors <- setNames(c( '#264653', '#2a9d8f', '#e9c46a'), c("GSM5460411", "GSM5460412", "GSM5460413"))

# Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
# sc_data_MDS5q_subseted <- subset(sc_data_MDS5q_subseted, idents = "del5q")
# Compare populations

pdf('./Plots/Comparison__pops_MDS5qVsElder.pdf')
sg_pal <- c(MDS_5q="#BF1363", Elder="#39A6A3")

populations <- rbind(setNames(as.data.frame(elder_integrated_subseted$Cluster_names), 'Cluster_names'), 
					 setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'Cluster_names'))

populations$Sample <- stringr::str_extract(rownames(populations),'(^[A-Z0-9]+)(?=_|-)')
pop_prop <- table(populations$Cluster_names, populations$Sample)
pop_prop <- apply(pop_prop, 2, FUN=function(x) (x/sum(x))*100)
pop_prop <- as.data.frame(reshape2::melt(pop_prop))
pop_prop$Case <- ifelse(stringr::str_detect(pop_prop$Var2, '^SMD'), 'MDS_5q', 'Elder')

print(ggplot(pop_prop, aes(x=Var1, y=value, fill=Case)) + 
geom_boxplot(outlier.shape = NA) + scale_fill_manual(values=sg_pal) + ggprism::theme_prism() + 
theme(axis.text.x= element_text(angle = 45, vjust = 1, hjust=1), legend.position='bottom', 
legend.title=element_blank(), axis.title.x=element_blank()) + ylab(paste('Percent of cells')) +
ggtitle('Elder Vs all cells in the MDS5q samples'))

meta <- setNames(cbind( as.data.frame(sc_data_MDS5q_subseted$Sample),
					cbind(
						as.data.frame(sc_data_MDS5q_subseted$Cluster_names),
						as.data.frame(sc_data_MDS5q_subseted$cell_5q))) , 
				c('Sample', 'Cluster_names', 'Cell_5q'))

plot_list <- list()
meta <- meta[meta$Cell_5q %in% c('normal', 'del5q'),]
sg_pal <- c(real_5q="#BF1363", normal="#39A6A3")
for (smp in unique(meta$Sample)){
	mt <- meta[meta$Sample == smp, ]
	mt <- table(mt$Cell_5q, mt$Cluster_names)
	mt <- apply(mt, 2, FUN=function(x) (x/sum(x))*100)
	mt <- as.data.frame(reshape2::melt(mt))
	mt$Var1 <- ifelse(mt$Var1 == 'del5q', 'real_5q', 'normal')

	plot_list[[smp]] <- ggplot(mt, aes(x=Var2, y=value, color=Var1, group=Var1)) + 
	geom_point() + geom_line() + scale_color_manual(values=sg_pal) + ggprism::theme_prism() + 
	theme(axis.text.x= element_text(angle = 45, vjust = 1, hjust=1), legend.position='none', 
	legend.title=element_blank(), axis.title.x=element_blank()) + ylab(paste('Percent of cells')) +
	ggtitle(smp)
}
lgd <- cowplot::get_legend(plot_list[[smp]] +theme(legend.position='bottom'))
print(cowplot::plot_grid(
			cowplot::plot_grid(plotlist=plot_list, nrow=2), 
			lgd, rel_heights=c(0.95, 0.05),  nrow=2))
dev.off()



################################################
# 5q Vs ELDER
################################################
library(Libra)

sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
sc_data_MDS5q_subseted <- filter_and_annotate(sc_data_MDS5q_subseted, annotation_lists[['5qSamples']])
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'


elder_integrated_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
elder_integrated_subseted <- filter_and_annotate(elder_integrated_subseted, annotation_lists[['elder_Samples']])



# LIBRA uses the raw counts
expression <- as.data.frame(elder_integrated_subseted@assays$RNA@counts)
sc_data_MDS5q_subseted <- subset(sc_data_MDS5q_subseted, idents = "del5q")
tmp <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
rownames(tmp) <- gsub('-', '\\.', rownames(tmp))
del_5qcells <- names(sc_data_MDS5q_subseted$cell_5q[sc_data_MDS5q_subseted$cell_5q  == 'del5q'])
tmp <- tmp[, colnames(tmp) %in% del_5qcells ]

# set dashes as dots
rownames(tmp)        <- gsub('-', '\\.', rownames(tmp))
rownames(expression) <- gsub('-', '\\.', rownames(expression))

expression <- merge(expression, tmp, by=0)

rownames(expression) <- expression$Row.names
expression <- expression[, -1]

# merged_seurat <- merge(elder_integrated_subseted, sc_data_MDS5q_subseted)

metadata <- rbind(setNames(data.frame(elder_integrated_subseted$Cluster_names),'cell_type') , 
				  setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type'))
metadata$label <- ifelse(stringr::str_detect(rownames(metadata), '^SMD'), 'MDS_5q', 'Elder')
metadata$replicate <- stringr::str_extract(rownames(metadata),'(^[A-Z0-9]+)(?=_|-)')

metadata$label <- factor(metadata$label, levels = c('MDS_5q', 'Elder'))
metadata <- metadata[rownames(metadata) %in% colnames(expression), ]

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_flt <- DE[DE$p_val_adj < 0.05,]

all_results_psuedo <- split( DE_flt , f = DE_flt$cell_type )
saveRDS(all_results_psuedo ,'./Data/Results_DE_5qVsElder.rds')

pdf('./Plots/Volcano_And_tornado_Pseudo_5qVsElder.pdf')
for (cell_type in sort(names(all_results_psuedo))){
	print(cell_type)
	tmp <- all_results_psuedo[[cell_type]]
	tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
	rownames(tt) <- tt$gene_name
	print(get_volcano(tt, cell_type))
	tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_logFC, tmp$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
	print(get_tornado(results_c2, cell_type))
	}
}
dev.off()


WriteXLS::WriteXLS(all_results_psuedo, ExcelFileName=paste0('./Plots/PseudoBulk_Results_5qVsElder.xlsx'), SheetNames = names(all_results_psuedo),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

matrices = to_pseudobulk(as.matrix(expression), meta = metadata)

saveRDS(matrices, './Data/matrices_pseudobulk_5qVsElder.rds')

# library(dplyr)
# pdf(paste0('./Plots/PseudoBulk_Results_5qVsElder_HM.pdf'))
# for (cell_type in sort(names(all_results_psuedo))){
# 	print(cell_type)
# 	tmp <- all_results_psuedo[[cell_type]]
# 	mat <- matrices[[cell_type]]
# 	mat <- mat[rownames(mat) %in% tmp$gene,]
# 	targets = data.frame(group_sample = colnames(mat)) %>%
# 		mutate(group = gsub(".*\\:", "", group_sample))
# 	group = gsub(".*\\:", "", colnames(mat))
# 	design = model.matrix(~ group, data = targets)
# 	dge = edgeR::DGEList(counts = mat, group =  group) %>%
# 			edgeR::calcNormFactors(method = 'TMM') %>%
# 			edgeR::estimateDisp(design)
# 	norm_mat <- edgeR::cpm(dge)
# 	norm_mat <- norm_mat[rownames(norm_mat) %in% tmp$gene,]
# 	norm_mat = t(scale(t(norm_mat)))
# 	column_ha = ComplexHeatmap::HeatmapAnnotation(Sample = ifelse(stringr::str_detect(colnames(norm_mat), ':Elder$'), 'Elder', 'MDS5q'))
# 	print(ComplexHeatmap::Heatmap(norm_mat, top_annotation = column_ha))
# }
# dev.off()



	
################################################
# Non5q Vs ELDER
################################################

expression <- as.data.frame(elder_integrated_subseted@assays$RNA@counts)


sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
sc_data_MDS5q_subseted <- filter_and_annotate(sc_data_MDS5q_subseted, annotation_lists[['5qSamples']])

tmp <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
rownames(tmp) <- gsub('-', '\\.', rownames(tmp))
Non_5q_cells <- names(sc_data_MDS5q_subseted$cell_5q[sc_data_MDS5q_subseted$cell_5q  == 'normal'])
tmp <- tmp[, colnames(tmp) %in% Non_5q_cells ]

# set dashes as dots
rownames(tmp)        <- gsub('-', '\\.', rownames(tmp))
rownames(expression) <- gsub('-', '\\.', rownames(expression))


expression <- merge(expression, tmp, by=0)
rownames(expression) <- expression$Row.names
expression <- expression[, -1]



# merged_seurat <- merge(elder_integrated_subseted, sc_data_MDS5q_subseted)

metadata <- rbind(setNames(data.frame(elder_integrated_subseted$Cluster_names),'cell_type') , 
				  setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type'))
metadata$label <- ifelse(stringr::str_detect(rownames(metadata), '^SMD'), 'MDS', 'Elder')
metadata$replicate <- stringr::str_extract(rownames(metadata),'(^[A-Z0-9]+)(?=_|-)')

metadata$label <- factor(metadata$label, levels = c('MDS', 'Elder'))

metadata <- metadata[rownames(metadata) %in% colnames(expression), ]

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_flt <- DE[DE$p_val_adj < 0.05,]

all_results_psuedo <- split( DE_flt , f = DE_flt$cell_type )
saveRDS(all_results_psuedo ,'./Data/Results_DE_Non5qVsElder.rds')



pdf('./Plots/Volcano_And_tornado_Pseudo_Non5qVsElder.pdf')
for (cell_type in sort(names(all_results_psuedo))){
	print(cell_type)
	tmp <- all_results_psuedo[[cell_type]]
	tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
	rownames(tt) <- tt$gene_name
	print(get_volcano(tt, cell_type, 'Non5q Vs Elder'))
	tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_logFC, tmp$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
	print(get_tornado(results_c2, cell_type))
	}
}
dev.off()

PLOT_PATH <- paste0(getwd(), '/Plots/')
WriteXLS::WriteXLS(all_results_psuedo, ExcelFileName=paste0('./Plots/PseudoBulk_Results_elderVsNon5q.xlsx'), SheetNames = names(all_results_psuedo),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


matrices = to_pseudobulk(as.matrix(expression), meta = metadata)

saveRDS(matrices, './Data/matrices_pseudobulk_Non5qVsElder.rds')





################################################
# 5q Vs Non5q
################################################

data_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

all_5q_selected_cells <- readRDS('./Data/CopyKat/all_5q_selected_cells.rds') # <- 
cells_selected_CASPER <- readRDS( './Data/cells_selected_CASPER.rds')

real_5q_cells <- intersect(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id)
cells_to_discard <- c(setdiff(all_5q_selected_cells$Cell_id, cells_selected_CASPER$Cell_id), 
					  setdiff(cells_selected_CASPER$Cell_id, all_5q_selected_cells$Cell_id))

# NONE OF THE CELLS TO DISCARD ARE IN THE SEURAT OBJ.

expression <- as.data.frame(data_5q@assays$RNA@counts)
# rownames(expression) <- gsub('-', '\\.', rownames(expression))



# merged_seurat <- merge(elder_integrated_subseted, sc_data_MDS5q_subseted)

metadata <- setNames(data_5q[[c('cell_5q', 'Sample', 'Cluster_names')]], 
			c('label', 'replicate', 'cell_type'))

library(Libra)

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_flt <- DE[DE$p_val_adj < 0.05,]

all_results_psuedo <- split( DE_flt , f = DE_flt$cell_type )

pdf('./Plots/Volcano_And_tornado_Pseudo_5qVsNon5q.pdf')
for (cell_type in sort(names(all_results_psuedo))){
	print(cell_type)
	tmp <- all_results_psuedo[[cell_type]]
	tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
	rownames(tt) <- tt$gene_name
	print(get_volcano(tt, cell_type, 'Non5q Vs Elder'))
	tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_logFC, tmp$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
	print(get_tornado(results_c2, cell_type))
	}
}
dev.off()

PLOT_PATH <- paste0(getwd(), '/Plots/')
WriteXLS::WriteXLS(all_results_psuedo, ExcelFileName=paste0('./Plots/PseudoBulk_Results_5qVsNon5q.xlsx'), 
SheetNames = names(all_results_psuedo),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


matrices = to_pseudobulk(as.matrix(expression), meta = metadata)

saveRDS(matrices, './Data/matrices_pseudobulk_Del5qVsNon5q.rds')



################################################
# 5q Vs Non 5q no cell type
################################################


Samples_5q <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')

DefaultAssay(Samples_5q) <- 'RNA'
Idents(Samples_5q) <- 'cell_5q'
Samples_5q[["RNA"]]@counts<-as.matrix(Samples_5q[["RNA"]]@counts)+1
de_all <- FindMarkers(Samples_5q, slot='counts', test.use = 'DESeq2', ident.1 = "del5q", ident.2 = "normal", verbose = FALSE)
saveRDS(de_all, paste0(getwd(), '/Plots/DE_all_del5qVsnon5q.rds'))
de_all <- na.omit(de_all)
results <- list()
results[['Raw']]      <- de_all
results[['Filtered']] <- de_all[de_all$p_val_adj < 0.05 & abs(de_all$avg_log2FC) > 0.2 ,]

WriteXLS::WriteXLS(results, ExcelFileName=paste0(PLOT_PATH,'/','Results_5qVsNon5q.xlsx'), SheetNames = names(results),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)


Samples_PrePost <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
Samples_Elder   <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')


pdf('Test.pdf')
VlnPlot(Samples_5q, 'HBB')
dev.off()



################################################
# 5q Vs Non 5q by cell type
################################################


sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
# sc_data_MDS5q_subseted <- filter_and_annotate(sc_data_MDS5q_subseted, annotation_lists[['5qSamples']])

cell_type <- sort(unique(sc_data_MDS5q_subseted$Cluster_names))[1]
DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
DE_per_CT <- list()
for (cell_type in sort(unique(sc_data_MDS5q_subseted$Cluster_names))){
	message(cell_type)
	Idents(sc_data_MDS5q_subseted) <- 'Cluster_names'
	tmp <- subset(sc_data_MDS5q_subseted,  idents = cell_type)
	Idents(tmp) <- 'cell_5q'
	tmp[["RNA"]]@counts <- as.matrix(tmp[["RNA"]]@counts)+1
	DE_per_CT[[cell_type]] <- na.omit(FindMarkers(tmp, slot='counts', test.use = 'DESeq2', ident.1 = "del5q", ident.2 = "normal", verbose =FALSE, max.cells.per.ident = 5000))
}

PLOT_PATH <- paste0(getwd(), '/Plots')

saveRDS(DE_per_CT, paste0(PLOT_PATH, '/DE_all_del5qVsnon5qByCelltype.rds'))

WriteXLS::WriteXLS(DE_per_CT, ExcelFileName=paste0('./Plots/DE_all_del5qVsnon5qByCelltypeALL.xlsx'), SheetNames = names(DE_per_CT),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

for (cell_type in names(DE_per_CT)){
	tmp <- DE_per_CT[[cell_type]]
	tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_log2FC) > 0.2,]
	DE_per_CT[[cell_type]] <- tmp
}
WriteXLS::WriteXLS(DE_per_CT, ExcelFileName=paste0('./Plots/DE_all_del5qVsnon5qByCelltype_filtered.xlsx'), SheetNames = names(DE_per_CT),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)





####### Plot genes for Nerea #######
library(Seurat)
library(ggplot2)


genes_2_plot <- c('FIGN', 'ID1', 'CD74', 'RPS14', 'EIF1AY',  'CLSPN', 'MYBL2', 'FBLN1', 'HES6', 'CD70', 'APOE')
genes_2_plot <- c('CSNK1A1', 'TP53', 'RUNX1', 'DDX41')

# sc_data_MDS5q_subseted <- readRDS('./Data/all_seurat_integrated_sct_subcluster_ClusterNames_5qNotation.rds')
elder_integrated_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS_std/all_seurat_Elder_integrated_sct_subcluster.rds')
elder_ann <- setNames(as.data.frame(elder_integrated_subseted$Cluster_names), 'cell_type')
data_5q_ann <- setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type')

DefaultAssay(sc_data_MDS5q_subseted) <- 'SCT'
data_5q <- as.data.frame(sc_data_MDS5q_subseted@assays$SCT@data)
data_5q <- data_5q[rownames(data_5q) %in% genes_2_plot,]
data_Non5q <- data_5q[, colnames(data_5q) %in% names(which(sc_data_MDS5q_subseted$cell_5q == 'normal'))]
data_5q    <- data_5q[, colnames(data_5q) %in% names(which(sc_data_MDS5q_subseted$cell_5q == 'del5q'))]

elder <- as.data.frame(elder_integrated_subseted@assays$SCT@data)
elder <- elder[rownames(elder) %in% genes_2_plot,]


elder$gene_name <- rownames(elder) 
elder <- reshape2::melt(elder)
elder <- merge(elder, elder_ann, by.x='variable', by.y=0)

data_5q$gene_name <- rownames(data_5q) 
data_5q <- reshape2::melt(data_5q)
data_5q <- merge(data_5q, data_5q_ann, by.x='variable', by.y=0)

data_Non5q$gene_name <- rownames(data_Non5q) 
data_Non5q <- reshape2::melt(data_Non5q)
data_Non5q <- merge(data_Non5q, data_5q_ann, by.x='variable', by.y=0)

elder$phenotype <- 'Elder'
data_5q$phenotype <- '5q_cells'
data_Non5q$phenotype <- 'non5q_cells'

all_data <- rbind(rbind(data_5q, data_Non5q), elder)
all_data$phenotype <- factor(all_data$phenotype, levels=c('Elder', 'non5q_cells', '5q_cells'))
all_data$cell_type <- factor(all_data$cell_type, levels=c(sort(unique(all_data$cell_type))))

all_data$sample <- stringr::str_extract(all_data$variable, '^[A-Z0-9]+')

pdf('./Plots/CheckGeneExpression3.pdf', 30, 30)
# pdf('./Plots/CheckGeneExpression.pdf', 20, 20)
ggplot(all_data, aes(x=gene_name, y=value, fill=phenotype)) + geom_boxplot() + 
	# scale_fill_manual(values=c(ggthemes::tableau_color_pal('Classic 20')(20)))+
	scale_fill_brewer(palette = "Dark2") + facet_wrap(~cell_type, scales='free_x') +
	theme_classic() + 
	theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom')

ggplot(all_data, aes(x=gene_name, y=value, fill=phenotype)) + geom_boxplot() + 
	# scale_fill_manual(values=c(ggthemes::tableau_color_pal('Classic 20')(20)))+
	scale_fill_brewer(palette = "Dark2") + facet_wrap(sample~cell_type, scales='free_x') +
	theme_classic() + 
	theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom')
dev.off()





################################################
# Pre Vs Post by cell type
################################################


library(Seurat)
library(future)
plan("multicore", workers = 64)


Samples_PrePost <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
Samples_PrePost$Case <- ifelse(Samples_PrePost$Sample == 'SMD211420', 'Pre', 'Post')


DefaultAssay(Samples_PrePost) <- 'RNA'
DE_per_CT <- list()
for (cell_type in sort(unique(Samples_PrePost$Cluster_names))){
	message(cell_type)
	Idents(Samples_PrePost) <- 'Cluster_names'
	tmp <- subset(Samples_PrePost,  idents = cell_type)
	Idents(tmp) <- 'Case'
	tmp[["RNA"]]@counts <- as.matrix(tmp[["RNA"]]@counts)+1
	DE_per_CT[[cell_type]] <- na.omit(FindMarkers(tmp, slot='counts', test.use = 'DESeq2', ident.1 = "Pre", ident.2 = "Post", verbose =FALSE, max.cells.per.ident = 5000))
}

PLOT_PATH <- paste0(getwd(), '/Plots')

saveRDS(DE_per_CT, paste0(PLOT_PATH, '/DE_all_PrePostByCelltype.rds'))

WriteXLS::WriteXLS(DE_per_CT, ExcelFileName=paste0('./Plots/DE_all_PrePostByCelltypeALL.xlsx'), SheetNames = names(DE_per_CT),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)

for (cell_type in names(DE_per_CT)){
	tmp <- DE_per_CT[[cell_type]]
	tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_log2FC) > 0.2,]
	DE_per_CT[[cell_type]] <- tmp
}
WriteXLS::WriteXLS(DE_per_CT, ExcelFileName=paste0('./Plots/DE_all_PrePostByCelltype_filtered.xlsx'), SheetNames = names(DE_per_CT),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)






################################################
# 5q Vs non5q with Bulk
################################################
library(Libra)

sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
# LIBRA uses the raw counts

expression <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
metadata <-	 setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type')
metadata$label <- sc_data_MDS5q_subseted$cell_5q
metadata$replicate <- sc_data_MDS5q_subseted$Sample

metadata$label <- factor(metadata$label, levels = c('del5q', 'normal'))
metadata <- metadata[rownames(metadata) %in% colnames(expression), ]

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE_plotter <- DE[DE$gene %in% c('PRSS21', 'ZCCHC10', 'FIGN', 'SIL1',  'CCL5', 'MAP3K7CL',  'IGKC'), c('cell_type', 'gene', 'avg_logFC', 'p_val_adj')]



DE_plotter_avg_logFC <- reshape2::dcast(DE_plotter[, c('cell_type', 'gene', 'avg_logFC')], cell_type ~ gene)
DE_plotter_p_val_adj <- reshape2::dcast(DE_plotter[, c('cell_type', 'gene', 'p_val_adj')], cell_type ~ gene)
rownames(DE_plotter_p_val_adj) <- DE_plotter_p_val_adj$cell_type
rownames(DE_plotter_avg_logFC) <- DE_plotter_avg_logFC$cell_type
DE_plotter_p_val_adj <- DE_plotter_p_val_adj[, -1]
DE_plotter_avg_logFC <- DE_plotter_avg_logFC[, -1]

quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n), na.rm =TRUE)
  breaks[!duplicated(breaks, na.rm=TRUE)]
}
mat_breaks <- quantile_breaks(-log(DE_plotter_p_val_adj), n = 11)

pdf('./Plots/DE_5qVsNon5q.pdf')
log_Pval <- -log(DE_plotter_p_val_adj)
log_Pval[log_Pval>10] <- 10
print(pheatmap::pheatmap(t(log_Pval),cellwidth = 20, cellheight = 20,
cluster_rows=FALSE, cluster_cols=FALSE, main='-log10(p_val_adj)', 
color = colorRampPalette(c("white", "yellow", "red"))(50)))

print(pheatmap::pheatmap(t(DE_plotter_avg_logFC),cellwidth = 20, cellheight = 20,
cluster_rows=FALSE, cluster_cols=FALSE, main='avg_logFC', 
color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="Spectral")))(100)))
dev.off()

DE_flt <- DE[DE$p_val_adj < 0.1,]
all_results_psuedo <- split( DE_flt , f = DE_flt$cell_type )
saveRDS(all_results_psuedo ,'./Data/Results_DE_5qVsNon5q.rds')

pdf('./Plots/Volcano_And_tornado_Pseudo_5qVsNon5q.pdf')
for (cell_type in sort(names(all_results_psuedo))){
	print(cell_type)
	tmp <- all_results_psuedo[[cell_type]]
	tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
	rownames(tt) <- tt$gene_name
	print(get_volcano(tt, cell_type))
	tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
	ranks <- setNames(tmp$avg_logFC, tmp$gene)
	results_c2 <- fgsea::fgsea(pathways_c2, ranks, minSize=15, maxSize = 600)
	results_c2 <- results_c2[results_c2$padj < 0.05,]
	results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
	results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
	results_c2$pathway <- factor(results_c2$pathway)
	if(nrow(results_c2)>0){
	print(get_tornado(results_c2, cell_type))
	}
}
dev.off()


WriteXLS::WriteXLS(all_results_psuedo, ExcelFileName=paste0('./Plots/PseudoBulk_Results_5qVsNon5q.xlsx'), SheetNames = names(all_results_psuedo),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)

matrices = to_pseudobulk(as.matrix(expression), meta = metadata)

saveRDS(matrices, './Data/matrices_pseudobulk_5qVsNon5q.rds')
matrices <- readRDS('./Data/matrices_pseudobulk_5qVsNon5q.rds')
genes_2_plot <- c('PRSS21', 'ZCCHC10', 'FIGN', 'SIL1',  'CCL5', 'MAP3K7CL',  'IGKC')

all_matrices <- NULL
for (mat in names(matrices)) {
	tmp <- matrices[[mat]]
	# tmp <- tmp[rownames(tmp) %in% genes_2_plot,]
	groups <- stringr::str_extract(colnames(tmp), '(?<=:)[\\w]+')
	colnames(tmp) <- paste0(colnames(tmp), '_', mat)
	dds <- edgeR::DGEList(counts = tmp, group = groups)
	dds <- edgeR::calcNormFactors(dds, method = 'TMM')
	norm_counts <- edgeR::cpm(dds,normalized=T)
	if (is.null(all_matrices)) {
		all_matrices <- norm_counts
	} else {
		all_matrices <- merge(all_matrices, norm_counts, by=0, all=TRUE)
		rownames(all_matrices) <- all_matrices$Row.names
		all_matrices <- all_matrices[,-1]
	}
}

genes_2_plot %in% rownames(all_matrices)

tmp <- all_matrices[rownames(all_matrices) %in% genes_2_plot,]

pdf('./Plots/Test.pdf')
cell_type <- stringr::str_extract(colnames(tmp), '(?<=_)[\\w-]+')
genotype <- stringr::str_extract(colnames(tmp), '(?<=:)[a-z0-9]+')
annotation_columns <-  data.frame(cell_type = cell_type, genotype = genotype)
rownames(annotation_columns) <- colnames(tmp)
genotype_colors <- c('black', 'red')
names(genotype_colors) <- unique(genotype)
cell_type_colors <- c(ggthemes::tableau_color_pal('Classic 20')(14))
names(cell_type_colors) <- unique(cell_type)
anno_colors <- list(genotype=genotype_colors, cell_type=cell_type_colors)
pheatmap::pheatmap(tmp, scale='row', annotation = annotation_columns, cluster_cols=F, annotation_colors = anno_colors)
dev.off()


################################################################################
# ADDITIONAL PLOTTING // INFORMATIVE PLOTS
################################################################################


Samples_5q      <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
Samples_PrePost <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
Samples_Elder   <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')

counts <- cbind(as.data.frame(GetAssayData(object = Samples_5q, slot = "scale.data")['CD69' , , drop=F]),
				as.data.frame(GetAssayData(object = Samples_Elder, slot = "scale.data")['CD69' , , drop=F]))

counts <- as.data.frame(t(counts))
counts$cell_id <- rownames(counts)

metadata <- Samples_Elder[['Cluster_names']]
metadata$cell_5q <- 'Elder'
metadata <-  rbind(metadata, Samples_5q[[c('Cluster_names', 'cell_5q')]])
counts <- merge(counts, metadata, by.x='cell_id', by.y=0)

pdf('./Plots/Expression_CD69.pdf', 10, 10)
ggplot(counts, aes(x=Cluster_names, y=CD69, fill=cell_5q)) + geom_boxplot() + 
	scale_fill_manual(values=c(Elder='#708d81', del5q='#8d0801', normal='#f6aa1c'))+
	theme_classic() + 
	theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom')

ggplot(counts[counts$cell_5q != 'Elder',], aes(x=cell_5q, y=CD69, fill=cell_5q)) + geom_boxplot() + 
	scale_fill_manual(values=c(del5q='#8d0801', normal='#f6aa1c'))+
	ggsignif::geom_signif(comparisons = list(c("del5q", "normal")), map_signif_level = TRUE, textsize=5) +
	theme_classic() + facet_wrap(~Cluster_names) + ylim(-7, 15)+
	theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom')
dev.off()




quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(as.matrix(xs), probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

library(viridis)
library(RColorBrewer)



Results_DE_Non5qVsElder          <- readRDS('./Data/Results_DE_Non5qVsElder.rds')
matrices_pseudobulk_Non5qVsElder <- readRDS('./Data/matrices_pseudobulk_Non5qVsElder.rds')

pdf('./Plots/Heatmaps_5qVsElder.pdf', 10, 10)
for (contrast in names(Results_DE_Non5qVsElder)){
	tmp <- Results_DE_Non5qVsElder[[contrast]]
	tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_logFC) > 2,]
	data_tmp <- matrices_pseudobulk_Non5qVsElder[[contrast]]
	data_tmp <- data_tmp[tmp$gene,]
	data_tmp <- data_tmp[,order(colnames(data_tmp))]
	data_tmp <- log(data_tmp+1)
	mat_breaks <- quantile_breaks(data_tmp, n = 11)
	annotation_col <- data.frame(Phenotype = stringr::str_extract(colnames(data_tmp), '(?<=:)[\\w]+'))
	rownames(annotation_col) <- colnames(data_tmp)
	anno_colors <- list()
	anno_colors[['Phenotype']] <- c(Elder='#708d81', MDS='#8d0801')
	print(pheatmap::pheatmap(data_tmp, color=rev(RColorBrewer::brewer.pal(11, 'Spectral')) ,
							 annotation_col = annotation_col, annotation_colors = anno_colors,
							 breaks = mat_breaks, fontsize_row =3, main=contrast))
}
dev.off()



Results_DE_5qVsElder          <- readRDS('./Data/Results_DE_5qVsElder.rds')
matrices_pseudobulk_5qVsElder <- readRDS('./Data/matrices_pseudobulk_5qVsElder.rds')
pdf('./Plots/Heatmaps_Non5qVsElder.pdf', 10, 10)
for (contrast in intersect(names(Results_DE_5qVsElder), names(matrices_pseudobulk_5qVsElder))){
	tmp <- Results_DE_5qVsElder[[contrast]]
	tmp <- tmp[tmp$p_val_adj < 0.05 & abs(tmp$avg_logFC) > 2,]
	data_tmp <- matrices_pseudobulk_5qVsElder[[contrast]]
	data_tmp <- data_tmp[tmp$gene,]
	data_tmp <- data_tmp[,order(colnames(data_tmp))]
	data_tmp <- na.omit(log(data_tmp+1))
	mat_breaks <- quantile_breaks(data_tmp, n = 11)
	annotation_col <- data.frame(Phenotype = stringr::str_extract(colnames(data_tmp), '(?<=:)[\\w]+'))
	rownames(annotation_col) <- colnames(data_tmp)
	anno_colors <- list()
	anno_colors[['Phenotype']] <- c(Elder='#708d81', MDS_5q='#8d0801')
	print(pheatmap::pheatmap(data_tmp, color=rev(RColorBrewer::brewer.pal(11, 'Spectral')) ,
							 annotation_col = annotation_col, annotation_colors = anno_colors,
							 breaks = mat_breaks, fontsize_row =3, main=contrast))
}
dev.off()







#### Comparative abundancies
mds_5q_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
elder_data  <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')


abundances <- rbind(mds_5q_data[[c('Cluster_names', 'Sample')]], elder_data[[c('Cluster_names', 'Sample')]])
abundances <- table(abundances$Cluster_names, abundances$Sample)
abundances <- unclass(abundances) 

library(edgeR)
# Attaching some column metadata.
extra.info <- colData(merged)[match(colnames(abundances), merged$sample),]
y.ab <- DGEList(abundances)
y.ab[['samples']]['batch']	<- seq(1, ncol(abundances))
y.ab[['samples']]['MDS']	<- ifelse(grepl('^SMD', rownames(y.ab[['samples']])), 'SMD', 'Elder')
design <- model.matrix(~factor(MDS), y.ab$samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)

fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)


pdf('./Plots/Abundances.pdf', 10, 10)
plotBCV(y.ab, cex=1)
plotQLDisp(fit.ab, cex=1)
dev.off()


res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)





abundances <- rbind(mds_5q_data[[c('Cluster_names', 'Sample')]], elder_data[[c('Cluster_names', 'Sample')]])
abundances <- table(abundances$Cluster_names, abundances$Sample)

abundances<-prop.table(abundances, margin = 2)

RES <- data.frame(CellType=NULL, p.value=NULL)
for (i in 1:nrow(abundances)){
	tmp <- abundances[i, ]
	res <- wilcox.test(tmp[grepl('GSM', names(tmp))], tmp[grepl('SMD', names(tmp))])
	var_LEDR <- var(tmp[grepl('GSM', names(tmp))])
	var_MDS <- var(tmp[grepl('SMD', names(tmp))])
	RES <- rbind(RES, data.frame(CellType=rownames(abundances)[i], p.value=res$p.value))
}







# Young senior and MDS
df <- reshape2::melt(abundances)
colnames(df)<-c("CellType", "Patient", "Proportion")
df$Condition <- ifelse(grepl('^SMD', df$Patient), 'MDS', 'Elder')

##### Group per condition #####
#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df2<-data_summary(df, varname = "Proportion", groupnames = c("Condition", "CellType"))
df2$sd[is.na(df2$sd)]<-0

df2 <- data.frame(Condition=NULL, CellType=NULL, Proportion=NULL, sd=NULL)
for (CT in unique(df$CellType)){
	tmp <- abundances[CT, ]
	df2 <- rbind(df2, data.frame(Condition= 'Elder', CellType = CT, Proportion = mean(tmp[grepl('GSM', names(tmp))]), sd = sd(tmp[grepl('GSM', names(tmp))])))
	df2 <- rbind(df2, data.frame(Condition= 'MDS'  , CellType = CT, Proportion = mean(tmp[grepl('SMD', names(tmp))]), sd = sd(tmp[grepl('SMD', names(tmp))])))
}

pdf("./Plots/Proportions_sd.pdf", width = 10, useDingbats = F)
ggplot(df2, aes(CellType, Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_errorbar(aes(ymin = Proportion - sd, ymax = Proportion + sd), width=.2, position=position_dodge(.9)) +
  theme_classic() + 
  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') +
  theme(axis.text=element_text(family = "Helvetica", size=8),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=8),
			legend.title=element_text(size=7),
			panel.grid.major.y = element_line( size=.1, color="grey" ),
			legend.position = 'top')
dev.off()


pdf("./Plots/Proportions_lines.pdf", width = 10, useDingbats = F)
ggplot(df, aes(x=CellType, y=Proportion*100, color = Condition)) +
 geom_point(aes(shape=Patient)) + geom_line(aes(group=Patient)) +
  theme_classic()  +  scale_color_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') +
  theme(axis.text=element_text(family = "Helvetica", size=8),
			axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
			axis.title.x = element_blank(),
			axis.title = element_text(size=7, family = "Helvetica"),
			axis.line = element_line(colour = 'black', size = 0.5),
			legend.text=element_text(size=8),
			legend.title=element_text(size=7),
			panel.grid.major.y = element_line( size=.1, color="grey" ),
			legend.position = 'top')
dev.off()




sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
tmp <- sc_data_MDS5q_subseted[[c('Cluster_names', 'cell_5q', 'S.Score', 'G2M.Score', 'Phase')]]

pdf('./Plots/Proliferation_ratios_5qVsnormal.pdf')
# ggplot(tmp, aes(x=S.Score, y= G2M.Score, color=cell_5q)) + geom_point(alpha=0.8) + 
# facet_wrap(~Cluster_names) +  scale_color_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
# theme_classic() + 
# theme(axis.text=element_text(family = "Helvetica", size=8),
# 		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
# 		axis.title.x = element_blank(),
# 		axis.title = element_text(size=7, family = "Helvetica"),
# 		axis.line = element_line(colour = 'black', size = 0.5),
# 		legend.text=element_text(size=8),
# 		legend.title=element_text(size=7),
# 		panel.grid.major.y = element_line( size=.1, color="grey" ),
# 		legend.position = 'top')

cowplot::plot_grid(
ggplot(tmp, aes(x=Cluster_names, y= S.Score, fill=cell_5q)) + geom_boxplot() + 
theme_classic() +  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top')
,
ggplot(tmp, aes(x=Cluster_names, y= G2M.Score, fill=cell_5q)) + geom_boxplot() + 
theme_classic() +  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top'),
nrow=2)

ggplot(tmp, aes(x=cell_5q, y= G2M.Score, fill=cell_5q)) + geom_boxplot() + 
facet_wrap(~Cluster_names) + ylim(round(min(tmp$G2M.Score), 1),2)+
ggsignif::geom_signif(comparisons = list(c("del5q", "normal")), y_position = 1.7,map_signif_level = TRUE) +
theme_classic() +  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top', strip.background = element_blank(), strip.placement = "outside")

ggplot(tmp, aes(x=cell_5q, y= S.Score, fill=cell_5q)) + geom_boxplot() + 
facet_wrap(~Cluster_names) + ylim(round(min(tmp$S.Score), 1),2)+
ggsignif::geom_signif(comparisons = list(c("del5q", "normal")), y_position = 1.7,map_signif_level = TRUE) +
theme_classic() +  scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top', strip.background = element_blank(), strip.placement = "outside")

dev.off()


sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
DefaultAssay(sc_data_MDS5q_subseted) <- 'SCT'
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
genes_2_check <- c('CD74', 'RPS14', 'BTF3', 'COX7C', 'HINT1', 'RPS23')

pdf('./Plots/CD74_RPS14.pdf')
VlnPlot(object = sc_data_MDS5q_subseted, features = genes_2_check, pt.size = 0,
split.by = 'cell_5q', group.by='Cluster_names')
dev.off()

plotter <- as.data.frame(sc_data_MDS5q_subseted@assays$SCT@data)
plotter <- t(plotter[genes_2_check,])
plotter <- merge(plotter, sc_data_MDS5q_subseted[[c('cell_5q', 'Cluster_names')]], by=0)

pdf('./Plots/CD74_RPS14.pdf', width=5, height=7)

# ggplot(reshape2::melt(plotter), aes(x=cell_5q, y=value , fill=cell_5q)) + 
# geom_boxplot() + scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
# theme_classic() + facet_wrap(variable~Cluster_names) + ylim(-0.1, 5)+
# ggsignif::geom_signif(comparisons = list(c("del5q", "normal")), y_position = 4,map_signif_level = TRUE) +
# theme(axis.text=element_text(family = "Helvetica", size=8),
# 		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
# 		axis.title.x = element_blank(),
# 		axis.title = element_text(size=7, family = "Helvetica"),
# 		axis.line = element_line(colour = 'black', size = 0.5),
# 		legend.text=element_text(size=8),
# 		legend.title=element_text(size=7),
# 		panel.grid.major.y = element_line( size=.1, color="grey" ),
# 		legend.position = 'top',  strip.placement = "outside")

ggplot(reshape2::melt(plotter), aes(x=cell_5q, y=value , fill=cell_5q)) + 
geom_boxplot() + scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme_classic() + facet_wrap(~variable) + ylim(-0.1, 7)+
ggsignif::geom_signif(comparisons = list(c("del5q", "normal")), y_position = 6,map_signif_level = TRUE) +
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top', strip.background = element_blank(), strip.placement = "outside")

ggplot(reshape2::melt(plotter), aes(x=Cluster_names, y=value , fill=cell_5q)) + 
geom_boxplot() + scale_fill_manual(values=c('#bb3e03', '#0a9396'), name='Genotype') + 
theme_classic() + facet_wrap(~variable, ncol=1) +
theme(axis.text=element_text(family = "Helvetica", size=8),
		axis.text.x = element_text( angle=90, vjust=0.5, hjust=1), 
		axis.title.x = element_blank(),
		axis.title = element_text(size=7, family = "Helvetica"),
		axis.line = element_line(colour = 'black', size = 0.5),
		legend.text=element_text(size=8),
		legend.title=element_text(size=7),
		panel.grid.major.y = element_line( size=.1, color="grey" ),
		legend.position = 'top', strip.background = element_blank(), strip.placement = "outside")
dev.off()







get_venn <- function(list_2_plot, title)
{
	venn <- ggvenn::ggvenn(list_2_plot,
		fill_color = destiny::cube_helix(length(unique(data_5q$Sample))),
		stroke_size = 0.4,
		show_percentage = TRUE,
		fill_alpha = 0.4,
		stroke_color = 'white',
		stroke_alpha = 1,
		stroke_linetype = 'solid',
		text_color = 'black',
		set_name_size = 4,
		text_size = 2)+ggtitle(title)
	return(venn)
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


all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')
all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')

all_results_psuedo_5qVsElder <- do.call("rbind", all_results_psuedo_5qVsElder)
all_results_psuedo_Non5qVsElder <- do.call("rbind", all_results_psuedo_Non5qVsElder)



Pos_5qVsElder    <- all_results_psuedo_5qVsElder[all_results_psuedo_5qVsElder$p_val_adj < 0.05 & all_results_psuedo_5qVsElder$avg_logFC >0, c('cell_type', 'gene')]
Pos_non5qVsElder <- all_results_psuedo_Non5qVsElder[all_results_psuedo_Non5qVsElder$p_val_adj < 0.05 & all_results_psuedo_Non5qVsElder$avg_logFC >0, c('cell_type', 'gene')]
Neg_5qVsElder    <- all_results_psuedo_5qVsElder[all_results_psuedo_5qVsElder$p_val_adj < 0.05 & all_results_psuedo_5qVsElder$avg_logFC <0, c('cell_type', 'gene')]
Neg_non5qVsElder <- all_results_psuedo_Non5qVsElder[all_results_psuedo_Non5qVsElder$p_val_adj < 0.05 & all_results_psuedo_Non5qVsElder$avg_logFC <0, c('cell_type', 'gene')]

Pos_5qVsElder <- split(Pos_5qVsElder, Pos_5qVsElder$cell_type)
Pos_5qVsElder <- lapply(Pos_5qVsElder, function(x) x$gene)

Pos_non5qVsElder <- split(Pos_non5qVsElder, Pos_non5qVsElder$cell_type)
Pos_non5qVsElder <- lapply(Pos_non5qVsElder, function(x) x$gene)

Neg_5qVsElder <- split(Neg_5qVsElder, Neg_5qVsElder$cell_type)
Neg_5qVsElder <- lapply(Neg_5qVsElder, function(x) x$gene)

Neg_non5qVsElder <- split(Neg_non5qVsElder, Neg_non5qVsElder$cell_type)
Neg_non5qVsElder <- lapply(Neg_non5qVsElder, function(x) x$gene)


pdf('./Plots/DE_Overlap_DE_Pos_5qVsElder.pdf')
	get_upset(Pos_5qVsElder)
dev.off()

pdf('./Plots/DE_Overlap_DE_Pos_non5qVsElder.pdf')
	get_upset(Pos_non5qVsElder)	
dev.off()

pdf('./Plots/DE_Overlap_DE_Neg_5qVsElder.pdf')
	get_upset(Neg_5qVsElder)
dev.off()

pdf('./Plots/DE_Overlap_DE_Neg_non5qVsElder.pdf')
	get_upset(Neg_non5qVsElder)
dev.off()






##################################################################################################
#     Tests for enrichment in DEG experiments using enrichR and fgsea 
##################################################################################################
library(fgsea)
library(ggplot2)
library(enrichR)

all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')

terms_2_plot <- c("DNA replication (GO:0006260)", "DNA metabolic process (GO:0006259)", "regulation of transcription involved in G1/S transition of mitotic cell cycle (GO:0000083)", "DNA replication initiation (GO:0006270", "mitotic sister chromatid segregation (GO:0000070)", "mitotic DNA replication (GO:1902969)", "DNA replication checkpoint signaling (GO:0000076)", "mitotic chromosome condensation (GO:0007076)", "mitotic nuclear division (GO:0140014)", "DNA-dependent DNA replication (GO:0006261)", "DNA strand elongation involved in DNA replication (GO:0006271)","peptide biosynthetic process (GO:0043043)", "cytoplasmic translation (GO:0002181)", "translation (GO:0006412)", "cellular protein metabolic process (GO:0044267)", "ribosome biogenesis (GO:0042254)", "rRNA processing (GO:0006364)", "mitochondrial translational elongation (GO:0070125)", "mitochondrial translational termination (GO:0070126)", "translational termination (GO:0006415)", "mitochondrial translation (GO:0032543)", "positive regulation of apoptotic process (GO:0043065)","intrinsic apoptotic signaling pathway (GO:0097193)", 'p53 signaling pathway','DNA repair (GO:0006281)')
# all_results_psuedo_5qVsElder <- do.call("rbind", all_results_psuedo_5qVsElder)
# all_results_psuedo_Non5qVsElder <- do.call("rbind", all_results_psuedo_Non5qVsElder)


setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human')

all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
for (CT in names(all_results_psuedo_5qVsElder)){
	message(CT)
	tmp <- all_results_psuedo_5qVsElder[[CT]]
	tmp_pos <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >0,'gene']
	tmp_neg <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC <0,'gene']
	enriched <- enrichr(tmp_pos, dbs)
	enriched <- do.call("rbind", enriched)
	enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
	if(nrow(enriched) > 1){
		all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
	}
	enriched <- enrichr(tmp_neg, dbs)
	enriched <- do.call("rbind", enriched)
	enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
	if(nrow(enriched) > 1){
		all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))
	}
}



pdf('./Plots/DotPlot_5qVsElder_enrichR.pdf')
tmp <- all_results
tmp[tmp$Combined.Score >200, 'Combined.Score'] <- 200
tmp$Combined.Score_size <- tmp$Combined.Score
tmp <- tmp[tmp$Term %in% terms_2_plot, ]
tmp[tmp$is_pos == FALSE, 'Combined.Score'] <- tmp[tmp$is_pos == FALSE, 'Combined.Score']*-1
ggplot() + 
geom_point(tmp, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
# scale_fill_distiller(palette = "Spectral", direction=-1, name='Score') +
scale_fill_gradient2(low =scales::muted("blue"),mid = "white",high = scales::muted("red"),midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(6, 1)) +
theme(axis.text.y=element_text(size=6),axis.text.x=element_text(angle=45, , hjust=1)) 
dev.off()




all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')

all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
for (CT in names(all_results_psuedo_Non5qVsElder)){
	message(CT)
	tmp <- all_results_psuedo_Non5qVsElder[[CT]]
	tmp_pos <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC >0,'gene']
	tmp_neg <- tmp[tmp$p_val_adj < 0.05 & tmp$avg_logFC <0,'gene']
	enriched <- enrichr(tmp_pos, dbs)
	enriched <- do.call("rbind", enriched)
	enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
	if(nrow(enriched) > 1){
		all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=TRUE))
	}
	enriched <- enrichr(tmp_neg, dbs)
	enriched <- do.call("rbind", enriched)
	enriched <- enriched[enriched$Adjusted.P.value < 0.05,]
	if(nrow(enriched) > 1){
		all_results <- rbind(all_results, data.frame(Term=enriched$Term, Adjusted.P.value=enriched$Adjusted.P.value, Combined.Score=enriched$Combined.Score, Genes=enriched$Genes, cell_type=CT, is_pos=FALSE))
	}
}

pdf('./Plots/DotPlot_Non5qVsElder_enrichR.pdf')
tmp <- all_results
tmp[tmp$Combined.Score >200, 'Combined.Score'] <- 200
tmp$Combined.Score_size <- tmp$Combined.Score
tmp <- tmp[tmp$Term %in% terms_2_plot, ]
tmp[tmp$is_pos == FALSE, 'Combined.Score'] <- tmp[tmp$is_pos == FALSE, 'Combined.Score']*-1
ggplot() + 
geom_point(tmp, mapping= aes(y=Term, x=cell_type, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
# scale_fill_distiller(palette = "Spectral", direction=-1, name='Score') +
scale_fill_gradient2(low =scales::muted("blue"),mid = "white",high = scales::muted("red"),midpoint = 0,)+
scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30)) +
scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(6, 1)) +
theme(axis.text.y=element_text(size=6),axis.text.x=element_text(angle=45, , hjust=1)) 
dev.off()




library(fgsea)
pathways_c2 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c2.cp.v7.5.1.symbols.gmt')
pathways_c5 <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/c5.go.v2023.1.Hs.symbols.gmt')
pathways_HM <- fgsea::gmtPathways('/home/tereshkova/data/gserranos/MDS/Data/Annotation/E-GEOD-100618-marker-genes-files/h.all.v2023.1.Hs.symbols.gmt')


get_gsea_results_perCT <- function(data_DE, dataset){
	results <- data.frame(pathway=NULL, pval=NULL, padj=NULL, log2err=NULL, ES=NULL, NES=NULL, size=NULL, leadingEdge=NULL, cell_type=NULL , is_pos=NULL)
	for (CT in names(data_DE)){
		message(CT)
		tmp <- data_DE[[CT]]
		tt <-  data.frame(gene_name = tmp$gene, p_val_adj=as.numeric(tmp$p_val_adj), avg_logFC=as.numeric(tmp$avg_logFC) )
		rownames(tt) <- tt$gene_name
		tmp <- tmp[order(tmp$avg_logFC, decreasing=TRUE),]
		ranks <- setNames(tmp$avg_logFC, tmp$gene)
		results_c2 <- fgsea::fgsea(dataset, ranks, minSize=15, maxSize = 600)
		# results_c2 <- results_c2[results_c2$padj < 0.05,]
		results_c2 <- results_c2[order(results_c2$NES, decreasing=TRUE),]
		results_c2$is_pos <- ifelse(results_c2$NES<=0, 'Neg', 'Pos')
		results_c2$pathway <- factor(results_c2$pathway)
		results_c2$cell_type <- CT
		results <- rbind(results, results_c2)
	}
	return(results)
}

get_gsea_dotplot <- function(results, terms_2_plot=NULL, threshold=3){
	plotter <- as.data.frame(results)
	plotter[plotter$NES >threshold , 'NES'] <- threshold
	plotter[plotter$NES <threshold*-1 , 'NES'] <- threshold*-1
	if(!is.null(terms_2_plot)){
		plotter <- plotter[plotter$pathway %in% terms_2_plot, ]
	}
	ggplot() + 
	geom_point(plotter, mapping= aes(y=pathway, x=cell_type, size =padj, fill=NES), color='black', shape=21) + 
	scale_fill_distiller(palette = "Spectral", direction=-1, name='NES') +
	scale_y_discrete(labels = function(x) stringr::str_wrap(tolower(gsub('_', ' ', x)), width = 30)) +
	scale_alpha(range = c(1, 0.1)) + theme_minimal() +   scale_size_continuous(range = c(6, 1)) +
	theme(axis.text.y=element_text(size=6),axis.text.x=element_text(angle=45, , hjust=1)) 
}

terms_2_plot <- c('GOBP_CELL_CYCLE', 'GOBP_CELL_CYCLE_CHECKPOINT_SIGNALING', 'GOBP_CELL_CYCLE_G2_M_PHASE_TRANSITION', 'GOBP_CELL_CYCLE_PHASE_TRANSITION', 'GOBP_CELL_CYCLE_PROCESS', 'GOBP_DNA_REPLICATION', 'GOBP_REGULATION_OF_CELL_CYCLE_PHASE_TRANSITION', 'WP_RETINOBLASTOMA_GENE_IN_CANCER', 'GOMF_TRANSCRIPTION_COREGULATOR_ACTIVITY', 'GOBP_DOUBLE_STRAND_BREAK_REPAIR', 'REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53', 'GOCC_CYTOSOLIC_RIBOSOME', 'GOCC_LARGE_RIBOSOMAL_SUBUNIT', 'GOCC_RIBOSOMAL_SUBUNIT', 'GOMF_STRUCTURAL_CONSTITUENT_OF_RIBOSOME', 'REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION', 'REACTOME_EUKARYOTIC_TRANSLATION_INITIATION', 'REACTOME_TRANSLATION', 'REACTOME_MITOCHONDRIAL_TRANSLATION')


all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')

results_c2 <- get_gsea_results_perCT(all_results_psuedo_5qVsElder, pathways_c2)
results_c5 <- get_gsea_results_perCT(all_results_psuedo_5qVsElder, pathways_c5)
results_HM <- get_gsea_results_perCT(all_results_psuedo_5qVsElder, pathways_HM)

results_c2 <- results_c2[results_c2$padj < 0.05, ]
results_c5 <- results_c5[results_c5$padj < 0.05, ]
results_HM <- results_HM[results_HM$padj < 0.05, ]

pdf('./Plots/DotPlot_5qVsElder_fgsea_ind.pdf')
get_gsea_dotplot(results_c2, terms_2_plot, 3) + ggtitle('C2')
get_gsea_dotplot(results_c5, terms_2_plot, 3) + ggtitle('C5')
get_gsea_dotplot(results_HM, terms_2_plot, 3) + ggtitle('HM')
dev.off()

pdf('./Plots/DotPlot_5qVsElder_fgsea.pdf')
all_results <- rbind(results_c2, results_c5, results_HM)
get_gsea_dotplot(all_results, terms_2_plot, 3) + ggtitle('C2+C5+HM')
dev.off()



all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')

results_c2 <- get_gsea_results_perCT(all_results_psuedo_Non5qVsElder, pathways_c2)
results_c5 <- get_gsea_results_perCT(all_results_psuedo_Non5qVsElder, pathways_c5)
results_HM <- get_gsea_results_perCT(all_results_psuedo_Non5qVsElder, pathways_HM)

results_c2 <- results_c2[results_c2$padj < 0.05, ]
results_c5 <- results_c5[results_c5$padj < 0.05, ]
results_HM <- results_HM[results_HM$padj < 0.05, ]

pdf('./Plots/DotPlot_Non5qVsElder_fgsea_ind.pdf')
get_gsea_dotplot(results_c2, terms_2_plot, 3) + ggtitle('C2')
get_gsea_dotplot(results_c5, terms_2_plot, 3) + ggtitle('C5')
get_gsea_dotplot(results_HM, terms_2_plot, 3) + ggtitle('HM')
dev.off()

pdf('./Plots/DotPlot_Non5qVsElder_fgsea.pdf')
all_results <- rbind(results_c2, results_c5, results_HM)
get_gsea_dotplot(all_results, terms_2_plot, 3) + ggtitle('C2+C5+HM')
dev.off()





results_4_Nerea <- list(C2_pathways=results_c2, GO_pathways=results_c5, HM_pathways=results_HM)
WriteXLS::WriteXLS(results_4_Nerea, ExcelFileName=paste0('/home/tereshkova/data/gserranos/MDS/Data/DE_results/Data_5qVsElder_fgsea.xlsx'), SheetNames = names(results_4_Nerea),  col.names=TRUE, row.names=FALSE, BoldHeaderRow=TRUE)


all_results_psuedo_5qVsElder    <- readRDS('./Data/Results_DE_5qVsElder.rds')
all_results_psuedo_Non5qVsElder <- readRDS('./Data/Results_DE_Non5qVsElder.rds')

all_results_psuedo_5qVsElder <- do.call("rbind", all_results_psuedo_5qVsElder)
all_results_psuedo_Non5qVsElder <- do.call("rbind", all_results_psuedo_Non5qVsElder)

DEG_5qVsElder <-  all_results_psuedo_5qVsElder[all_results_psuedo_5qVsElder$p_val_adj <0.05,]
DEG_Non5qVsElder <-  all_results_psuedo_Non5qVsElder[all_results_psuedo_Non5qVsElder$p_val_adj <0.05,]

DEG_5qVsElder$is_pos <- ifelse(DEG_5qVsElder$avg_logFC >0, 'UP', 'Neg')
DEG_5qVsElder$Comparison <- 'del(5q) cells Vs Older healthy cells'
DEG_Non5qVsElder$is_pos <- ifelse(DEG_Non5qVsElder$avg_logFC >0, 'UP', 'Neg')
DEG_Non5qVsElder$Comparison <- 'WT cells Vs Older healthy cells'


colors_SP <- RColorBrewer::brewer.pal(9, "Spectral")
plotter <- rbind(DEG_Non5qVsElder, DEG_5qVsElder)
plotter <- as.data.frame(table(plotter$is_pos, plotter$Comparison))
plotter$Freq <- ifelse(plotter$Var1 == 'UP', plotter$Freq, plotter$Freq*-1)

p1 <- ggplot(plotter, aes(x=Var2, y=Freq,  fill=Var1)) + 
geom_hline(yintercept = 0) + 
geom_bar(stat = "identity", position = "stack") + theme_classic() + 
theme(axis.text.x=element_text(angle=45, hjust=1),  axis.title.x=element_blank()) + 
scale_fill_manual(values=c(colors_SP[9],colors_SP[1])) + labs(y='Number of DEG', x = 'Comparisons') + guides(fill=guide_legend(title="FDR<0.05"))


plotter <- rbind(DEG_Non5qVsElder, DEG_5qVsElder)
plotter <- plotter[abs(plotter$avg_logFC) >=2,]
plotter <- as.data.frame(table(plotter$is_pos, plotter$Comparison))
plotter$Freq <- ifelse(plotter$Var1 == 'UP', plotter$Freq, plotter$Freq*-1)
p2 <- ggplot(plotter, aes(x=Var2, y=Freq,  fill=Var1)) + 
geom_hline(yintercept = 0) + 
geom_bar(stat = "identity", position = "stack") + theme_classic() + 
theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x=element_blank()) + 
scale_fill_manual(values=c(colors_SP[9],colors_SP[1])) + labs(y='Number of DEG') + guides(fill=guide_legend(title="FDR<0.05\n|logFC|>=2"))

pdf('./Plots/Supplementary_figure_DEG.pdf')
cowplot::plot_grid(p1, p2, ncol=2)
dev.off()



pdf('./Plots/Supplementary_figure_DEG_CellType.pdf', height=10)
colors_SP <- RColorBrewer::brewer.pal(9, "Spectral")
plotter <- rbind(DEG_Non5qVsElder, DEG_5qVsElder)
plotter <- as.data.frame(table(plotter$is_pos, plotter$Comparison, plotter$cell_type))
plotter$Freq <- ifelse(plotter$Var1 == 'UP', plotter$Freq, plotter$Freq*-1)

ggplot(plotter, aes(x=Var2, y=Freq,  fill=Var1)) + 
geom_hline(yintercept = 0) + 
geom_bar(stat = "identity", position = "stack") + theme_classic() + 
theme(axis.text.x=element_text(angle=90, hjust=1),  axis.title.x=element_blank(), legend.position='bottom') + 
scale_fill_manual(values=c(colors_SP[9],colors_SP[1])) + labs(y='Number of DEG') + 
guides(fill=guide_legend(title="FDR<0.05"))  + scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) + facet_wrap(~Var3, scales='free')


plotter <- rbind(DEG_Non5qVsElder, DEG_5qVsElder)
plotter <- plotter[abs(plotter$avg_logFC) >=2,]
plotter <- as.data.frame(table(plotter$is_pos, plotter$Comparison, plotter$cell_type))
plotter$Freq <- ifelse(plotter$Var1 == 'UP', plotter$Freq, plotter$Freq*-1)
ggplot(plotter, aes(x=Var2, y=Freq,  fill=Var1)) + 
geom_hline(yintercept = 0) + 
geom_bar(stat = "identity", position = "stack") + theme_classic() + 
theme(axis.text.x=element_text(angle=90, hjust=1), axis.title.x=element_blank(), legend.position='bottom') + 
scale_fill_manual(values=c(colors_SP[9],colors_SP[1])) + labs(y='Number of DEG') + 
guides(fill=guide_legend(title="FDR<0.05\n|logFC|>=2"))+ scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 20)) + facet_wrap(~Var3, scales='free')
dev.off()

