library(ggplot2)
library(enrichR)

get_enrichr_results_SimiC <- function(data_set, dbs=dbs){
	all_results <- data.frame(Term=NULL, Adjusted.P.value=NULL, Combined.Score=NULL,Genes=NULL, cell_type=NULL, is_pos =NULL)
	for (driver in sort(unique(as.character(data_set$driver)))){
		message(driver)
		target_genes <- as.character(data_set[data_set$driver == driver & data_set$value!=0,'target'])
		enriched <- enrichr(target_genes, dbs)
		enriched <- do.call("rbind", enriched)
		enriched <- enriched[enriched$Adjusted.P.value < 0.5,]
		if(nrow(enriched) > 1){
				all_results <- rbind(all_results, 
								data.frame(Term=enriched$Term, 
										   Adjusted.P.value=enriched$Adjusted.P.value, 
										   Combined.Score=enriched$Combined.Score, 
										   Genes=enriched$Genes, 
										   driver=driver))
		}
	}
	return(all_results)
}


get_dotplot <- function(data){
		ggplot() + 
	geom_point(data, mapping= aes(y=Term, x=driver, size =Adjusted.P.value, fill=Combined.Score), color='black', shape=21) + 
	scale_fill_gradient2(low = '#5788c9', high = '#db3e25', mid ='white', midpoint = 0,)+
	scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
	scale_alpha(range = c(1, 0.1)) + theme_minimal() +   
	scale_size_continuous(name = 'Adjusted\n Pvalue', range = c(5, 1), breaks=c(0.01, 0.05, 0.1, 0.5), limits=c(0, 0.5)) +
	facet_wrap(~Comparison, nrow=1) +  
	guides(fill=guide_colorbar(title='Combined\nScore'), size=guide_legend(title='Adjusted\nP.value')) + 
	theme(
		axis.text.y=element_text(size=3),
		axis.title.y=element_blank(),
		axis.title.x=element_blank(),
		axis.text.x=element_text(size = 8, angle=90, vjust=1, hjust=1),
		legend.position='right',
		legend.key.size = unit(0.3, "cm"),
		legend.title=element_text(size=8), 
		legend.text=element_text(size=6))
}


df_weigths_ThreeWay <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_Weigths_ThreeWay.rds') # del5q Vs nondel Vs healthy
SimiC_weigths_PostNormal <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_Weigths_PostNormal.rds') # normal cells Pre vs Post
SimiC_weigths_PostDel3 <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SimiC_Weigths_PostDel3.rds') # MDS(x4) -> NR(pre) -> PR ->  -> NR(post)



setEnrichrSite("Enrichr")
dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021", 'KEGG_2021_Human', 'Reactome_2022')


ora_ThreeWay   <- get_enrichr_results_SimiC(df_weigths_ThreeWay, dbs)
ora_PostNormal <- get_enrichr_results_SimiC(SimiC_weigths_PostNormal, dbs)
ora_PostDel3   <- get_enrichr_results_SimiC(SimiC_weigths_PostDel3, dbs)

ora_ThreeWay$Comparison   <- 'del(5q)MDSVsnormalMDSVshealthy'
ora_PostNormal$Comparison <- 'HealthyVsnormalMDSVsPRVsCR'
ora_PostDel3$Comparison   <- 'del(5q)MDSVsdl(5q)NRVsdel(5q)PR'

all_ora <- rbind(ora_ThreeWay, ora_PostNormal, ora_PostDel3)

tmp <- split(all_ora , all_ora$Comparison)
WriteXLS::WriteXLS(tmp, ExcelFileName='/home/tereshkova/data/gserranos/MDS/Plots/PaperFigures/Supplementary_Table_Regulon_ORA_analysis.xlsx',  SheetNames= names(tmp), col.names = TRUE, row.names = FALSE, verbose = FALSE)

pdf('./Plots/PaperFigures/FigS4.pdf')

get_dotplot(rbind(ora_ThreeWay, ora_PostNormal, ora_PostDel3))
dev.off()





ora_ThreeWay_GO   <- get_enrichr_results_SimiC(df_weigths_ThreeWay, 'GO_Biological_Process_2021')
ora_PostNormal_GO <- get_enrichr_results_SimiC(SimiC_weigths_PostNormal, 'GO_Biological_Process_2021')
ora_PostDel3_GO   <- get_enrichr_results_SimiC(SimiC_weigths_PostDel3, 'GO_Biological_Process_2021')
all_ora_GO <- rbind(ora_ThreeWay, ora_PostNormal, ora_PostDel3, 'GO_Biological_Process_2021')



all_ora_only_GOs <- all_ora_GO[grepl('(GO:)', all_ora_GO$Term),]
all_ora_only_GOs$GOs <- stringr::str_extract(all_ora_only_GOs$Term, 'GO:[\\d]+')

mat = GO_similarity(all_ora_only_GOs$GOs, ont='BP')
df = simplifyGO(mat, plot = FALSE)

library(simplifyEnrichment)


GO:0034310 GO:0042574 GO:0016655 GO:0071398 GO:0050810
GO:0051055 GO:0033764 GO:0034035 GO:0050660 GO:0005929
GO:0050427 GO:0033619 GO:0004623 GO:0008374 GO:0052689
GO:0031253 GO:0046475 GO:0030263 GO:0051783 GO:0032395
GO:0051998 GO:0005021 GO:0071617 GO:0004862 GO:0008603
GO:0034236 GO:0042805 GO:0004435 GO:0042613 GO:0014731
GO:0002484 GO:0002486 GO:0046635 GO:0001916 GO:0060219
GO:0090207 GO:0030202 GO:1905605 GO:0032831 GO:0045345
GO:0007171 GO:0051133 GO:0045622 GO:0072311 GO:0061318
GO:0002692 GO:1902237 GO:2000698 GO:0001996 GO:0042270
GO:0060372 GO:1905603 GO:0014842 GO:0002664 GO:0002320
GO:0044091 GO:2000514 GO:2000480 GO:0043372 GO:0080154
GO:0072182 GO:0001912 GO:0043029 GO:2000344 GO:0090257
GO:0001754 GO:1904752 GO:0043116 GO:0043117 GO:0010676
GO:0060046 GO:2000479 GO:0016045 GO:0010907 GO:0010560
GO:0045953 GO:0062013 GO:0034199 GO:0010460 GO:0036149
GO:0071377 GO:0033762 GO:0043949 GO:2000352 GO:0045214
GO:0043552 GO:0003091 GO:0016747 GO:0030258 GO:0036018
GO:0071294 GO:0071280 GO:0071276 GO:0006882 GO:0046686
GO:0010043 GO:0055069 GO:0032415 GO:2000564 GO:1902099
GO:0001913 GO:2001185 GO:2001187 GO:0006825 GO:0035774
GO:0033617 GO:0032673 GO:0008535 GO:0061178 GO:0004691
GO:0004690 GO:0097413 GO:1903203 GO:0034380 GO:0007127
GO:0031089 GO:0016198 GO:0033007 GO:0033004 GO:0002887
GO:1905598 GO:1900272 GO:0072698 GO:0043301 GO:0035235
GO:1990535 GO:0016322 GO:0002281 GO:1900273 GO:0002695
GO:0042167 GO:0006787 GO:1905606 GO:0099174 GO:0042551
GO:0034142 GO:0007026 GO:0045821 GO:1900544 GO:1901880
GO:0045913 GO:0031114 GO:0010757 GO:0048711 GO:0090331
GO:0034111 GO:0048710 GO:0010544 GO:0006665 GO:0004181
GO:0070593 GO:0002003 GO:0002002 GO:0060048 GO:0050145
GO:0050699 GO:0009262 GO:0006564 GO:0006555 GO:0006563
GO:0009070 GO:1905268 GO:0031936 GO:0060969 GO:0031935
GO:0009165 GO:0009066 GO:1902895 GO:0004385 GO:0098837
GO:0005828 GO:0051648 GO:0098887 GO:0099639 GO:0098969
GO:0032402 GO:0051904 GO:0032401 GO:0032400 GO:1900271
GO:0046703 GO:0042583 GO:0002715 GO:0002727 GO:1901077
GO:0010459 GO:0002729 GO:0045899 GO:1903115 GO:0046135
GO:0002716 GO:0045822 GO:0009164 GO:0008253 GO:0008252
GO:0001405 GO:0072529 GO:0005796 GO:0051963 GO:0002040
GO:0050806 GO:0006024 GO:0009071 GO:0006544 GO:0032753
GO:0045954 GO:0002717 GO:0042269 GO:0010959 GO:1901385
GO:0010918 GO:0031639 GO:0045838 GO:0018193 GO:0018206
GO:0046394 GO:0031365 GO:0000338 GO:0010906 GO:0019318
GO:0004177 GO:0008238 GO:0008235 GO:0005887 GO:0098916
GO:0043922 GO:2000104 GO:0000077 GO:0016301 