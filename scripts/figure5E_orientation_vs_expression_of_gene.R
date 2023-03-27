#this script plots figure 5E: where x axis is category of gene (either a gene that has a downstream same de novo ORF or a gene that has a upstream same de novo ORF )
#y axis is median expression level of the conserved gene across all rna seq samples where the gene is detected
# only considers genes that are neighboring (within 500bp) de novo ORFs that are in one orientation

#tif_info data generated in get_tif_info.R
#expression_info generated in get_median_tpm.R
#OrfsInOneOrientation generated in get_orf_pairs_in_one_orientation.R
#orfsFarAway generated in get_far_away_orfs.R


library(dplyr)
library(ggpubr)
library(magrittr)
library(forcats)
library(scales)
library(RMariaDB)
library(stringr)
library(EnvStats)
library(effsize)

orientation_levels <- c("far\naway", "down\nopposite", "down\nsame", "overlap", "up\nopposite", "up\nsame")
significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

#get ORF data
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
expression_info<-dbGetQuery(conn, "SELECT * FROM april.coexpression_median_tpm")
tif_info<-dbGetQuery(conn, "SELECT * FROM april.TIF_info")
OrfsInOneOrientation<-dbGetQuery(conn, "SELECT * FROM april.OrfsInOneOrientation_500bp")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0)
tif_info<-filter(tif_info,primary_orf %in% rownames(rho))
tif_info<-filter(tif_info,neighbor_orf %in% rownames(rho))
tif_info$rho<-NA
for(i in 1:nrow(tif_info)){
  tif_info$rho[i]<-rho[tif_info$primary_orf[i],tif_info$neighbor_orf[i]]
}

tif_info<-filter(tif_info, !is.na(rho))



x<-table(tif_info$neighbor_orf) #236 genes which have more than one de novo ORF within 500 bp 

genesWithdORFs<-filter(tif_info, orientation=='downsame')$neighbor_orf
genesWithuORFs<-filter(tif_info, orientation=='upsame')$neighbor_orf
geneswithBoth<-intersect(genesWithuORFs, genesWithdORFs)
tif_info<-filter(tif_info, !(neighbor_orf %in% geneswithBoth))

tif_info<-left_join(tif_info, expression_info,by=c('neighbor_orf'='transcript'))

genes<-tif_info[,c('neighbor_orf','orientation','tpm')]
genes<-unique(genes) #1508


file_path<-'./20230106_figures/'

pdf(sprintf('%sdORFvsuORF_gene_expression.pdf',file_path),width=2.2,height=3)
genes %>% filter(neighbor_orf %in% filter(tif_info, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf) %>%
  mutate(orientation=case_when(orientation=='downsame' ~'genes with\ndownsame\nORFs', 
                               orientation =='upsame' ~ 'genes with\nupsame\nORFs')) %>%
  ggboxplot(x="orientation",y="tpm",notch=TRUE, "ylab"="median expression",yscale = "log10")+
  stat_n_text(size = 2.75)+
  stat_compare_means(label = "p.signif", label.y = 4.25, label.x = 1.45, size=significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("legend")
dev.off()



# 
wilcox.test(filter(genes,orientation=='downsame' & neighbor_orf %in% filter(tif_info, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='upsame' & neighbor_orf %in% filter(tif_info, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
cliff.delta(filter(genes,orientation=='downsame' & neighbor_orf %in% filter(tif_info, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='upsame' & neighbor_orf %in% filter(tif_info, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
# #cliffs delta d=0.1877 (small) p= 0.00428
