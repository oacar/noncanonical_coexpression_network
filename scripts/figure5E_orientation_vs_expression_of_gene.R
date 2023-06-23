#this script plots figure 5E: where x axis is category of gene (either a gene that has a downstream same de novo ORF or a gene that has a upstream same or opposite de novo ORF )
#y axis is median expression level of the conserved gene across all rna seq samples where the gene is detected
# only considers genes that are neighboring (within 500bp) de novo ORFs that are in one orientation

#tif_info data generated in get_tif_info.R
#overlap_TSS_info generated in get_overlap_TSS_info.R
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


significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
terminator_info<-dbGetQuery(conn, "SELECT * FROM april.Terminator_info")
expression_info<-dbGetQuery(conn, "SELECT * FROM april.coexpression_median_tpm")
tif_info<-dbGetQuery(conn, "SELECT * FROM april.TIF_info")
OrfsInOneOrientation<-dbGetQuery(conn, "SELECT * FROM april.OrfsInOneOrientation_500bp")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
gsea<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_similarity')
background_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_bp_noexponent_similarity_random')
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
dbDisconnect(conn)

rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

orf_pairs<-filter(orf_pairs, distance <=500 & orientation %in% c('downsame','upsame','upopposite') & primary_orf %in% rownames(rho) & neighbor_orf %in% rownames(rho))


overlap_TSS_info<-readRDS('TSS_overlap_info_prop.RDS')
overlap_TSS_info<-filter(overlap_TSS_info,distance <=500 & !is.na(num_gene_tifs))

orf_pairs<-left_join(orf_pairs, insulator_info)
orf_pairs<-left_join(orf_pairs,overlap_TSS_info)


tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0)
orf_pairs<-left_join(orf_pairs,tif_info)


orf_pairs$rho<-NA
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho[orf_pairs$primary_orf[i],orf_pairs$neighbor_orf[i]]
}


orf_pairs<-filter(orf_pairs, !is.na(rho))


orf_pairs<-orf_pairs %>% mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                                       orientation =='upsame' ~'up same',
                                                       orientation =='upopposite' ~'up opposite')) %>% 
  mutate(share_promoter = case_when((orientation =='down same' | orientation =='up same') & tif_intersect_counts >0 ~ 'yes',
                                    (orientation =='down same' | orientation=='up same' ) & tif_intersect_counts==0 ~'no',
                                    orientation == 'up opposite'  & num_time_share_promoter>=1 ~ 'yes',
                                    orientation == 'up opposite' & prop_time_share_promoter==0 ~ 'no'))

orf_pairs<-filter(orf_pairs, share_promoter=='yes')

x<-table(orf_pairs$neighbor_orf) 
genesInMoreThanOneOrientation<-names(which(x>1))

genes2remove<-c()
for(i in 1:length(genesInMoreThanOneOrientation)){
  xx<-filter(orf_pairs, neighbor_orf ==genesInMoreThanOneOrientation[i])
  if(length(unique(xx$orientation))>1){genes2remove<-c(genes2remove,genesInMoreThanOneOrientation[i])}
}

orf_pairs<-filter(orf_pairs, !(neighbor_orf %in% genes2remove))


orf_pairs<-left_join(orf_pairs, expression_info,by=c('neighbor_orf'='transcript'))


genes<-orf_pairs[,c('neighbor_orf','orientation','tpm')]
genes<-unique(genes) 

orientation_levels <- c("down same", "up same", "up opposite")
file_path<-'./20230106_figures/'

pdf(sprintf('%sgene_expression.pdf',file_path),width=2.8,height=3)
genes %>% filter(neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf) %>%
  mutate(orientation = factor(orientation, levels = orientation_levels)) %>%
  # mutate(orientation=case_when(orientation=='down same' ~'genes with\ndown same\nORFs', 
  #                              orientation =='up same' ~ 'genes with\nup same\nORFs',
  #                              orientation =='up opposite' ~ 'genes witth\nup opposite\nORFs')) %>%
  ggboxplot(x="orientation",y="tpm",notch=TRUE, ylab="median expression",yscale = "log10", fill= "orientation",
            palette = c('#86B040','#BBB0D7', '#0B6C9A'))+
  stat_n_text(size = 2.75)+
  stat_compare_means(label = "p.signif", ref.group = "down same", label.y = 4, size = significance_size) +
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("legend")
dev.off()


wilcox.test(filter(genes,orientation=='down same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='up same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
cliff.delta(filter(genes,orientation=='down same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='up same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
# cliffs delta d=0.2 (small) p= 5.4e-3

wilcox.test(filter(genes,orientation=='down same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='up opposite' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
cliff.delta(filter(genes,orientation=='down same' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm,
            filter(genes, orientation =='up opposite' & neighbor_orf %in% filter(orf_pairs, primary_orf %in% OrfsInOneOrientation$transcript)$neighbor_orf)$tpm)
# cliffs delta d=0.342 (medium) p= 6.5e-4