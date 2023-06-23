#this script plots figure 5A: where x axis is whether an ORF shares a promoter with neighboring conserved gene 
# (either by sharing at least one transcript with neighboring cORF for up and down same, or by not having overlapping TSS for up opposite orfs)
# TIF seq data from Pelechano et al 2013.
#y axis is median expression level of de novo ORF across all samples where the ORF is detected
# graph is faceted by orientation 
# only considers de novo ORFs that are in one orientation

#tif_info data generated in get_tif_info.R
#expression_info generated in get_median_tpm.R
#OrfsInOneOrientation generated in get_orf_pairs_in_one_orientation.R
#orfsFarAway generated in get_far_away_orfs.R
#overlap_TSS_info generated in get_overlap_TSS_info.R

library(dplyr)
library(ggpubr)
library(magrittr)
library(forcats)
library(scales)
library(RMariaDB)
library(stringr)
library(EnvStats)
library(effsize)
orientation_levels <- c("down same", "up same", "up opposite")
significance_size <- 4 # default is 3.88 and does not scale with other size parameters, e.g., element_text(size=4) is different size
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
orfsFarAway<-dbGetQuery(conn,'SELECT * FROM april.OrfsFarAway_500bp')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
dbDisconnect(conn)

far_away_orfs<-filter(orf_pairs, primary_orf %in% orfsFarAway$transcript)
far_away_orfs<-left_join(far_away_orfs, expression_info, by=c('primary_orf'='transcript')) 
far_away_orfs<-filter(far_away_orfs, !is.na(tpm))
far_away_orfs<-far_away_orfs[,c('primary_orf','tpm')]
far_away_orfs<-unique(far_away_orfs)



rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

orf_pairs<-filter(orf_pairs, distance <=500 & orientation %in% c('downsame','upsame','upopposite') & primary_orf %in% rownames(rho) & neighbor_orf %in% rownames(rho))



overlap_TSS_info<-readRDS('TSS_overlap_info_prop.RDS')

orf_pairs<-left_join(orf_pairs,overlap_TSS_info)


tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0)
orf_pairs<-left_join(orf_pairs,tif_info)


orf_pairs$rho<-NA
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho[orf_pairs$primary_orf[i],orf_pairs$neighbor_orf[i]]
}


orf_pairs<-filter(orf_pairs, !is.na(rho))



#expression:
orf_pairs<-left_join(orf_pairs, expression_info, by=c('primary_orf'='transcript')) 

#path to save plots
file_path<-'./20230106_figures/'

plot_df<-orf_pairs %>%  filter(primary_orf %in% OrfsInOneOrientation$transcript & orientation %in% c('downsame', 'upsame', 'upopposite')) %>% 
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same',
                                orientation =='upopposite' ~'up opposite')) %>% 
  mutate(share_promoter = case_when((orientation =='down same' | orientation =='up same') & tif_intersect_counts >0 ~ 'yes',
                                    (orientation =='down same' | orientation=='up same' ) & tif_intersect_counts==0 ~'no',
                                    orientation == 'up opposite'  & num_time_share_promoter>=1 ~ 'yes',
                                    orientation == 'up opposite' & prop_time_share_promoter==0 & num_gene_tifs >0 ~ 'no')) %>%
  mutate(orientation=factor(orientation,levels=orientation_levels)) %>%
  filter(!is.na(share_promoter))     

p<-  plot_df %>% 
  ggboxplot(x="share_promoter",y="tpm", fill="orientation", xlab="share a promoter",ylab ="median expression (tpm)", 
            notch=TRUE,yscale = "log10", palette=c('#86B040','#BBB0D7', '#0B6C9A'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  stat_compare_means(aes(group = share_promoter),label = "p.signif", label.x = 1.45, label.y=2.15, size=significance_size)+
  geom_hline(yintercept=median(far_away_orfs$tpm),linetype='dotted')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")
pdf(sprintf('%sExpressionVssharePromoter_boxplot_1Tif_500.pdf',file_path),width=4.5,height=2.6)
ggpar(p)
dev.off()



wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='down same' &share_promoter=='no')$tpm)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='down same' & share_promoter=='no')$tpm)
#Cliffs delta: 0.75 large, p=1.06e-8

wilcox.test(filter(plot_df,orientation =='up same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up same' &share_promoter=='no')$tpm)
cliff.delta(filter(plot_df,orientation =='up same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up same' & share_promoter=='no')$tpm)
#Cliffs delta d=0.3799 medium p = 1.23e-7

wilcox.test(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up opposite' &share_promoter=='no')$tpm)
cliff.delta(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up opposite' & share_promoter=='no')$tpm)
#cliffs delta d=0.303 small p = 0.002982


#down same vs up same both share a promoter
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up same' &share_promoter=='yes')$tpm)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up same' &share_promoter=='yes')$tpm)
#Cliffs delta: 0.577 large, p < 2.2e-16

#down same vs up opposite both share a promoter
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up opposite' &share_promoter=='yes')$tpm)
cliff.delta(filter(plot_df,orientation =='up same' & share_promoter=='yes')$tpm, filter(plot_df,orientation =='up opposite' & share_promoter=='no')$tpm)
#Cliffs delta d=0.545 large p < 2.2e-16
