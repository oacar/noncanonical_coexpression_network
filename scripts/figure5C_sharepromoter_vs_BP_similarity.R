#this script plots figure 5C: where x axis is whether an ORF shares at promoter with neighboring cORF
#graph is faceted by orientation ie either for de novo ORFs that are upstream of a conserved gene or downstream of a conserved gene
#y axis is BP similarity of de novo ORF with neighboring conserved gene 
# TIF seq data from Pelechano et al 2013.


#tif_info data generated in get_tif_info.R
#overlap_TSS_info generated in get_overlap_TSS_info.R
#BP enrichment similarity (omer.fgsea_slim_similarity) generated in gsea_relevance_similarity.R
#background BP enrichment similarity (omer.fgsea_slim_bp_noexponent_similarity_random) generated in gsea_jaccard_random.r

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
overlap_TSS_info<-filter(overlap_TSS_info,distance<=500 & !is.na(num_gene_tifs))

orf_pairs<-left_join(orf_pairs, insulator_info)
orf_pairs<-left_join(orf_pairs,overlap_TSS_info)


tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0)
orf_pairs<-left_join(orf_pairs,tif_info)


orf_pairs$rho<-NA
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho[orf_pairs$primary_orf[i],orf_pairs$neighbor_orf[i]]
}


orf_pairs<-filter(orf_pairs, !is.na(rho))


file_path<-'./20230106_figures/'

#BP similarity

orf_pairs<-left_join(orf_pairs,gsea[,c('primary_orf','neighbor_orf','rel')])
idx<- which(is.na(orf_pairs$rho) & !is.na(orf_pairs$rel))
orf_pairs$rel[idx]<-NA
idx<-which(is.na(orf_pairs$rel) & !is.na(orf_pairs$rho))
orf_pairs$rel[idx]<-0


plot_df<-orf_pairs %>%  filter(orientation %in% c('downsame', 'upsame', 'upopposite')) %>% 
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same',
                                orientation =='upopposite' ~'up opposite')) %>% 
  mutate(share_promoter = case_when((orientation =='down same' | orientation =='up same') & tif_intersect_counts >0 ~ 'yes',
                                    (orientation =='down same' | orientation=='up same' ) & tif_intersect_counts==0 ~'no',
                                    orientation == 'up opposite'  & num_time_share_promoter>=1 ~ 'yes',
                                    orientation == 'up opposite' & prop_time_share_promoter==0 ~ 'no')) %>%
  mutate(orientation=factor(orientation,levels=orientation_levels)) %>%
  filter(!is.na(share_promoter))     

p<-plot_df %>% 
  ggboxplot(x="share_promoter",y="rel",xlab="share promoter",ylab ="BP similarity",fill="orientation",
            notch=TRUE, palette=c('#86B040','#BBB0D7', '#0B6C9A'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  stat_compare_means(aes(group = share_promoter),label = "p.signif", label.x = 1.45, label.y=.9, size=significance_size)+
  geom_hline(yintercept=median(background_pairs$rel),linetype='dotted')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")
pdf(sprintf('%sBPSimilarityVssharePromoter_boxplot_500.pdf',file_path),width=4.5,height=2.6)
ggpar(p)
dev.off()

wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='down same' & share_promoter=='no')$rel)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='down same' & share_promoter=='no')$rel)
#Cliffs delta: 0.2392547  small, p= 5.468e-07

wilcox.test(filter(plot_df,orientation =='up same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up same' &share_promoter=='no')$rel)
cliff.delta(filter(plot_df,orientation =='up same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up same' & share_promoter=='no')$rel)
#Cliffs delta d=0.1087984  negligible p =0.003781

wilcox.test(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up opposite' &share_promoter=='no')$rel)
cliff.delta(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up opposite' & share_promoter=='no')$rel)
#cliffs delta d= 0.235 small p =  6.056e-06


#down same vs up same both share a promoter
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up same' & share_promoter=='yes')$rel)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up same' & share_promoter=='yes')$rel)
#Cliffs delta: 0.369 medium, p= <2.2e-16

#down same vs up opposite both share a promoter
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rel)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rel, filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rel)
#Cliffs delta: 0.445 medium, p= <2.2e-16






