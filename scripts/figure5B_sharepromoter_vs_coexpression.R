#this script plots figure 5B: where x axis is whether an ORF shares at promoter with neighboring cORF
#y axis is coexpression of de novo ORF with neighboring conserved gene 
#graph is faceted by orientation ie either for de novo ORFs that are upstream of a conserved gene or downstream of a conserved gene
# TIF seq data from Pelechano et al 2013.
#tif_info data generated in get_tif_info.R
#overlap_TSS_info generated in get_overlap_tss_info.R

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



#path to save plots
file_path<-'./20230106_figures/'

background<-readRDS('rhoForDCPairsNotOnSameChr.RDS')
colnames(background)[3]<-'rho'
background$orientation<-'background'



plot_df<-orf_pairs %>%  filter(orientation %in% c('downsame', 'upsame', 'upopposite')) %>% 
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same',
                                orientation =='upopposite' ~'up opposite')) %>% 
  mutate(share_promoter = case_when((orientation =='down same' | orientation =='up same') & tif_intersect_counts >0 ~ 'yes',
                                    (orientation =='down same' | orientation=='up same' ) & tif_intersect_counts==0 ~'no',
                                    orientation == 'up opposite'  & num_time_share_promoter>=1 ~ 'yes',
                                    orientation == 'up opposite' & prop_time_share_promoter==0 & num_gene_tifs >0 ~ 'no')) %>%
  mutate(orientation=factor(orientation,levels=orientation_levels)) %>%
  filter(!is.na(share_promoter))     

p<-plot_df %>% 
  ggboxplot(x="share_promoter",y="rho", fill="orientation", xlab="share a promoter",ylab ="coexpression with neighbor", 
            notch=TRUE, palette=c('#86B040','#BBB0D7', '#0B6C9A'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  stat_compare_means(aes(group = share_promoter),label = "p.signif", label.x = 1.45, label.y=.9, size=significance_size)+
  geom_hline(yintercept=median(background$rho),linetype='dotted')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")
pdf(sprintf('%sCoexpressionVssharePromoter_boxplot_500.pdf',file_path),width=4.5,height=2.6)
ggpar(p)
dev.off()


wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='down same' & share_promoter=='no')$rho)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='down same' & share_promoter=='no')$rho)
#Cliffs delta: 0.283 small, p=2.994e-09

wilcox.test(filter(plot_df,orientation =='up same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up same' &share_promoter=='no')$rho)
cliff.delta(filter(plot_df,orientation =='up same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up same' & share_promoter=='no')$rho)
#Cliffs delta d=0.309 small p < 2.2e-16

wilcox.test(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up opposite' &share_promoter=='no')$rho)
cliff.delta(filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up opposite' & share_promoter=='no')$rho)
#cliffs delta d=0.27 small p = 2.105e-07

#down same vs up same both share promoter
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up same' & share_promoter=='yes')$rho)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up same' & share_promoter=='yes')$rho)
#Cliffs delta: 0.294 small, p <2.2e-16

#down same vs up oppsote both share a promoter:
wilcox.test(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rho)
cliff.delta(filter(plot_df,orientation =='down same' & share_promoter=='yes')$rho, filter(plot_df,orientation =='up opposite' & share_promoter=='yes')$rho)
#Cliffs delta: 0.377 medium, p <2.2e-16
