#this script plots figure 5A: where x axis is whether an ORF shares at least one transcript with neighboring
# conserved gene as determined by TIF seq data from Pelechano et al 2013.
#y axis is median expression level of de novo ORF across all samples where the ORF is detected
#graph is faceted by orientation ie either for de novo ORFs that are upstream of a conserved gene or downstream of a conserved gene
# only considers de novo ORFs that are in one orientation

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
orfsFarAway<-dbGetQuery(conn,'SELECT * FROM april.OrfsFarAway_500bp')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
dbDisconnect(conn)

rho<-readRDS('/home/aar75/coexpression/20221110_rho/spqn_raw5_sample400.RDS')
num_obs<-readRDS('/home/aar75/coexpression/20221110_rho/numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0)
tif_info<-filter(tif_info,primary_orf %in% rownames(rho))
tif_info<-filter(tif_info,neighbor_orf %in% rownames(rho))


tif_info$rho<-NA
for(i in 1:nrow(tif_info)){
  tif_info$rho[i]<-rho[tif_info$primary_orf[i],tif_info$neighbor_orf[i]]
}


tif_info<-filter(tif_info, !is.na(rho))

far_away_orfs<-filter(orf_pairs, primary_orf %in% orfsFarAway$transcript)


far_away_orfs<-left_join(far_away_orfs, expression_info, by=c('primary_orf'='transcript')) 

far_away_orfs<-filter(far_away_orfs, !is.na(tpm))
far_away_orfs<-far_away_orfs[,c('primary_orf','tpm')]
far_away_orfs<-unique(far_away_orfs)


#expression:
tif_info<-left_join(tif_info, expression_info, by=c('primary_orf'='transcript')) 

#path to save plots
file_path<-'./20230106_figures/'


p<-tif_info %>%  filter(primary_orf %in% OrfsInOneOrientation$transcript & orientation %in% c('downsame', 'upsame')) %>% 
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same')) %>% 
  mutate(share_tif= case_when( tif_intersect_counts>0 ~'yes', tif_intersect_counts==0 ~'no')) %>% 
  ggboxplot(x="share_tif",y="tpm", fill="orientation", xlab="share > 1 transcript",ylab ="median expression (tpm)", 
            notch=TRUE,yscale = "log10", palette=c('#86B040','#BBB0D7'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  stat_compare_means(aes(group = share_tif),label = "p.signif", label.x = 1.45, label.y=2.3, size=significance_size)+
  geom_hline(yintercept=median(far_away_orfs$tpm),linetype='dotted')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")
pdf(sprintf('%sExpressionVsshareTIF_boxplot_1Tif_500.pdf',file_path),width=2.75,height=3.2)
ggpar(p)
dev.off()



x<-tif_info %>%  filter(primary_orf %in% OrfsInOneOrientation$transcript & orientation =='downsame') %>% 
  mutate(share_tif= case_when( tif_intersect_counts>0 ~'yes', tif_intersect_counts==0 ~'no'))
wilcox.test(filter(x,share_tif=='yes')$tpm, filter(x,share_tif=='no')$tpm)
cliff.delta(filter(x,share_tif=='yes')$tpm, filter(x,share_tif=='no')$tpm)
#Cliffs delta: 0.75 large, p=1.06e-8

x<-tif_info %>%  filter(primary_orf %in% OrfsInOneOrientation$transcript & orientation =='upsame') %>% 
  mutate(share_tif= case_when( tif_intersect_counts>0 ~'yes', tif_intersect_counts==0 ~'no'))
wilcox.test(filter(x,share_tif=='yes')$tpm, filter(x,share_tif=='no')$tpm)
cliff.delta(filter(x,share_tif=='yes')$tpm, filter(x,share_tif=='no')$tpm)
