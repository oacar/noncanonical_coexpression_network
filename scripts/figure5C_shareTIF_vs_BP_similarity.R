#this script plots figure 5C: where x axis is whether an ORF shares at least one transcript with neighboring
# conserved gene as determined by TIF seq data from Pelechano et al 2013.
#y axis is BP similarity of de novo ORF with neighboring conserved gene 
#graph is faceted by orientation ie either for de novo ORFs that are upstream of a conserved gene or downstream of a conserved gene

#tif_info data generated in get_tif_info.R
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
tif_info<-dbGetQuery(conn, "SELECT * FROM april.TIF_info")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
gsea<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_similarity')
background_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_bp_noexponent_similarity_random')
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

file_path<-'./20230106_figures/'

#BP similarity

tif_info<-left_join(tif_info,gsea[,c('primary_orf','neighbor_orf','rel')])
idx<- which(is.na(tif_info$rho) & !is.na(tif_info$rel))
tif_info$rel[idx]<-NA
idx<-which(is.na(tif_info$rel) & !is.na(tif_info$rho))
tif_info$rel[idx]<-0

wilcox.test(filter(tif_info,orientation=='downsame' & tif_intersect_counts>0)$rel, filter(tif_info, orientation =='downsame' & tif_intersect_counts==0)$rel)
cliff.delta(filter(tif_info,orientation=='downsame' & tif_intersect_counts>0)$rel, filter(tif_info, orientation =='downsame' & tif_intersect_counts==0)$rel)
#down same: cliffs delta d=0.239 (small) p=5.468e-07

wilcox.test(filter(tif_info,orientation=='upsame' & tif_intersect_counts>0)$rel, filter(tif_info, orientation =='upsame' & tif_intersect_counts==0)$rel)
cliff.delta(filter(tif_info,orientation=='upsame' & tif_intersect_counts>0)$rel, filter(tif_info, orientation =='upsame' & tif_intersect_counts==0)$rel)
#upsame: cliffs delta d=0.108 (negligible) p=0.00378


p<-tif_info %>% filter(orientation %in% c('downsame','upsame')) %>%
  mutate(share_tif= case_when( tif_intersect_counts>0 ~'yes', tif_intersect_counts==0 ~'no')) %>%
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same')) %>% 
  ggboxplot(x="share_tif",y="rel",xlab="share > 1 transcript",ylab ="BP similarity",fill="orientation",
            notch=TRUE, palette=c('#86B040','#BBB0D7'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  stat_compare_means(aes(group = share_tif),label = "p.signif", label.x = 1.45, size=significance_size)+
  geom_hline(yintercept=median(background_pairs$rel),linetype='dotted')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")
pdf(sprintf('%sBPSimilarityVsshareTIF_boxplot_1TIF_500.pdf',file_path),width=2.5,height=3.2)
ggpar(p)
dev.off()


