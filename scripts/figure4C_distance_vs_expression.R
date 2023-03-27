#script to figure 4C: x axis distance, y axis median expression across all samples, faceted by orientation
#only considers ORFs that are in one orientation (ie are only within 500bp of one gene)
#list of ORFs that are in one orientation from script get_orf_pairs_in_one_orientation.R
#list of independent orfs (ie orfs that are not within 500 bp of a conserved gene) from script get_far_away_orfs.R
#median tpm for each ORF from script get_median_tpm.R
library(dplyr)
library(ggpubr)
library(magrittr)
library(forcats)
library(scales)
library(RMariaDB)
library(stringr)

orientation_levels <- c("far\naway", "down\nopposite", "down\nsame", "overlap", "up\nopposite", "up\nsame")
significance_size <- 4
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

conn <- dbConnect(MariaDB(), 
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
orfsInOneOrientation_500<-dbGetQuery(conn,'SELECT * FROM april.OrfsInOneOrientation_500bp')
median_tpm<-dbGetQuery(conn,'SELECT * FROM april.coexpression_median_tpm')
dbDisconnect(conn)

orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down same',
                                                   orientation=='upsame' ~'up same',
                                                   orientation =='antisense-overlap' ~'overlap',
                                                   orientation=='downopposite' ~'down opposite',
                                                   orientation=='upopposite' ~'up opposite'))
#path to save plots
file_path<-'./20230106_figures/'

orf_pairs<-left_join(orf_pairs, median_tpm[,c('transcript','tpm')], by=c('primary_orf'='transcript'))

#500 bp away
pairs_df_500<-filter(orf_pairs, distance <= 500 & primary_orf %in% orfsInOneOrientation_500$transcript)

pdf(sprintf("%sexpression_distance_lineplot_500bp.pdf", file_path), width=8,height=2.5)
pairs_df_500 %>% filter(orientation %in% c('down same','up same','up opposite','down opposite') & !is.na(tpm)) %>%
  group_by(orientation) %>%
  mutate(orientation = paste0( orientation, "\nn=", n())) %>%
ggplot(aes(x=distance, y=tpm))+geom_point(alpha=0.8,color = "black", shape = 21, size = 1) + facet_wrap(~orientation,ncol=4)+ theme_minimal()+
  labs(y='median expression (tpm)', x='distance between de novo ORF and closest conserved ORF')+
  theme(axis.text=element_text(size=axis_text_size), axis.title=element_text(size=axis_title_size))+
  geom_smooth(method='lm',color="#76ABDD")+ #adds line with error region
  ggpubr::stat_cor(size=2.75,label.x = 150, method="spearman")+
  scale_y_continuous(trans = pseudo_log_trans(base= 10),breaks=c(1,10,100))+
  scale_x_continuous(breaks=c(0,250,500))
dev.off()
