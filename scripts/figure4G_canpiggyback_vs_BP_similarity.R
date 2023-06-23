#script to plot figure 4G: boxplots of orientation type (ie either can or cannot piggyback) vs biological process similarity
#omer.gsea_similarity_random data is generated in gsea_jaccard_random.r
# omer.fgseas_slim_similarity data is generated in gsea_relevance_similarity.R

library(dplyr)
library(magrittr)
library(stringr)
library(ggpubr)
library(forcats)
library(scales)
library(EnvStats)
library(RMariaDB)
library(effsize)
#path to save plots
file_path<-'./20230106_figures/'


orientation_levels <- c("background", "down\nopposite", "antisense\noverlap", "up\nopposite", "up\nsame", "down\nsame")
significance_size <- 4 # default is 3.88 and does not scale with other size parameters, e.g., element_text(size=4) is different size
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

conn <- dbConnect(MariaDB(), 
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
background_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_bp_noexponent_similarity_random')
orf_pairs<-dbGetQuery(conn,'SELECT * FROM omer.fgsea_slim_similarity')
dbDisconnect(conn)
background_pairs<-data.frame('orientation'='background', 'rel'=background_pairs)

orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down same',
                                                   orientation=='upsame' ~'up same',
                                                   orientation =='antisense-overlap' ~'antisense overlap',
                                                   orientation=='downopposite' ~'down opposite',
                                                   orientation=='upopposite' ~'up opposite'))
rho_matrix<-readRDS('spqn_raw5_sample400.RDS')
numObs<-readRDS('numobs_raw5_sample400.RDS')
rho_matrix[numObs<400]<-NA
orf_pairs$rho<-NA
orf_pairs<-filter(orf_pairs, distance <=500)
orf_pairs<-filter(orf_pairs, primary_orf %in% rownames(rho_matrix))
orf_pairs<-filter(orf_pairs, neighbor_orf %in% rownames(rho_matrix))
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho_matrix[orf_pairs$neighbor_orf[i],orf_pairs$primary_orf[i]]
}
idx<- which(is.na(orf_pairs$rho) & !is.na(orf_pairs$rel))
orf_pairs$rel[idx]<-NA
idx<-which(is.na(orf_pairs$rel) & !is.na(orf_pairs$rho))
orf_pairs$rel[idx]<-0
orf_pairs<-orf_pairs[,c('orientation','rel')]

orf_pairs<-rbind(orf_pairs,background_pairs)
orf_pairs<-filter(orf_pairs, !is.na(rel))


#can piggyback vs cannot piggyback:
plot_df<-orf_pairs %>% filter(orientation !='background') %>% 
  mutate(orientation_type=case_when(orientation %in% c('down same','up same','up opposite') ~ 'can\npiggyback',
                                    orientation %in% c('down opposite','antisense overlap') ~ 'cannot\npiggyback'))
p<- plot_df%>%
  ggboxplot(x = "orientation_type", y = "rel", ylab = "biological process similarity with neighbor", notch = TRUE, fill = "orientation_type",
            palette = c('#F6DFE0','#DBE8C6')) +
  geom_hline(yintercept=median(filter(orf_pairs,orientation=='background')$rel),linetype='dotted')+
  stat_n_text(size = 2.75) +
  stat_compare_means(label = "p.signif", label.y = 1, size = significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  #rotate_x_text(45)+
  rremove("xlab")+
  rremove("legend")

pdf(sprintf('%sBPsimilarityVsOrientation_type.pdf',file_path),width=2.25,height=2.5)
ggpar(p)
dev.off()

wilcox.test(filter(plot_df,orientation_type =='can\npiggyback')$rel, filter(plot_df,orientation_type =='cannot\npiggyback')$rel)
cliff.delta(filter(plot_df,orientation_type =='can\npiggyback')$rel, filter(plot_df,orientation_type =='cannot\npiggyback')$rel)
#Cliffs delta: , p<

