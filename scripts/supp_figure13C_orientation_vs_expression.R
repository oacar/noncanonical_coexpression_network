#script to plot supp figure 13C: x axis orientation of de novo ORF relative to conserved gene,
#y axis median expression of denovo ORF across all samples 
#only considers ORFs that are in one orientation (ie are only within 500bp of one gene)
#list of ORFs that are in one orientation from script get_orf_pairs_in_one_orientation.R
#list of independent orfs (ie orfs that are not within 500 bp of a conserved gene) from script get_far_away_orfs.R
#median tpm for each ORF from script get_median_tpm.R

library(dplyr)
library(magrittr)
library(stringr)
library(ggpubr)
library(forcats)
library(scales)
library(RMariaDB)
library(EnvStats)
#path to save plots
file_path<-'./20221110_rho/20230106_figures/'
orientation_levels <- c("independent", "down\nopposite", "down\nsame", "antisense\noverlap", "up\nopposite", "up\nsame")
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
median_tpm<-dbGetQuery(conn,'SELECT * FROM april.coexpression_median_tpm')
orfsInOneOrientation<-dbGetQuery(conn,'SELECT * FROM april.OrfsInOneOrientation_500bp')
orfsFarAway<-dbGetQuery(conn,'SELECT * FROM april.OrfsFarAway_500bp')
dbDisconnect(conn)

orfsFarAway<-orfsFarAway$transcript
orfsInOneOrientation<-orfsInOneOrientation$transcript

orf_pairs<-filter(orf_pairs, orientation !='sense-overlap')
orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down same',
                                                   orientation=='upsame' ~'up same',
                                                   orientation =='antisense-overlap' ~'antisense overlap',
                                                   orientation=='downopposite' ~'down opposite',
                                                   orientation=='upopposite' ~'up opposite'))


orf_pairs<-left_join(orf_pairs, median_tpm[,c('transcript','tpm')], by=c('primary_orf'='transcript'))

pairs_df_500<-filter(orf_pairs, distance <= 500 & primary_orf %in% orfsInOneOrientation)
far_away_orfs<-filter(orf_pairs, primary_orf %in% orfsFarAway)
far_away_orfs<-far_away_orfs[,c('primary_orf','tpm')]
far_away_orfs<-unique(far_away_orfs)
far_away_orfs$orientation<-'independent'

pairs_df<-rbind(pairs_df_500[,c('primary_orf','tpm','orientation')], far_away_orfs)
pairs_df<-filter(pairs_df, !is.na(tpm))



pdf(sprintf('%sorientationVsExpression.pdf',file_path),width=3.6, height = 3)
pairs_df %>% 
  mutate(orientation = factor(str_replace(orientation, "\ ", "\n"), levels = orientation_levels)) %>%
  ggboxplot(x="orientation",y="tpm",ylab ="median expression (tpm)", notch=TRUE,fill='#76ABDD',yscale = "log10")+
  geom_hline(yintercept=median(filter(pairs_df,orientation=='independent')$tpm),linetype='dotted')+
  stat_n_text(size = 2.75) + 
  stat_compare_means(label="p.signif", ref.group="independent",label.y=2.5,size=significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("xlab")
dev.off()


#down same
o<-'down same'
wilcox.test(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
cliff.delta(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
#down same: large (Cliff’s Delta d=0.58, Mann-Whitney U-test, p< 2.2e-16)

o<-'up same'
wilcox.test(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
cliff.delta(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
#up same: negligible (Cliff’s Delta d=0.077, Mann-Whitney U-test, p=0.109)

o<-'antisense overlap'
wilcox.test(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
cliff.delta(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
#antisense overlap: large (Cliff’s Delta d=-0.52, Mann-Whitney U-test, p < 2.2e-16)

o<-'up opposite'
wilcox.test(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
cliff.delta(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
#up opposite:small (Cliff’s Delta d=-0.22, Mann-Whitney U-test, p=2.616e-4)

o<-'down opposite'
wilcox.test(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
cliff.delta(filter(pairs_df,orientation==o)$tpm, filter(pairs_df,orientation=='independent')$tpm)
#down opposite: negligible (Cliff’s Delta d=0.027, Mann-Whitney U-test, p=0.673)




