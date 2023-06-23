
#script to plot figure 4E: x axis is wether ORFs can piggyback or cannot piggyback 
#y axis median expression of denovo ORF across all samples 
#only considers ORFs that are in one orientation (ie are only within 500bp of one gene)
#list of ORFs that are in one orientation from script get_orf_pairs_in_one_orientation.R
#dashed line is median of independent orfs (ie orfs that are not within 500 bp of a conserved gene) and is calculated in script get_far_away_orfs.R
#median tpm for each ORF from script get_median_tpm.R

library(dplyr)
library(magrittr)
library(stringr)
library(ggpubr)
library(forcats)
library(scales)
library(RMariaDB)
library(EnvStats)
library(effsize)
#path to save plots
file_path<-'./20221110_rho/20230106_figures/'
orientation_levels <- c("independent", "down\nopposite", "antisense\noverlap", "up\nopposite", "up\nsame", "down\nsame")
significance_size <- 4 # default is 3.88 and does not scale with other size parameters, e.g., element_text(size=4) is different size
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


conn <- dbConnect(MariaDB(), 
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
median_tpm<-dbGetQuery(conn,'SELECT * FROM april.coexpression_median_tpm')
orfsInOneOrientation<-dbGetQuery(conn,'SELECT * FROM april.OrfsInOneOrientation_500bp')
orfsFarAway<-dbGetQuery(conn,'SELECT * FROM april.OrfsFarAway_500bp')
dbDisconnect(conn)

orfsFarAway<-orfsFarAway$transcript
orfsInOneOrientation<-orfsInOneOrientation$transcript

orf_pairs<-filter(orf_pairs, orientation !='sense-overlap')
orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down\nsame',
                                                   orientation=='upsame' ~'up\nsame',
                                                   orientation =='antisense-overlap' ~'antisense\noverlap',
                                                   orientation=='downopposite' ~'down\nopposite',
                                                   orientation=='upopposite' ~'up\nopposite'))


orf_pairs<-left_join(orf_pairs, median_tpm[,c('transcript','tpm')], by=c('primary_orf'='transcript'))

pairs_df<-filter(orf_pairs, distance <= 500 & primary_orf %in% orfsInOneOrientation)
far_away_orfs<-filter(orf_pairs, primary_orf %in% orfsFarAway)
far_away_orfs<-far_away_orfs[,c('primary_orf','tpm')]
far_away_orfs<-unique(far_away_orfs)
far_away_orfs<-filter(far_away_orfs, !is.na(tpm))
far_away_orfs$orientation<-'independent'

# pairs_df<-rbind(pairs_df_500[,c('primary_orf','tpm','orientation')], far_away_orfs)
pairs_df<-filter(pairs_df, !is.na(tpm))
# pairs_df$logtpm<-log(pairs_df$tpm+1)


#can piggyback vs cannot piggyback:
plot_df<-pairs_df %>% mutate(orientation_type=case_when(orientation %in% c('down\nsame','up\nsame','up\nopposite') ~ 'can\npiggyback',
                                                        orientation %in% c('down\nopposite','antisense\noverlap') ~ 'cannot\npiggyback'))
p<- plot_df%>% mutate(orientation_type=factor(orientation_type,levels=c('cannot\npiggyback','can\npiggyback'))) %>%
  ggboxplot(x = "orientation_type", y = "tpm", ylab = "median expression (TPM)", notch = TRUE, fill = "orientation_type",
            palette = c('#F6DFE0','#DBE8C6'),yscale = "log10") +
  geom_hline(yintercept=median(far_away_orfs$tpm),linetype='dotted')+
  stat_n_text(size = 2.75) +
  stat_compare_means(label = "p.signif", size = significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  # rotate_x_text(45)+
  rremove("xlab")+
  rremove("legend")

pdf(sprintf('%sExpressionVsOrientation_type.pdf',file_path),width=2.25,height=2.5)
ggpar(p)
dev.off()

wilcox.test(filter(plot_df,orientation_type =='can\npiggyback')$tpm, filter(plot_df,orientation_type =='cannot\npiggyback')$tpm)
cliff.delta(filter(plot_df,orientation_type =='can\npiggyback')$tpm, filter(plot_df,orientation_type =='cannot\npiggyback')$tpm)
#Cliffs delta:0.395 medium , p<2.2e-16



