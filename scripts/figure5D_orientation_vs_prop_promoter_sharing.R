#this script plots figure 5D: where x axis orientation of de novo ORF relative to conserved gene 
#y axis is proportion of time a de novo ORF shares a promoter with neighboring conserved gene as determined by TIF seq data from Pelechano et al 2013.
#ie y axis is number of transcripts containing de novo ORF that shares promoter with neighboring gene / number of tifs containing de novo ORF 
#tif_info data generated in get_tif_info.R
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
significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

#get ORF data
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

orf_pairs<-left_join(orf_pairs,overlap_TSS_info)

tif_info<-tif_info %>% filter(distance <= 500 & gene_tif_total_counts>0 )
orf_pairs<-left_join(orf_pairs,tif_info)


orf_pairs$rho<-NA
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho[orf_pairs$primary_orf[i],orf_pairs$neighbor_orf[i]]
}


orf_pairs<-filter(orf_pairs, !is.na(rho))


plot_df<-orf_pairs %>%  filter(orientation %in% c('downsame', 'upsame', 'upopposite')) %>% 
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same',
                                orientation =='upopposite' ~'up opposite')) %>% 
  mutate(prop_promoter_sharing = case_when((orientation =='down same' | orientation =='up same') & tif_intersect_counts ==0  ~ 0,
                                           (orientation =='down same' | orientation=='up same' ) & tif_intersect_counts>0 ~ denovo_tif_ratio_counts,
                                           orientation == 'up opposite'  & num_time_share_promoter>=1 ~ prop_time_share_promoter,
                                           orientation == 'up opposite' & prop_time_share_promoter==0 & num_gene_tifs >0 ~ 0)) %>%
  mutate(orientation=factor(orientation,levels=orientation_levels)) %>%
  filter(!is.na(prop_promoter_sharing))     


#path to save plots
file_path<-'./20230106_figures/'

pdf(sprintf('%sshareTIFProp_Orientation.pdf',file_path),width=2.8,height=3)
plot_df %>% mutate(orientation = factor(str_replace(orientation, " ", "\n"), levels = c('down\nsame','up\nsame','up\nopposite'))) %>%
  ggboxplot(x="orientation",y="prop_promoter_sharing", notch=TRUE,ylab="proportion of de novo ORF\ntranscripts that share promoter",
            fill="orientation",palette=c('#86B040','#BBB0D7', '#0B6C9A')) +
  stat_n_text(size = 2.75)+
  stat_compare_means(label = "p.signif", ref.group = "down\nsame", size = significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("legend")
dev.off()


wilcox.test(filter(plot_df,orientation=='down same')$prop_promoter_sharing, filter(plot_df,orientation=='up same')$prop_promoter_sharing)
cliff.delta(filter(plot_df,orientation=='down same')$prop_promoter_sharing, filter(plot_df,orientation=='up same')$prop_promoter_sharing)
#Cliffs delta d= 0.26 small, p<2.2e-16


wilcox.test(filter(plot_df,orientation=='down same')$prop_promoter_sharing, filter(plot_df,orientation=='up opposite')$prop_promoter_sharing)
cliff.delta(filter(plot_df,orientation=='down same')$prop_promoter_sharing, filter(plot_df,orientation=='up opposite')$prop_promoter_sharing)
