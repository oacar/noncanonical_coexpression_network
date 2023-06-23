#script to plot supplementary figure 14: Proportion of de novo ORFs that share a promoter with their neighboring conserved ORF.

library(dplyr)
library(ggpubr)
library(magrittr)
library(forcats)
library(scales)
library(RMariaDB)
library(stringr)
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
tif_info<-dbGetQuery(conn, "SELECT * FROM april.TIF_info")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
dbDisconnect(conn)




rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

orf_pairs<-filter(orf_pairs, distance <=500 & orientation %in% c('downsame','upsame','upopposite') 
                  & primary_orf %in% rownames(rho) & neighbor_orf %in% rownames(rho))



overlap_TSS_info<-readRDS('TSS_overlap_info_prop.RDS')

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


percentage <- plot_df %>%
  group_by(orientation) %>%
  summarise(prop = sum(share_promoter == "yes") / sum(!is.na(share_promoter)),
            se = sqrt(prop * (1 - prop) / sum(!is.na(share_promoter))))


pdf(sprintf('%sprop_share_promoter_counts.pdf',file_path),width=3,height=3)
ggplot(percentage, aes(x = orientation, y = prop, fill=orientation)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "black") +
  labs(x = "orientation", y = "proportion of de novo ORFs that share a promoter with neighboring conserved gene") + 
  scale_fill_manual(values=c('#86B040','#BBB0D7', '#0B6C9A'))+
  theme_minimal()+ theme(legend.position = "none")+ ylim(c(0,1))
dev.off()


