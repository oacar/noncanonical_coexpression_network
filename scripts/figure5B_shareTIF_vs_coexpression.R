#this script plots figure 5B: where x axis is whether an ORF shares at least one transcript with neighboring
# conserved gene as determined by TIF seq data from Pelechano et al 2013.
#y axis is coexpression of de novo ORF with neighboring conserved gene 
#graph is faceted by orientation ie either for de novo ORFs that are upstream of a conserved gene or downstream of a conserved gene

#tif_info data generated in get_tif_info.R

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
OrfsInOneOrientation<-dbGetQuery(conn, "SELECT * FROM april.OrfsInOneOrientation_500bp")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
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



#path to save plots
file_path<-'./20230106_figures/'

wilcox.test(filter(tif_info,orientation=='downsame' & tif_intersect_counts>0)$rho, filter(tif_info, orientation =='downsame' & tif_intersect_counts==0)$rho)
cliff.delta(filter(tif_info,orientation=='downsame' & tif_intersect_counts>0)$rho, filter(tif_info, orientation =='downsame' & tif_intersect_counts==0)$rho)
#down same: cliffs delta d=0.283 (small) p=2.994e-09

wilcox.test(filter(tif_info,orientation=='upsame' & tif_intersect_counts>0)$rho, filter(tif_info, orientation =='upsame' & tif_intersect_counts==0)$rho)
cliff.delta(filter(tif_info,orientation=='upsame' & tif_intersect_counts>0)$rho, filter(tif_info, orientation =='upsame' & tif_intersect_counts==0)$rho)
#upsame: cliffs delta d=0.309 (small) p< 2.2e-16

#background is the coexpression values for 1 million random de novo-conserved ORF pairs that are on separate chromosomes
background<-readRDS('rhoForDCPairsNotOnSameChr.RDS')
colnames(background)[3]<-'rho'
background$orientation<-'background'

p<-tif_info %>% filter(orientation %in% c('downsame','upsame')) %>%
  mutate(share_tif= case_when( tif_intersect_counts>0 ~'yes',
                               tif_intersect_counts==0 ~'no')) %>%
  mutate(orientation= case_when(orientation =='downsame' ~'down same',
                                orientation =='upsame' ~'up same')) %>% 
  ggboxplot(x="share_tif",y="rho",xlab="share > 1 transcript",ylab ="coexpression with neighbor",
            fill="orientation",notch=TRUE, palette=c('#86B040','#BBB0D7'),facet.by="orientation") +
  stat_n_text(size = 2.75)+
  geom_hline(yintercept=median(background$rho),linetype='dotted')+
  stat_compare_means(aes(group = share_tif),label = "p.signif", label.x = 1.45, size=significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)+
  rremove("legend")

pdf(sprintf('%sCoxpressionVsshareTIF_boxplot_1Tif_500.pdf',file_path),width=2.75,height=3.2)
ggpar(p)
dev.off()

