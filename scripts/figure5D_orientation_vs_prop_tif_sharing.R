#this script plots figure 5D: where x axis orientation of de novo ORF relative to conserved gene 
#y axis is proportion of time a de novo ORF shares a transcript with neighboring conserved gene as determined by TIF seq data from Pelechano et al 2013.
#ie y axis is number of transcripts containing both de novo ORF and conserved gene / number of tifs containing de novo ORF 


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

pdf(sprintf('%sshareTIFProp_Orientation.pdf',file_path),width=2.2,height=3)
tif_info %>%  mutate(denovo_tif_ratio_counts2= case_when(tif_intersect_counts==0 ~ 0, TRUE ~ denovo_tif_ratio_counts)) %>% 
  ggboxplot(x="orientation",y="denovo_tif_ratio_counts2", notch=TRUE,ylab="proportion of de novo ORF\ntranscripts that also contain gene") +
  stat_n_text(size = 2.75)+
  stat_compare_means(label = "p.signif", label.x = 1.45, size=significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("legend")
dev.off()

tif_info<-tif_info %>%  mutate(denovo_tif_ratio_counts2= case_when(tif_intersect_counts==0 ~ 0, TRUE ~ denovo_tif_ratio_counts))
wilcox.test(filter(tif_info,orientation=='downsame')$denovo_tif_ratio_counts2, filter(tif_info,orientation=='upsame')$denovo_tif_ratio_counts2)
cliff.delta(filter(tif_info,orientation=='downsame')$denovo_tif_ratio_counts2, filter(tif_info,orientation=='upsame')$denovo_tif_ratio_counts2)
#Cliffs delta d= 0.26 small, p<2.2e-16