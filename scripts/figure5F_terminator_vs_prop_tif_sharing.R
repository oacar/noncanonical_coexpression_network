#this script plots figure 5F: where x axis is whether an downsame de novo ORF has a transcription terminator present between it 
# and upstream neighboring conserved gene. terminator binding determined from chip exo data of nrd1 and pcf11 binding from rossi et al 2021
# conserved gene as determined by TIF seq data from Pelechano et al 2013.
#y axis is proportion of time down same de novo ORF shares a transcript with neighboring conserved gene, determined by TIF seq from Pelechano et al 2013.


#tif_info data generated in get_tif_info.R
#terminator_info data generated in get_terminator_info.R


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
terminator_info<-dbGetQuery(conn, "SELECT * FROM april.Terminator_info")
tif_info<-dbGetQuery(conn, "SELECT * FROM april.TIF_info")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


terminator_info<-mutate(terminator_info, terminator_present=case_when(nrd1>=1 | pcf11>=1 ~'yes',
                                                                      nrd1==0 & pcf11==0 ~'no'))

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

tif_info<-left_join(tif_info, terminator_info[,c('primary_orf','neighbor_orf','nrd1','pcf11','terminator_present')])

#path to save plots
file_path<-'./20230106_figures/'

p<-tif_info %>% filter(orientation %in% c("downsame")) %>%
  mutate(terminator_present=factor(terminator_present,levels=c('yes','no'))) %>% 
  mutate(denovo_tif_ratio_counts=case_when(tif_intersect_counts==0 ~ 0,
                                           tif_intersect_counts >0 ~ denovo_tif_ratio_counts)) %>%
  ggboxplot(x="terminator_present",y="denovo_tif_ratio_counts",ylab ="proportion of de novo\ncontaining transcripts shared", xlab="terminator present",
            notch=TRUE, fill= "terminator_present",
            palette=c('#FF4236','#C03830'))+
  stat_compare_means(aes(group = terminator_present), label = "p.signif", label.x = 1.45, size=significance_size)+
  stat_n_text(size = 2.75)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  rremove("legend")
pdf(sprintf('%sterminatorVsshareTIF_boxplot_500.pdf',file_path),width=2,height=3)
ggpar(p,legend.title = "terminator present",font.legend = axis_text_size)
dev.off()


x<-tif_info %>% filter(orientation %in% c("downsame")) %>%
  mutate(terminator_present=factor(terminator_present,levels=c('yes','no')))%>%
  mutate(denovo_tif_ratio_counts=case_when(tif_intersect_counts==0 ~ 0,
                                           tif_intersect_counts >0 ~ denovo_tif_ratio_counts))

wilcox.test(filter(x,terminator_present=='yes')$denovo_tif_ratio_counts, filter(x, terminator_present=='no')$denovo_tif_ratio_counts)
cliff.delta(filter(x,terminator_present=='yes')$denovo_tif_ratio_counts, filter(x, terminator_present=='no')$denovo_tif_ratio_counts)
#medium (Cliff’s Delta d=-0.386, Mann-Whitney U-test, p=1.591e-10)


