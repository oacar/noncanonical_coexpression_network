
#script to plot figure 4F: x axis is wether ORFs can piggyback or cannot piggyback 
#y axis is coexpression of de novo ORF with neighboring conserved gene 

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
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
dbDisconnect(conn)



#get background distribution ie all other orf pairs not on the same chromosome 
rho_matrix<-readRDS('spqn_raw5_sample400.RDS')
numObs<-readRDS('numobs_raw5_sample400.RDS')
flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}
rho_matrix[numObs<400]<-NA
rho<-flattenCorrMatrix(rho_matrix)
rho<-filter(rho,!is.na(cor))
rho<-left_join(rho,orf_info[,c('transcript','chr_num','classification')],by=c('row'='transcript'))
rho<-left_join(rho,orf_info[,c('transcript','chr_num','classification')],by=c('column'='transcript'))

background<-filter(rho, chr_num.x != chr_num.y & ( (classification.x=='denovo' & classification.y=='conserved') | (classification.x=='conserved' & classification.y=='denovo')))
saveRDS(background, 'rhoForDCPairsNotOnSameChr.RDS')
colnames(background)[3]<-'rho'
background$orientation<-'background'


background<-background[,c('row','orientation','column','rho')]
colnames(background)<-c('primary_orf','orientation','neighbor_orf','rho') #do this so u can rbind with orf_pairs to make plot

#subset background pairs bc otherwise 34 million datapoints... 
background_sub<-background[sample(1:nrow(background),size=1E6, replace =F),]
quantile(background_sub$rho) 
# > quantile(background_sub$rho)
# 0%        25%        50%        75%       100%
# -0.4126164  0.3210812  0.4753267  0.6134138  0.9607631

#^^ which is the same basically as the full 34 mill observations
# > quantile(background$rho)
# 0%        25%        50%        75%       100%
# -0.4126164  0.3214636  0.4754028  0.6135944  0.9607631




orf_pairs<-filter(orf_pairs, distance <=500)
orf_pairs$rho<-NA
orf_pairs<-filter(orf_pairs, primary_orf %in% rownames(rho_matrix))
orf_pairs<-filter(orf_pairs, neighbor_orf %in% rownames(rho_matrix))
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho_matrix[orf_pairs$neighbor_orf[i],orf_pairs$primary_orf[i]]
}

orf_pairs<-filter(orf_pairs, orientation !='sense-overlap')
orf_pairs<-mutate(orf_pairs, orientation=case_when(orientation=='downsame' ~'down\nsame',
                                                   orientation=='upsame' ~'up\nsame',
                                                   orientation =='antisense-overlap' ~'antisense\noverlap',
                                                   orientation=='downopposite' ~'down\nopposite',
                                                   orientation=='upopposite' ~'up\nopposite'))



#can piggyback vs cannot piggyback:
plot_df<-orf_pairs %>% mutate(orientation_type=case_when(orientation %in% c('down\nsame','up\nsame','up\nopposite') ~ 'can\npiggyback',
                                                         orientation %in% c('down\nopposite','antisense\noverlap') ~ 'cannot\npiggyback'))
p<- plot_df%>%
  ggboxplot(x = "orientation_type", y = "rho", ylab = "coexpression with neighbor", notch = TRUE, fill = "orientation_type",
            palette = c('#F6DFE0','#DBE8C6')) +
  geom_hline(yintercept=median(background_sub$rho),linetype='dotted')+
  stat_n_text(size = 2.75) +
  stat_compare_means(label = "p.signif", label.y = 1, size = significance_size)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  # rotate_x_text(45)+
  rremove("xlab")+
  rremove("legend")

pdf(sprintf('%sCoexpressionVsOrientation_type.pdf',file_path),width=2.25,height=2.5)
ggpar(p)
dev.off()

wilcox.test(filter(plot_df,orientation_type =='can\npiggyback')$rho, filter(plot_df,orientation_type =='cannot\npiggyback')$rho)
cliff.delta(filter(plot_df,orientation_type =='can\npiggyback')$rho, filter(plot_df,orientation_type =='cannot\npiggyback')$rho)
#Cliffs delta: 0.426 medium, p<2.2e-16

