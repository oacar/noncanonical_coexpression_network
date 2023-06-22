#script to plot supp figure 13A: x axis orientation of de novo ORF relative to conserved gene,
#y axis coexpression with neighbor

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
orientation_levels <- c("background", "down\nopposite", "down\nsame", "antisense\noverlap", "up\nopposite", "up\nsame")
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
dbDisconnect(conn)

flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}

#get background distribution ie all other denovo-conserved orf pairs not on the same chromosome 
rho_matrix<-readRDS('spqn_raw5_sample400.RDS')
numObs<-readRDS('numobs_raw5_sample400.RDS')

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
colnames(background)<-c('primary_orf','orientation','neighbor_orf','rho') 

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


figure_data<-rbind(orf_pairs[,c('primary_orf','orientation','neighbor_orf','rho')],background_sub) 
saveRDS(figure_data,'coexpression_data_background_dcpairs.RDS')

fig_data<-figure_data
stat.test <- data.frame(
  group1 = c('background','background', "background",'background', "background"),
  group2= c('down\nopposite','down\nsame', "antisense\noverlap",'up\nopposite', "up\nsame")
)
stat.test$effect.size<- sapply(stat.test$group2, function(x) {
  fig_data %>%
    filter(orientation==x | orientation == "background") %>%
    cliff.delta(rho ~ orientation, data = .)%$%magnitude
})
effect.size.labels = c(" "="negligible", "+"="small", "++"="medium", "+++"="large" )
stat.test$label = names(effect.size.labels)[match(stat.test$effect.size,effect.size.labels)]
p<-fig_data %>% mutate(orientation = factor(str_replace(orientation, "\n ", "\n"), levels = orientation_levels)) %>%
  ggboxplot(x = "orientation", y = "rho", ylab = "coexpression", notch = TRUE, fill = "#BBCC32") +
  geom_hline(yintercept=median(filter(fig_data,orientation=='background')$rho),linetype='dotted')+
  stat_n_text(size = 2.75) +
  stat_compare_means(label = "p.signif", ref.group = "background", label.y = 1, size = significance_size) +
  theme(axis.text = element_text(size = axis_text_size), axis.title.y = element_text(size = axis_title_size)) +
  rremove("xlab")+
  geom_text(data = stat.test%>%
              mutate(group2= factor(str_replace(group2, "\n", "\n"), levels = orientation_levels)) , aes(x = group2, y = 1.15, label = label), size = significance_size) 
ggsave(filename = "./20230106_figures/orientationVsCoexp.pdf", plot = p, width = 3.6, height = 3)




#down same
o<-'down\nsame'
wilcox.test(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
cliff.delta(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
#down same: large (Cliff’s Delta d=0.644, Mann-Whitney U-test, p< 2.2e-16)

o<-'up\nsame'
wilcox.test(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
cliff.delta(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
#up same: small (Cliff’s Delta d=0.318, Mann-Whitney U-test, p< 2.2e-16)

o<-'antisense\noverlap'
wilcox.test(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
cliff.delta(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
#antisense overlap: negligible (Cliff’s Delta d=0.082, Mann-Whitney U-test, p=0.0013)

o<-'up\nopposite'
wilcox.test(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
cliff.delta(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
#up opposite:small (Cliff’s Delta d=0.281, Mann-Whitney U-test, p< 2.2e-16))

o<-'down\nopposite'
wilcox.test(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
cliff.delta(filter(fig_data,orientation==o)$rho, filter(fig_data,orientation=='background')$rho)
#down opposite: negligible (Cliff’s Delta d=-0.0222, Mann-Whitney U-test, p=0.305)




