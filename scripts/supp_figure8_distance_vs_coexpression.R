#this script plots supplementary figure 8: looking at association between closest distance between a denovo orf and its neighboring conserved genes
#in all orientations 
#x axis is closest distance between a denovo ORF and neighboring conserved gene
#y axis is coexpression of de novo ORF and neighboring conserved gene
#graph is faceted by orientation 

library(RMariaDB)
library(dplyr)
library(ggplot2)
library(ggpubr)


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
orf_info<-dbGetQuery(conn,"SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

#load coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

orf_pairs<-filter(orf_pairs, primary_orf %in% rownames(rho))
orf_pairs<-filter(orf_pairs, neighbor_orf %in% rownames(rho))
orf_pairs$rho<-NA
for(i in 1:nrow(orf_pairs)){
  orf_pairs$rho[i]<-rho[orf_pairs$primary_orf[i],orf_pairs$neighbor_orf[i]]
}
file_path<-'./20230106_figures/'



pdf(sprintf("%scoexpression_distance_lineplot.pdf", file_path), width = 6, height = 4.5)
orf_pairs %>% filter(orientation %in% c('downsame','upsame','upopposite','downopposite') & !is.na(rho) & distance <10000) %>%
  mutate(orientation= case_when(orientation=='downsame' ~'down same', 
                                orientation=='upsame' ~'up same',
                                orientation=='upopposite' ~'up opposite',
                                orientation=='downopposite' ~'down opposite')) %>% 
  # group_by(orientation) %>%
  # mutate(orientation = paste0( orientation, "\nn=", n())) %>%
  ggplot(aes(x=distance, y=rho))+geom_point(alpha=0.5,color = "black", shape = 21, size = 1) + facet_wrap(~orientation)+ theme_minimal()+
  theme(axis.text=element_text(size=axis_text_size))+
  geom_smooth(method='lm')+ #adds line with error region
  ggpubr::stat_cor(label.x =4000, label.y=-0.3,method="spearman")
dev.off()

