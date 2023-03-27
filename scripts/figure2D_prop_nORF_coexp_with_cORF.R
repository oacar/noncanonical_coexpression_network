#script to look at the proportion of noncanonical ORFs (nORFs) that are coexpressed with at least one canonical ORF (cORF) in the expanded (ie real network) compared to 
#1000 randomized networks, where the network randomization was done in rewire_network.py
#edge list for expanded network generated in file get_networks_for_plotting.R

library(reshape2)
library(dplyr)
library(RMariaDB)
library(ggplot2)
library(ggpubr)
library(EnvStats)
library(rstatix)
library(tibble)

significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

#import network
expanded_network<-read.csv('expanded_network.csv')

nORFs_connected_to_cORFs<-function(network, orf_info_df){
  cORFs<-filter(orf_info_df, is_canonical=='canonical')$transcript
  nORFs<-filter(orf_info_df, is_canonical=='noncanonical')$transcript
  
  return(unique(c(filter(network, row %in% nORFs & column %in% cORFs)$row,filter(network, row %in% cORFs & column %in% nORFs)$column)))
}

nORFs_not_connect_to_cORFs<-function(network, orf_info_df,nORFsConnected2cORFs){
  
  nORFs<-filter(orf_info_df, is_canonical=='noncanonical' & transcript %in% c(network$row, network$column))$transcript
  
  return(setdiff(nORFs,nORFsConnected2cORFs))
  
}


#for each random network, get the number of nORFs coexpressed with atleast 1 cORF and get the number of nORFs that are not coexpressed with any cORF
num_networks<-length(list.files(path = "./randomized_networks/"))

random_df<-data.frame('canonical_degree_zero'=NA, 'canonical_degree_nonzero'=NA,'network_type'=rep('random',num_networks))
for(i in 1:num_networks){
  random_network<-read.csv(sprintf('./randomized_networks/%i.csv',i-1),header=FALSE)
  colnames(random_network)[2:3]<-c('row','column')
  
  random_nORFs_with_cORFs<-nORFs_connected_to_cORFs(random_network,orf_info)
  random_nORFs_without_cORFs<-nORFs_not_connect_to_cORFs(random_network,orf_info,random_nORFs_with_cORFs)
  
  random_df$canonical_degree_zero[i]<-length(random_nORFs_without_cORFs)
  random_df$canonical_degree_nonzero[i]<-length(random_nORFs_with_cORFs)
  
}

saveRDS(random_df,'nORFsWithcORFs_random.RDS')

#calculate proportion of nORFs not coexpressed with any cORFs
random_df$prop_no_cORF<-random_df$canonical_degree_zero/(random_df$canonical_degree_zero+random_df$canonical_degree_nonzero)

#calculate proportion of nORFs coexpressed with at least one cORF
random_df$prop_cORF<-random_df$canonical_degree_nonzero/(random_df$canonical_degree_zero+random_df$canonical_degree_nonzero)

#get mean proportion of nORFs coexpressed with atleast one cORF over the 1000 randomzied networks
mean_prop<-mean(random_df$prop_cORF)

prop_upper_ci = mean_prop + sd(random_df$prop_cORF)/sqrt(length(random_df$prop_cORF))
prop_lower_ci = mean_prop - sd(random_df$prop_cORF)/sqrt(length(random_df$prop_cORF))

mean(random_df$canonical_degree_zero)
mean(random_df$canonical_degree_nonzero)


expanded_nORFs_with_cORFs<-nORFs_connected_to_cORFs(expanded_network,orf_info)
expanded_nORFs_without_cORFs<-nORFs_not_connect_to_cORFs(expanded_network,orf_info,expanded_nORFs_with_cORFs)


cm<-matrix(nrow=2,ncol=2)
rownames(cm)<-c('nORFsWithoutCORFs','nORFsWithCORFs')
colnames(cm)<-c('random','real')
cm['nORFsWithCORFs','real']<-length(expanded_nORFs_with_cORFs)
cm['nORFsWithCORFs','random']<-round(mean(random_df$canonical_degree_nonzero))
cm['nORFsWithoutCORFs','real']<-length(expanded_nORFs_without_cORFs)
cm['nORFsWithoutCORFs','random']<-round(mean(random_df$canonical_degree_zero))

getPlotdf<-function(cm){
  real_prop = 100*(cm[2,2]/(cm[1,2]+cm[2,2]))
  random_prop = 100*(cm[2,1]/(cm[1,1]+cm[2,1]))
  plot_df<-data.frame('networkType'=c('real','random'),'prop'=c(real_prop, random_prop))
  return(plot_df)
}

#convert contingency matrix into data frame that contains the proportion of nORFs in real network that are coexpressed with cORFs and 
#the proportion of nORFs in randomized networks that are coexpressed with cORFs 
plotdf<-getPlotdf(cm)

#calculate CI for randomized networks
plotdf$upper_ci<-NA
plotdf$lower_ci<-NA
plotdf$upper_ci[(which(plotdf$networkType=='real'))]<-NA
plotdf$lower_ci[(which(plotdf$networkType=='real'))]<-NA
plotdf$upper_ci[(which(plotdf$networkType=='random'))]<-prop_upper_ci*100
plotdf$lower_ci[(which(plotdf$networkType=='random'))]<-prop_lower_ci*100

plotdf<- plotdf %>% mutate(networkType=case_when(networkType=='real' ~'expanded\nnetwork\nn=1',
                                                 networkType=='random'~'random\nnetworks\nn=1000'))

complex_stats<-pairwise_fisher_test(cm) %>% add_column('networkType'='expanded\nnetwork\nn=1')

pdf("./20230106_figures/nORFsCoexpWithcORFs.pdf", width=2,height=3)
plotdf %>%
  ggbarplot( x="networkType", y="prop", xlab="cORF pair type",ylab ='% of nORFs coexpressed with cORFs',
             fill = "networkType", color=NA, palette = c('#F26322','#FCAF17'),
             label = FALSE)+
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci), width = 0.25)+
  stat_pvalue_manual(complex_stats, x='networkType', label="p.adj.signif",y.position=89.9, size = significance_size)+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + grids(axis = c("xy"), color = "grey92", size = NULL, linetype = NULL)+
  rremove("x.title")+
  rremove("legend")
dev.off()


