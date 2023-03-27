#script to plot the number of cORFs that are coexpressed with atleast 1 nORF and number of cORFs that arent coexpressed with any nORFs

library(dplyr)
library(RMariaDB)
library(ggplot2)
library(ggpubr)
library(EnvStats)
significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

flattenCorrMatrix <- function(cormat){
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

getNetworkCutoff<-function(rho_matrix,quant_cutoff){
  x<-flattenCorrMatrix(rho_matrix)
  cutoff<-quantile(x$cor,quant_cutoff,na.rm=TRUE)
  return(cutoff)
}

getAdj<-function(rho_matrix,threshold){
  rho_matrix[is.na(rho_matrix)]<-0
  rho_matrix[rho_matrix <= threshold ]<-0 
  rho_matrix[rho_matrix > threshold]<-1
  diag(rho_matrix)<-0
  return(rho_matrix)
}

#load orf info
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

#import coexpression matrix
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
#keep only pairs where there is atleast 400 samples expressing both ORFs
rho[num_obs<400]<-NA

#get rho cutoff corresponding to 99.8th percentile
network_cutoff<-getNetworkCutoff(rho,0.998)
#use rho>99.8th percentile cutoff to create adjacency matrix
rho_adj<-getAdj(rho,network_cutoff)

#get list of canonical ORFs
cORFs<-filter(orf_info, is_canonical=='canonical' & transcript %in% rownames(rho))$transcript
#get list of noncanonical ORFs
nORFs<-filter(orf_info, is_canonical=='noncanonical' & transcript %in% rownames(rho))$transcript

#for all canonical ORFs, get their noncanonical degree (ie how many noncanonical ORFs they are coexpressed with)
cORFs_nDegree<-rowSums(rho_adj[cORFs,nORFs],na.rm = TRUE)

#for all canonical ORFs, get their canonical degree (ie how many canonical ORFs they are coexpressed with)
cORFs_cDegree<-rowSums(rho_adj[cORFs,cORFs],na.rm=TRUE)
all(names(cORFs_cDegree)==names(cORFs_nDegree))
plot_df<-data.frame('transcript'=names(cORFs_cDegree), 'canonical_degree'=cORFs_cDegree, 'noncanonical_degree'=cORFs_nDegree)

#keep only ORFs that are in the network (ie are coexpressed with at least one other ORF at rho>99.8th percentile)
plot_df<-filter(plot_df,(canonical_degree>0 | noncanonical_degree>0))

#how many cORFs are coexpressed with at least one nORF? how many are coexpressed with no nORFs 
filter(plot_df, noncanonical_degree==0) %>% nrow()
filter(plot_df, noncanonical_degree>0) %>% nrow()


pdf("./cORFsCoexpWithnORFs.pdf", width=2.25,height=3)
plot_df %>% mutate(noncanonical_degree =case_when(noncanonical_degree==0 ~'noncanonical\ndegree = 0 ',
                                                  noncanonical_degree >0 ~'noncanonical\ndegree > 0')) %>%
  count(noncanonical_degree) %>%
  ggbarplot(., x="noncanonical_degree",y="n",ylab ='cORF count',
            fill='#7570b3', color=NA,
            label = FALSE)+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + 
  grids(axis = c("xy"), color = "grey92", size = NULL, linetype = NULL)+
  rremove("legend")+
  rremove("x.title")

dev.off()

