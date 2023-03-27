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

flattenCorrMatrix <- function(cormat){#, pmat) {
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
rho[num_obs<400]<-NA


network_cutoff<-getNetworkCutoff(rho,0.998)
#get adjacency matrix
rho_adj<-getAdj(rho,network_cutoff)

#get degree for each orf
degrees<-rowSums(rho_adj,na.rm=TRUE)
degree_df<-data.frame('transcript'=names(degrees), 'degree'=degrees)
degree_df<-filter(degree_df,degree>0)

#add orf information (ie is orf canonical or noncanonical)
degree_df<-left_join(degree_df,orf_info[,c('transcript','is_canonical')])

pdf("./20230106_figures/degree_density_plot.pdf", width=3,height=3)
ggdensity(degree_df, x = "degree",
          add = "median", 
          color = "is_canonical", fill = "is_canonical",xlab='number of coexpressed partners', 
          ylab='ORF density', legend.title='ORF type',
          palette = c("#7570b3", "#1b9e77"),xscale = "log10")+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)
dev.off()


