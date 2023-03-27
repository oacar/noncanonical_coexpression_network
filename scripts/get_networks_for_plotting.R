#get expanded and canonical networks for plotting
#saves edge list for canonical and expanded network into canonical_network.csv and expanded_network.csv

library(dplyr)
library(RMariaDB)

flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}

getNetworkCutoff<-function(rho_matrix,quant_cutoff){
  x<-flattenCorrMatrix(rho_matrix)
  cutoff<-quantile(x$cor,quant_cutoff,na.rm=TRUE)
  return(cutoff)
}

getAdj<-function(rho_matrix,threshold){
  rho_matrix[rho_matrix <= threshold ]<-0 
  rho_matrix[rho_matrix > threshold]<-1
  diag(rho_matrix)<-0
  return(rho_matrix)
}


#load coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

#load orf info
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

#subset rows and columns to be only canonical genes
idx<-which(colnames(rho) %in% filter(orf_info, is_canonical =='canonical')$transcript)
rho_canonical<-rho[idx,idx]

#get 99.8 percentile cutoff 
whole_network_cutoff<-getNetworkCutoff(rho,0.998)

#save canonical network:
rho_canonical_flat<-flattenCorrMatrix(rho_canonical)
rho_canonical_flat<-filter(rho_canonical_flat,cor>whole_network_cutoff & !is.na(cor))
write.csv(rho_canonical_flat[,c('row','column')],'canonical_network.csv')

#save expanded network:
rho_whole_flat<-flattenCorrMatrix(rho)
rho_whole_flat<-filter(rho_whole_flat, cor>whole_network_cutoff & !is.na(cor))
write.csv(rho_whole_flat[,c('row','column')],'expanded_network.csv')
