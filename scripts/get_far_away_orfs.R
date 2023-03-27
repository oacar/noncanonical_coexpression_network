#identify de novo orfs that are not near any conserved genes (ie orfs that are further than 500 bp away from a gene in any orientation)
library(dplyr)
library(RMariaDB)
library(ggplot2)

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
dbDisconnect(conn)

getFarAwayOrfs<-function(pair_df,distance_cutoff){
  x<-data.frame('transcript'=unique(pair_df$primary_orf),'far_away'=NA)
  for(i in 1:nrow(x)){
    sub<-filter(pair_df, primary_orf==x$transcript[i])
    if(all(sub$distance > distance_cutoff)){
      x$far_away[i]<-'yes'
    }else{x$far_away[i]<-'no'}
  }
  x<-filter(x,far_away=='yes')
  return(x$transcript)
}


nearPol13<-readRDS('orfs_near_pol3_pol1_targets.RDS')

orfs_500<-getFarAwayOrfs(orf_pairs,distance_cutoff = 500)
orfs_500<-setdiff(orfs_500,nearPol13$primary_orf)
saveRDS(orfs_500,'OrfsFarAway_500bp.RDS')
