library(dplyr)
library(RMariaDB)
library(ggplot2)

conn <- dbConnect(MariaDB(), dbname = 'omer',
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
dbDisconnect(conn)

getOrfsInOneOrientation<-function(pair_df,distance_cutoff){
  #filter orf pairs to only consider those that are within 500 bp of each other
  pairs<-filter(pair_df, distance <= distance_cutoff)
  x<-table(pairs$primary_orf)
  orfsInOneOrientation<-names(which(x==1))
  return(orfsInOneOrientation)
}

#find ORFs that are within 500 bp of a conserved gene in only one orientation 
orfs_500<-getOrfsInOneOrientation(orf_pairs,distance_cutoff = 500)
saveRDS(orfs_500,'./OrfsInOneOrientation_500bp.RDS')

