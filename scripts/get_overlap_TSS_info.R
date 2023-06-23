#script to get information about if the TSS of up opposite ORFs overlap with TSS of neighboring cORF
#uses tif seq data from Pelechano et al 2013.

library(dplyr)
library(RMariaDB)
library(data.table)

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")

dbDisconnect(conn)

rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

orf_pairs<-filter(orf_pairs,orientation=='upopposite' & distance <1000)

#import tif seq data (from Pelechano et al 2013)
tif<-fread('./tif_seq/S1_TIFs.txt',stringsAsFactors = F,fill=TRUE)
tif<-as.data.frame(tif)
idx<-which(tif$strand=='-')
temp<-tif$t5[idx]
tif$t5[idx]<-tif$t3[idx]
tif$t3[idx]<-temp


#get TSS for both de novo ORF and neighboring gene 

getTSS<-function(tifs,strand){
  if(strand == '+'){
    TSSs<-c(rep(tifs$t5,tifs$ypd),rep(tifs$t5,tifs$gal))
  } else{ #ie negative strand
    TSSs<-c(rep(tifs$t3,tifs$ypd),rep(tifs$t3,tifs$gal))
  }
  return(TSSs)
}

orf_pairs$num_time_share_promoter<-NA
orf_pairs$prop_time_share_promoter<-NA
orf_pairs$num_denovo_tifs<-NA
orf_pairs$num_gene_tifs<-NA
for(i in 1:nrow(orf_pairs)){
  orf<-filter(orf_info, transcript==orf_pairs$primary_orf[i])
  gene<-filter(orf_info, transcript==orf_pairs$neighbor_orf[i])
  orf_tifs<-filter(tif,chr==orf$chr_num & strand==orf$strand & t5<= orf$coor1 & t3>=orf$coor2)
  gene_tifs<-filter(tif,chr==gene$chr_num & strand==gene$strand & t5<= gene$coor1 & t3>=gene$coor2)
  if(nrow(gene_tifs)>0){
    gene_tss<-getTSS(gene_tifs,strand=gene$strand)
    orf_pairs$num_gene_tifs[i]<-length(gene_tss)
  }
  if(nrow(orf_tifs)>0){
    orf_tss<-getTSS(orf_tifs,strand=orf$strand)
    orf_pairs$num_denovo_tifs[i]<-length(orf_tss)
  }
  
  
  if(nrow(orf_tifs)>0 & nrow(gene_tifs)>0){
    if(gene$strand =='+'){
      orf_pairs$num_time_share_promoter[i]<-length(which(orf_tss < median(gene_tss)))
     
    }else{ #negative strand
      orf_pairs$num_time_share_promoter[i]<-length(which(orf_tss > median(gene_tss)))
    }
    orf_pairs$prop_time_share_promoter[i]<-orf_pairs$num_time_share_promoter[i]/orf_pairs$num_denovo_tifs[i]
  } else{orf_pairs$prop_time_share_promoter[i]<-0}
  
}

saveRDS(orf_pairs,'TSS_overlap_info_prop.RDS')