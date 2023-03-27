#for each de novo ORF that is within 500 bp of a conserved gene, either upstream or downstream on the same strand
#find out how many times that ORF shares a transcript with neighboring gene

library(dplyr)
library(data.table)
library(ggplot2)
library(RMariaDB)

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
dbDisconnect(conn)


orf_pairs<-filter(orf_pairs, orientation %in% c('downsame', 'upsame') & distance<=5000)

#import tif seq data (from Pelechano et al 2013)
tif<-fread('S1_TIFs.txt',stringsAsFactors = F,fill=TRUE)
tif<-as.data.frame(tif)
idx<-which(tif$strand=='-')
temp<-tif$t5[idx]
tif$t5[idx]<-tif$t3[idx]
tif$t3[idx]<-temp


orf_pairs$denovo_tif_total_counts<-NA #number of TIFs containing de novo ORF
orf_pairs$gene_tif_total_counts<-NA #number of TIFs containing conserved gene
orf_pairs$tif_intersect_counts<-NA #number of TIFs containing both de novo ORF and neighboring conserved genes
orf_pairs$denovo_tif_ratio_counts<-NA # number of tifs containing both de novo ORF & neighboring gene / number of tifs containing de novo ORF 
orf_pairs$gene_tif_ratio_counts<-NA # number of tifs containing both de novo ORF & neighboring gene / number of tifs containing conserved gene 

for(i in 1:nrow(orf_pairs)){
  
  orf<-filter(orf_info, transcript==orf_pairs$primary_orf[i])
  gene<-filter(orf_info, transcript==orf_pairs$neighbor_orf[i])
  
  orf_tifs<-filter(tif,chr==orf$chr_num & strand==orf$strand & t5<= orf$coor1 & t3>=orf$coor2)
  gene_tifs<-filter(tif,chr==gene$chr_num & strand==gene$strand & t5<= gene$coor1 & t3>=gene$coor2)
  
  
  orf_pairs$denovo_tif_total_counts[i]<-sum(orf_tifs$ypd)+sum(orf_tifs$gal)
  orf_pairs$gene_tif_total_counts[i]<-sum(gene_tifs$ypd)+sum(gene_tifs$gal)
  
  orf_gene_tif<-intersect(orf_tifs,gene_tifs)
  orf_pairs$tif_intersect_counts[i]<-sum(orf_gene_tif$ypd)+sum(orf_gene_tif$gal)
  
  orf_pairs$denovo_tif_ratio_counts[i]<-orf_pairs$tif_intersect_counts[i]/orf_pairs$denovo_tif_total_counts[i]
  orf_pairs$gene_tif_ratio_counts[i]<-orf_pairs$tif_intersect_counts[i]/orf_pairs$gene_tif_total_counts[i]
}


saveRDS(orf_pairs,'TIFinfo.RDS')

