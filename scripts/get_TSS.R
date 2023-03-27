# get TSS using TIF seq data for all ORFs in coexpression network 

library(dplyr)
library(RMariaDB)
library(data.table)

#load orf info
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


#import tif seq data
tif<-fread('./tif_seq/S1_TIFs.txt',stringsAsFactors = F,fill=TRUE)
tif<-as.data.frame(tif)

idx<-which(tif$strand=='-')
temp<-tif$t5[idx]
tif$t5[idx]<-tif$t3[idx]
tif$t3[idx]<-temp
# tif$tot_count<-tif$ypd+tif$gal

#load coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')

#keep only orfs in coexpression network
orf_info<-filter(orf_info, transcript %in% rownames(rho))


orf_info$TSS<-NA

#for each ORF find all transcripts containing the given ORF, get median transcription start site (TSS)
for(i in 1:nrow(orf_info)){
  tif_sub<-filter(tif,chr==orf_info$chr_num[i] & strand==orf_info$strand[i] & t5<= orf_info$coor1[i] & t3>=orf_info$coor2[i])
  if(nrow(tif_sub)>0)
  {
    if(orf_info$strand[i]=='+')
    {
      orf_info$TSS[i]<-median(tif_sub$t5)
    }else
    {
      orf_info$TSS[i]<-median(tif_sub$t3)
    }
  }
}

saveRDS(orf_info,'orf_info_TSS.RDS')


