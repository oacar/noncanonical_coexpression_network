#get median TPM across all samples for each ORF
#only include ORF which have atleast 400 samples that the orf is detected in
#defining detected as raw count >5

library(dplyr)
library(RMariaDB)

getTPMdf<-function(raw_cut, sample_cut, NA_cut, raw_data,coexpression_orf_list){
  plot_df<-data.frame('transcript'=rownames(raw_data),
                      'sample_count'=apply(raw_data,1,function(x){length(which(x>raw_cut))}),stringsAsFactors = F)
  
  transcripts2keep<-filter(plot_df,sample_count>=sample_cut)$transcript
  
  # coexpression_orf_list<-filter(orf_data, transcript %in% transcripts2keep)
  
  raw_data<-raw_data[transcripts2keep,]
  
  #convert 0s to 0
  raw_data[raw_data<NA_cut]<-0
  
  
  length_df<-data.frame('transcript'=rownames(raw_data))
  length_df<-left_join(length_df,coexpression_orf_list[,c('transcript','length')])
  lengths_matrix<-matrix(length_df$length, nrow=length(length_df$length), ncol=ncol(raw_data), byrow=FALSE)
  lengths_matrix<-lengths_matrix/1000
  
  #calculate TPM
  rpk<-raw_data/lengths_matrix
  scaling<-colSums(rpk)/1E6
  scaling_matrix<-matrix(scaling, nrow=nrow(raw_data),ncol=length(scaling), byrow=TRUE)
  tpm<-rpk/scaling_matrix
  tpm[raw_data==0]<-NA
  median_tpm<-apply(tpm,1,median,na.rm=TRUE)
  
  
  mean_tpm<-rowMeans(tpm,na.rm=TRUE)
  median_tpm<-data.frame('transcript'=names(median_tpm),'tpm'=median_tpm,'sample_count'=apply(tpm,1,function(x){length(which(!is.na(x)))}))
  median_tpm<-left_join(median_tpm,coexpression_orf_list[,c('transcript','orf_id','classification')])
  median_tpm<-median_tpm[,c('orf_id','transcript','sample_count','tpm','classification')]
  
  return(median_tpm)
}

raw<-readRDS('./raw_counts.RDS')
raw[raw<5]<-0


#get ORF data
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_list<-dbGetQuery(conn, "SELECT * FROM april.expression_transcriptome")
dbDisconnect(conn)

orf_list$length<-orf_list$coor2-orf_list$coor1


tpm_5_400<-getTPMdf(raw_cut=5, sample_cut=400, NA_cut=5, raw_data=raw, coexpression_orf_list=orf_list)



