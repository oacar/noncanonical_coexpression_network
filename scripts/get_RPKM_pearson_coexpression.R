#get RPKM + pearson correlation coexpression network


library(dplyr)
library(RMariaDB)

getRpkmMatrix<-function(raw_cut, sample_cut, NA_cut, raw_data,coexpression_orf_list){
  plot_df<-data.frame('transcript'=rownames(raw_data),
                      'sample_count'=apply(raw_data,1,function(x){length(which(x>raw_cut))}),stringsAsFactors = F)
  
  transcripts2keep<-filter(plot_df,sample_count>=sample_cut)$transcript
  
  raw_data<-raw_data[transcripts2keep,]
  
  #convert 0s to 0
  raw_data[raw_data<NA_cut]<-0
  
  
  length_df<-data.frame('transcript'=rownames(raw_data))
  length_df<-left_join(length_df,coexpression_orf_list[,c('transcript','length')])
  lengths_matrix<-matrix(length_df$length, nrow=length(length_df$length), ncol=ncol(raw_data), byrow=FALSE)
  lengths_matrix<-lengths_matrix/1000
  
  # rpk<-raw_data/lengths_matrix
  scale_factor<-colSums(raw_data)/1E6
  scaling_matrix<-matrix(scale_factor, nrow=nrow(raw_data),ncol=length(scale_factor), byrow=TRUE)
  rpkm<-raw_data/scaling_matrix
  rpkm<-rpkm/lengths_matrix
  rpkm[raw_data==0]<-NA
  
  return(rpkm)
}

raw<-readRDS('raw_counts.RDS')
raw[raw<5]<-0


#get ORF data

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_list<-dbGetQuery(conn, "SELECT * FROM april.expression_transcriptome")
dbDisconnect(conn)

orf_list$length<-orf_list$coor2-orf_list$coor1

rpkm_5_400_matrix<-getRpkmMatrix(raw_cut=5, sample_cut=400, NA_cut=5, raw_data=raw, coexpression_orf_list=orf_list)
saveRDS(rpkm_5_400_matrix,'rpkm_5_400_matrix.RDS')


coexpression_pearson_rpkm<-cor(t(rpkm_5_400_matrix), use="pairwise.complete.obs")

#import the number of samples where both ORFs in a pair were detected
numObs<-readRDS('numobs_raw5_sample400.RDS')
all(rownames(coexpression_pearson_rpkm)==rownames(numObs))
all(colnames(coexpression_pearson_rpkm)==colnames(numObs))

#if less than 400 samples where both ORFs are detected, set coexp value to NA
coexpression_pearson_rpkm[numObs<400]<-NA

saveRDS(coexpression_pearson_rpkm,'rpkm_raw5_sample400_coexp_pearson.RDS')


