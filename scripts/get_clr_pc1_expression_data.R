#get clr + remove top PC1 batch corrected expression data

library(WGCNA)
library(dplyr)
clr<-function(counts)
{ 
  counts<-t(counts)
  # ivar=='clr'
  # use <- ivar2index(counts, ivar)
  logX <- log(counts)
  # logSet <- logX[, use, drop = FALSE]
  ref <- rowMeans(logX,na.rm = T)
  lr <- sweep(logX, 1, ref, "-")
  lr<-t(lr)
  return(lr)
}

#get expression matrices and do batch correction on them
raw_data<-readRDS('raw_counts.RDS')

plot_df<-data.frame('transcript'=rownames(raw_data),
                    'sample_count'=apply(raw_data,1,function(x){length(which(x>5))}),stringsAsFactors = F)

transcripts2keep<-filter(plot_df,sample_count>=400)$transcript


raw_data<-raw_data[transcripts2keep,]
clr_data<-clr(raw_data+0.01)

saveRDS(clr_data,'clr_5_400_noNAs.RDS') #save this so it can be used for applying spqn

ave_exp <- rowMeans(clr_data,na.rm = TRUE)
clr_data <- clr_data - ave_exp # mean centering
clr_data  <- clr_data / matrixStats::rowSds(clr_data,na.rm=TRUE) # variance scaling


clr_pc1<-removePrincipalComponents(as.matrix(clr_data), n = 1)
saveRDS(clr_pc1,'clr_pc1.RDS')


