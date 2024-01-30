#generate coexpression matrix that includes lowly expressed rna seq observations (ie does not set raw counts <5 to NA) 
#then calculate coexpression as previously done, by clr normalization of expression matrix then 
#rho to calculate coexpression followed by spatial quantile normalization

library(dplyr)
library(RMariaDB)
library(spqn)

#get clr expression matrix
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

rho<-function(gene1,gene2)
{
  idx2keep<-intersect(which(!is.na(gene1)),which(!is.na(gene2)))
  gene1<-gene1[idx2keep]
  gene2<-gene2[idx2keep]
  return(1- (var(gene1-gene2,na.rm = T)/(var(gene1,na.rm = T)+var(gene2,na.rm = T) )))
}

getCLRdata<-function(raw_cut, sample_cut, raw_data){
  plot_df<-data.frame('transcript'=rownames(raw_data),
                      'sample_count'=apply(raw_data,1,function(x){length(which(x>raw_cut))}),stringsAsFactors = F)
  
  transcripts2keep<-filter(plot_df,sample_count>=sample_cut)$transcript
  raw_data<-raw_data[transcripts2keep,]
  #add a psuedocount of 1 to all values
  raw_data<-raw_data+1
  #do clr transformation
  #clr_data<-apply(raw,2,clr)
  clr_data<-clr(raw_data)
  return(clr_data)
}


raw<-readRDS('raw_counts.RDS')
clr_nozeros_same_genes_as_paper<-getCLRdata(raw_cut=5, sample_cut=400,  raw_data=raw)

saveRDS(clr_nozeros_same_genes_as_paper, 'clr_raw5_sample400_noNA.RDS')

dim(clr_nozeros_same_genes_as_paper) #11630

#calculate coexpression, do this for all columns seperately and then rejoin
rho_output<-data.frame('transcript'=rownames(clr_nozeros_same_genes_as_paper), 'rho'=NA)
colnames(rho_output)[2]<-rownames(clr_nozeros_same_genes_as_paper)[orf_idx]

rho_output<-data.frame('transcript'=rownames(clr_nozeros_same_genes_as_paper), 'rho'=NA)
colnames(rho_output)[2]<-rownames(clr_nozeros_same_genes_as_paper)[orf_idx]

rho<-function(gene1,gene2)
{
  return(1- (var(gene1-gene2)/(var(gene1)+var(gene2) )))
}

for(i in 1:nrow(clr_nozeros_same_genes_as_paper)){
  rho_output[i,2]<-rho(clr_nozeros_same_genes_as_paper[orf_idx,], clr_nozeros_same_genes_as_paper[i,])
}

write.table(rho_output,file = sprintf('%s/rho_%d_%s.txt', output_file_path, orf_idx,rownames(clr_nozeros_same_genes_as_paper)[orf_idx]),row.names = FALSE,quote=FALSE)

files<-list.files(path = input_files_path, full.names = TRUE, recursive = FALSE)
x<-read.delim(files[1],stringsAsFactors = FALSE,sep = ' ')
for(i in 2:length(files)){
  temp<-read.delim(files[i],stringsAsFactors = FALSE,sep = ' ')
  x<-left_join(x,temp)
}

rownames(x)<-x$transcript
x$transcript<-NULL
saveRDS(x,'rho_raw5_sample400_noNA.RDS')


#normalize coexpression matrix using SpQN
#read in coexpression matrix
rho_matrix<-readRDS('rho_raw5_sample400_noNA.RDS')
#read in normalized expression data
clr_data<-readRDS('clr_raw5_sample400_noNA.RDS')

#get median clr expression level across all samples for each ORF
exp_df<-data.frame(transcript=rownames(clr_data),mean_value=rowMeans(clr_data,na.rm = T),stringsAsFactors = F)
rownames(exp_df)<-exp_df$transcript
all(rownames(rho_matrix)==rownames(clr_data))

cor_m_spqn <- normalize_correlation(as.matrix(rho_matrix), ave_exp=exp_df$mean_value, ngrp=20, size_grp=300, ref_grp=18)
rownames(cor_m_spqn)<-rownames(rho_matrix)
colnames(cor_m_spqn)<-colnames(rho_matrix)
saveRDS(cor_m_spqn,'spqn_rho_raw5_sample400_noNA.RDS')
