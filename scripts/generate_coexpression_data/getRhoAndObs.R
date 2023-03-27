#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
orf_idx<-as.integer(args[1])
clr_data_path<-args[2]
output_file_path<-args[3]

start_time <- Sys.time()
clr_data<-readRDS(clr_data_path)
rho_output<-data.frame('transcript'=rownames(clr_data), 'rho'=NA)
colnames(rho_output)[2]<-rownames(clr_data)[orf_idx]

num_obs_output<-data.frame('transcript'=rownames(clr_data),'num_obs'=NA)
colnames(num_obs_output)[2]<-rownames(clr_data)[orf_idx]

orf1_keep_idx<-which(!is.na(clr_data[orf_idx,]))
clr_data<-clr_data[,orf1_keep_idx]

getSamplesInCommon<-function(gene2){
  idx2keep<-which(!is.na(gene2))
  return(idx2keep)
}

rho<-function(gene1,gene2,idx2keep)
{
  gene1<-gene1[idx2keep]
  gene2<-gene2[idx2keep]
  return(1- (var(gene1-gene2)/(var(gene1)+var(gene2) )))
}


for(i in 1:nrow(clr_data)){
  keep_idx<-getSamplesInCommon(clr_data[i,])
  num_obs_output[i,2]<-length(keep_idx)
  rho_output[i,2]<-rho(clr_data[orf_idx,], clr_data[i,],keep_idx)
}

write.table(rho_output,file = sprintf('%s/rho/rho_%d_%s.txt', output_file_path, orf_idx,rownames(clr_data)[orf_idx]),row.names = FALSE,quote=FALSE)
write.table(num_obs_output,file = sprintf('%s/num_obs/numobs_%d_%s.txt', output_file_path, orf_idx,rownames(clr_data)[orf_idx]),row.names = FALSE,quote=FALSE)
end_time <- Sys.time()
print(end_time - start_time)


