library(dplyr)
library(spqn)

#read in coexpression matrix
rho_matrix_5_400<-readRDS('rho_raw5_sample400.RDS')
#read in normalized expression data
clr_data_5_400<-readRDS('clr_raw5_sample400.RDS')

#get median clr expression level across all samples for each ORF
exp_df_5_400<-data.frame(transcript=rownames(clr_data_5_400),mean_value=rowMeans(clr_data_5_400,na.rm = T),stringsAsFactors = F)
rownames(exp_df_5_400)<-exp_df_5_400$transcript
all(rownames(rho_matrix_5_400)==rownames(clr_data_5_400))

#run spqn on coexpression matrix
cor_m_spqn_5_400 <- normalize_correlation(as.matrix(rho_matrix_5_400), ave_exp=exp_df_5_400$mean_value, ngrp=20, size_grp=300, ref_grp=18)
rownames(cor_m_spqn_5_400)<-rownames(rho_matrix_5_400)
colnames(cor_m_spqn_5_400)<-colnames(rho_matrix_5_400)
saveRDS(cor_m_spqn_5_400,'spqn_raw5_sample400.RDS')

