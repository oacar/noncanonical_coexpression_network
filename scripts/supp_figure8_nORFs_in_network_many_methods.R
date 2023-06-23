# script to plot supp figure 8: proportion of nORFs in the network at each threshold for the different methods to calculate coexpression
# data to generate these plots is from:
# coexpression matrix for 'TPM + pearson correlation' is from: get_TPM_pearson_coexpression.R
# coexpression matrix for 'RPKM + pearson correlation' is from: get_RPKM_pearson_coexpression.R
# coexpression matrix for 'clr + rho + spqn' is from: /generate_coexpression_data/runSpqn.R
# coexpression matrix for 'clr + batch correction + rho + spqn' is from: get_clr_pc1_expression_data.R


library(dplyr)
library(ggpubr)
library(RMariaDB)
library(forcats)

#load orf annotation info
#get ORF data

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}

rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA


rho_pc1<-readRDS('spqn_raw5_sample400_clr_pc1.RDS')
rho_pc1[num_obs<400]<-NA

pearson<-readRDS('tpm_raw5_sample400_coexp_pearson.RDS')
pearson[num_obs<400]<-NA

rpkm<-readRDS('rpkm_raw5_sample400_coexp_pearson.RDS')
rpkm[num_obs<400]<-NA

nORFs<-filter(orf_info, transcript %in% rownames(rho) & is_canonical=='noncanonical')$transcript

get_prop_nORFs<-function(coexp_df, nORFs, orf_df){
  
  rf<-flattenCorrMatrix(coexp_df)
  cutoff_rf999<-quantile(rf$cor,0.999,na.rm=TRUE)
  cutoff_rf998<-quantile(rf$cor,0.998,na.rm=TRUE)
  cutoff_rf99<-quantile(rf$cor, 0.99, na.rm=TRUE)
  cutoff_rf95<-quantile(rf$cor,0.95,na.rm=TRUE)
  cutoff_rf90<-quantile(rf$cor, 0.90, na.rm=TRUE)
  
  rf<-left_join(rf,orf_df[,c('transcript','is_canonical')],by=c('row'='transcript'))
  rf<-left_join(rf,orf_df[,c('transcript','is_canonical')],by=c('column'='transcript'))
  df<-data.frame('cutoffs'=c('0.9','0.95','0.99','0.998','0.999'),'prop_nORFs'=NA)
  
  df$prop_nORFs[which(df$cutoff=='0.9')]<-length(unique(c(filter(rf, cor > cutoff_rf90 & is_canonical.x=='noncanonical')$row, filter(rf, cor > cutoff_rf90 & is_canonical.y=='noncanonical')$column)))/length(nORFs)
  df$prop_nORFs[which(df$cutoff=='0.95')]<-length(unique(c(filter(rf, cor > cutoff_rf95 & is_canonical.x=='noncanonical')$row, filter(rf, cor > cutoff_rf95 & is_canonical.y=='noncanonical')$column)))/length(nORFs)
  df$prop_nORFs[which(df$cutoff=='0.99')]<-length(unique(c(filter(rf, cor > cutoff_rf99 & is_canonical.x=='noncanonical')$row, filter(rf, cor > cutoff_rf99 & is_canonical.y=='noncanonical')$column)))/length(nORFs)
  df$prop_nORFs[which(df$cutoff=='0.998')]<-length(unique(c(filter(rf, cor > cutoff_rf998 & is_canonical.x=='noncanonical')$row, filter(rf, cor > cutoff_rf998 & is_canonical.y=='noncanonical')$column)))/length(nORFs)
  df$prop_nORFs[which(df$cutoff=='0.999')]<-length(unique(c(filter(rf, cor > cutoff_rf999 & is_canonical.x=='noncanonical')$row, filter(rf, cor > cutoff_rf999 & is_canonical.y=='noncanonical')$column)))/length(nORFs)
  return(df)
}

rho_df<-get_prop_nORFs(rho,nORFs, orf_info)
rho_df$network_type<-'clr normalization + rho (similarity metric) + SpQN our method'
rho_pc1_df<-get_prop_nORFs(rho_pc1,nORFs, orf_info)
rho_pc1_df$network_type<-'clr normalization + batch correction + rho (similarity metric) + SpQN'
pearson_df<-get_prop_nORFs(pearson,nORFs, orf_info)
pearson_df$network_type<-'TPM normalization + pearson correlation (similarity metric)'
rpkm_df<-get_prop_nORFs(rpkm,nORFs, orf_info)
rpkm_df$network_type<-'RPKM normalization + pearson correlation (similarity metric)'


plot_df<-rbind(rho_df,rho_pc1_df)
plot_df<-rbind(plot_df, pearson_df)
plot_df<-rbind(plot_df, rpkm_df)

significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

method_levels <- c('clr normalization + rho (similarity metric) + SpQN our method', 'TPM normalization + pearson correlation (similarity metric)', 
                   'clr normalization + batch correction + rho (similarity metric) + SpQN', 'RPKM normalization + pearson correlation (similarity metric)')

pdf('./20230106_figures/numberOFSinNetwork_plots.pdf',width=4.5,height=3.2)
plot_df %>%   mutate(network_type = factor(network_type, levels= method_levels)) %>%
  ggline(x="cutoffs", y="prop_nORFs",color="network_type",
         xlab='coexpression percentile cutoff for defining network', ylab='proportion of nORFs in the network')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)
dev.off()


