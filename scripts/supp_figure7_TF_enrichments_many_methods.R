#script to plot supplementary figure 7: odds ratio of TF enrichments (ie is there an association between an 
#ORF being coexpressed and having their promoters bound by a common TF) using different methods to construct coexpression network

#data to generate these plots is from:
# contingency matrices for method 'TPM + pearson correlation' is from: get_TF_enrichment_pearson_tpm.R
# contingency matrices for method 'RPKM + pearson correlation' is from: get_TF_enrichment_pearson_rpkm.R
# contingency matrices for method 'clr + rho + spqn' is from: get_TF_enrichment_clr_rho_spqn.R
# contingency matrices for method 'clr + batch correction + rho + spqn' is from: get_TF_enrichment_clr_PC1_rho_spqn.R

library(ggplot2)
library(dplyr)
library(ggpubr)
library(forcats)

method_levels <- c('clr normalization + rho (similarity metric) + SpQN', 'TPM normalization + pearson correlation (similarity metric)', 
                   'clr normalization + batch correction + rho (similarity metric) + SpQN', 'RPKM normalization + pearson correlation (similarity metric)')

get_stats_df<-function(cutoff,coexp_type){
  if(coexp_type=='pearson_tpm'){
    ccm<-readRDS(sprintf('pearson_tpm_cc_pairswithTFBS_200_%i.RDS',cutoff))
    cnm<-readRDS(sprintf('pearson_tpm_cn_pairswithTFBS_200_%i.RDS',cutoff))
    nnm<-readRDS(sprintf('pearson_tpm_nn_pairswithTFBS_200_%i.RDS',cutoff))
    plot_df<-data.frame('cutoff'=cutoff,'odds_ratio'=NA,'pval'=NA,'method'='TPM normalization + pearson correlation (similarity metric)',
                        'orf_type'=c('cc','cn','nn'),'upper_ci'=NA,'lower_ci'=NA)
  }else if(coexp_type=='rpkm'){
    ccm<-readRDS(sprintf('pearson_rpkm_cc_pairswithTFBS_200_%i.RDS',cutoff))
    cnm<-readRDS(sprintf('pearson_rpkm_cn_pairswithTFBS_200_%i.RDS',cutoff))
    nnm<-readRDS(sprintf('pearson_rpkm_nn_pairswithTFBS_200_%i.RDS',cutoff))
    plot_df<-data.frame('cutoff'=cutoff,'odds_ratio'=NA,'pval'=NA,'method'='RPKM normalization + pearson correlation (similarity metric)',
                        'orf_type'=c('cc','cn','nn'),'upper_ci'=NA,'lower_ci'=NA)
  }else if(coexp_type=='pc1'){
    ccm<-readRDS(sprintf('clr_rho_pc1_spqn_cc_pairswithTFBS_200_%i.RDS',cutoff))
    cnm<-readRDS(sprintf('clr_rho_pc1_spqn_cn_pairswithTFBS_200_%i.RDS',cutoff))
    nnm<-readRDS(sprintf('clr_rho_pc1_spqn_nn_pairswithTFBS_200_%i.RDS',cutoff))
    plot_df<-data.frame('cutoff'=cutoff,'odds_ratio'=NA,'pval'=NA,'method'='clr normalization + batch correction + rho (similarity metric) + SpQN',
                        'orf_type'=c('cc','cn','nn'),'upper_ci'=NA,'lower_ci'=NA)
  }else{
    ccm<-readRDS(sprintf('rho_cc_pairswithTFBS_200_%i.RDS',cutoff))
    cnm<-readRDS(sprintf('rho_cn_pairswithTFBS_200_%i.RDS',cutoff))
    nnm<-readRDS(sprintf('rho_nn_pairswithTFBS_200_%i.RDS',cutoff))
    plot_df<-data.frame('cutoff'=cutoff,'odds_ratio'=NA,'pval'=NA,'method'='clr normalization + rho (similarity metric) + SpQN', 
                        'orf_type'=c('cc','cn','nn'),'upper_ci'=NA,'lower_ci'=NA)
    
  }
  ccm_test<-fisher.test(ccm)
  plot_df$odds_ratio[which(plot_df$orf_type=='cc')]<-ccm_test$estimate
  plot_df$pval[which(plot_df$orf_type=='cc')]<-ccm_test$p.value
  plot_df$lower_ci[which(plot_df$orf_type=='cc')]<-ccm_test$conf.int[1]
  plot_df$upper_ci[which(plot_df$orf_type=='cc')]<-ccm_test$conf.int[2]
  
  cnm_test<-fisher.test(cnm)
  plot_df$odds_ratio[which(plot_df$orf_type=='cn')]<-cnm_test$estimate
  plot_df$pval[which(plot_df$orf_type=='cn')]<-cnm_test$p.value
  plot_df$lower_ci[which(plot_df$orf_type=='cn')]<-cnm_test$conf.int[1]
  plot_df$upper_ci[which(plot_df$orf_type=='cn')]<-cnm_test$conf.int[2]
  
  nnm_test<-fisher.test(nnm)
  plot_df$odds_ratio[which(plot_df$orf_type=='nn')]<-nnm_test$estimate
  plot_df$pval[which(plot_df$orf_type=='nn')]<-nnm_test$p.value
  plot_df$lower_ci[which(plot_df$orf_type=='nn')]<-nnm_test$conf.int[1]
  plot_df$upper_ci[which(plot_df$orf_type=='nn')]<-nnm_test$conf.int[2]
  
  return(plot_df)
}
df<-NA
for(cutoff in c(90,95,99,998,999)){
  df<-rbind(df,get_stats_df(cutoff,'pearson_tpm'))
}
for(cutoff in c(90,95,99,998,999)){
  df<-rbind(df,get_stats_df(cutoff,'clr_rho'))
}
for(cutoff in c(90,95,99,998,999)){
  df<-rbind(df,get_stats_df(cutoff,'pc1'))
}
for(cutoff in c(90,95,99,998,999)){
  df<-rbind(df,get_stats_df(cutoff,'rpkm'))
}

df<-filter(df, !is.na(cutoff))
df$cutoff<-as.factor(df$cutoff)

significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

pdf('./20230106_figures/TFBS_manymethodsVsRho.pdf',width=4.5,height=3.2)
df %>% mutate(method = factor(method, levels= method_levels)) %>% filter(pval<0.05 & !is.na(cutoff)) %>% 
  ggline(x="cutoff", y="odds_ratio",color="method",facet.by = "orf_type",
         xlab='percentile cutoff for defining coexpression', ylab='odds ratio for association\nof coexpression and sharing TF')+
  geom_errorbar(aes(ymax = upper_ci, ymin = lower_ci), width = 0.1)+
  geom_hline(yintercept=1,linetype='dashed')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)
dev.off()
