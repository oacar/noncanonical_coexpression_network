#script to generate the TF enrichments and plot for supplementary figure 10 looking at the odds ratio for association
#between being expressed and having promoter bound by a common transcription factor at various network cutoffs that lead to 
#different proportions of nORFs being present in the network
#comparing constructing networks using 'clr + rho + spqn' vs networks constructed using 'clr + rho' 

#first get the prop of nORFs present in the network at various cutoffs

library(RMariaDB)
library(dplyr)

#load orf annotation info
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
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
nORF_df<-filter(orf_info, transcript %in% rownames(rho) & is_canonical=='noncanonical')$transcript


get_prop_nORFs<-function(coexp_df, nORFs, orf_df){
  rf<-flattenCorrMatrix(coexp_df)
  rf<-filter(rf, !is.na(cor))
  rf<-left_join(rf,orf_df[,c('transcript','is_canonical')],by=c('row'='transcript'))
  rf<-left_join(rf,orf_df[,c('transcript','is_canonical')],by=c('column'='transcript'))
  rf<-filter(rf, !(is_canonical.x=='canonical' & is_canonical.y=='canonical'))
  cutoffs<-seq(from=0.4, to=1, by=0.01)
  
  df<-data.frame('cor_cutoffs'=cutoffs,'prop_nORFs'=NA)
  for(i in 1:length(cutoffs)){
    rf<-filter(rf, cor>cutoffs[i])
    df$prop_nORFs[i]<-length(unique(c(filter(rf, is_canonical.x=='noncanonical')$row, filter(rf, is_canonical.y=='noncanonical')$column)))/length(nORFs)
  }
  df$prop_nORFs_rounded<-round(df$prop_nORFs,digits=2)
  return(df)
}



#clr + rho + SpQN (method used in manuscript)
rho<-readRDS('/spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA
clr_rho_SpQN<-get_prop_nORFs(rho, nORF_df, orf_info)
clr_rho_SpQN$method<-'clr normalization + rho (similarity metric) + SpQN our method'
saveRDS(clr_rho_SpQN, 'clr_rho_SpQN_prop_nORF.RDS')

#clr + rho 
clr_rho<-readRDS('rho_raw5_sample400.RDS')
clr_rho[num_obs<400]<-NA
clr_rho<-get_prop_nORFs(clr_rho, nORF_df, orf_info)
clr_rho$method<-'clr normalization + rho (similarity metric)'
saveRDS(clr_rho, 'clr_rho_prop_nORF.RDS')


#get TFBS enrichment at the various cutoffs
TF_distance<-readRDS('ORF_TF_distances_TSS_all.RDS')
TF_distance[TF_distance<=200]<-1
TF_distance[TF_distance>200]<-0
TF_distance[is.na(TF_distance)]<-0
TF_count<-rowSums(TF_distance,na.rm=TRUE)
orfsWithTFBS<-names(which(TF_count>0))
orfsWithTFBS<-intersect(orfsWithTFBS,rownames(rho))

TFsincommon<-readRDS('TFsincommonMatrix_TSS_200_all.rds')
TFsincommon<-TFsincommon[orfsWithTFBS,orfsWithTFBS]

get_enrichment_stats<-function(flatten_coexp_data, cutoff, comp_type){
  res<-mutate(flatten_coexp_data,coexpressed=case_when(cor> cutoff ~ 1,
                                                       TRUE ~ 0))
  cm<-table(res$TFsincommon,res$coexpressed)
  if(dim(cm)[2]!=2){
    cm<-matrix(nrow = 2, ncol=2)
    cm[1,1]<- filter(res, coexpressed==0 & TFsincommon==0) %>% nrow()
    cm[2,1]<-filter(res, coexpressed==0 & TFsincommon==1) %>% nrow()
    cm[2,2]<-filter(res, coexpressed==1 & TFsincommon==1) %>% nrow()
    cm[1,2]<-filter(res, coexpressed==1 & TFsincommon==0) %>% nrow()
  }
  x<-fisher.test(cm)
  return(data.frame('cutoff'=cutoff,'odds_ratio'=x$estimate, 'p_value'=x$p.value, 'upper_ci'=x$conf.int[1], 'lower_ci'=x$conf.int[2], 'comparison_type'=comp_type))
  
}

getTFBSenrichment<-function(TF_data,coexp_data,orf_df){
  coexp_data<-coexp_data[rownames(TF_data),colnames(TF_data)]
  coexp_data<-flattenCorrMatrix(coexp_data)
  temp<-flattenCorrMatrix(TF_data)
  colnames(temp)[3]<-'TFsincommon'
  coexp_data<-left_join(temp,coexp_data)
  coexp_data<-filter(coexp_data,!is.na(cor))
  coexp_data<-left_join(coexp_data,orf_df[,c('transcript','is_canonical')],by=c('row'='transcript'))
  coexp_data<-left_join(coexp_data,orf_df[,c('transcript','is_canonical')],by=c('column'='transcript'))
  cc_coexp_data<-filter(coexp_data,is_canonical.x=='canonical' & is_canonical.y=='canonical')
  cn_coexp_data<-filter(coexp_data, (is_canonical.x=='canonical' & is_canonical.y=='noncanonical') |
                          (is_canonical.x=='noncanonical' & is_canonical.y=='canonical'))
  nn_coexp_data<-filter(coexp_data,is_canonical.x=='noncanonical' & is_canonical.y=='noncanonical')
  
  cutoffs<-seq(from=0.4, to=1, by=0.01)
  df<-NA
  for(i in 1:length(cutoffs)){
    df<-rbind(df,get_enrichment_stats(cc_coexp_data, cutoff = cutoffs[i], comp_type = 'cc'))
    df<-rbind(df,get_enrichment_stats(cn_coexp_data, cutoff = cutoffs[i], comp_type = 'cn'))
    df<-rbind(df,get_enrichment_stats(nn_coexp_data, cutoff = cutoffs[i], comp_type = 'nn'))
  }
  df<-filter(df, !is.na(cutoff))
  return(df)
}

#clr + rho + SpQN 
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

rho_spqn_df<-getTFBSenrichment(TFsincommon,rho,orf_info)
rho_spqn_df$method<-'clr normalization + rho (similarity metric) + SpQN our method'
saveRDS(rho_spqn_df, 'clr_rho_SpQN_TFBS.RDS')

#clr + rho 
clr_rho<-readRDS('rho_raw5_sample400.RDS')
clr_rho[num_obs<400]<-NA
clr_rho<-getTFBSenrichment(TFsincommon,clr_rho,orf_info)
clr_rho$method<-'clr normalization + rho (similarity metric)'
saveRDS(clr_rho, 'clr_rho_TFBS.RDS')

#plot TFBS odds ratio vs prop nORF in network
library(ggplot2)
library(dplyr)
library(ggpubr)
library(forcats)


#read in nORF prop data
temp<-rbind(readRDS('clr_rho_SpQN_prop_nORF.RDS'),readRDS('clr_rho_prop_nORF.RDS'))

#read in TFBS odds ratio data
df<-rbind(readRDS('clr_rho_SpQN_TFBS.RDS'),readRDS('clr_rho_TFBS.RDS'))

df<-left_join(df,temp, by=c('cutoff'='cor_cutoffs','method'))

significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

df<-mutate(df, comparison_type=case_when(comparison_type=='cc' ~'cORF-cORF pairs',
                                         comparison_type=='cn'~ 'cORF-nORF pairs',
                                         comparison_type=='nn' ~ 'nORF-nORF pairs'))

pdf('nORFprop_vs_TFBSenrichment_spqnvsnospqn.pdf',width=6.5,height=3.2)
df %>% filter(method %in% c('clr normalization + rho (similarity metric) + SpQN our method', 'clr normalization + rho (similarity metric)')) %>% 
  mutate(method = factor(method, levels= c('clr normalization + rho (similarity metric) + SpQN our method','clr normalization + rho (similarity metric)'))) %>%
  filter(p_value<0.05 & !is.na(cutoff) & !is.infinite(odds_ratio)) %>% mutate(odds_ratio=log10(odds_ratio)) %>%
  mutate(upper_ci=log10(upper_ci)) %>% mutate(lower_ci=log10(lower_ci)) %>%
  ggline(x="prop_nORFs_rounded", y="odds_ratio",color="method",facet.by = "comparison_type", numeric.x.axis = TRUE,
         xlab='proportion of nORFs in the network', ylab='log odds ratio for association\nof coexpression and having\npromoter bound by a common TF',
         add = "loess", conf.int = TRUE)+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("legend.text", size=legend_text_size)+
  font("legend.title",size=legend_title_size)
dev.off()