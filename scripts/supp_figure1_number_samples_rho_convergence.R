#script to find the number of samples needed when caluclating coexpression for the coexpression values to converge/stabalize within +/-0.5 or +/-0.1 upon addition of 10 more samples
library(dplyr)
library(ggplot2)
library(RMariaDB)
library(RColorBrewer)

#read in rho
rho<-readRDS('rho_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA


#read in rna seq data 
raw<-readRDS('raw_counts.RDS')

#read in orf data
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
coexpression_orf_list<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


flattenCorrMatrix <- function(cormat){
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]
  )
}

addAnnoInfo<-function(flatcorrmat,coexpression_orf_list){
  flatcorrmat<-left_join(flatcorrmat,coexpression_orf_list[,c('transcript','is_canonical')],by=c("row"="transcript"))
  flatcorrmat<-left_join(flatcorrmat,coexpression_orf_list[,c('transcript','is_canonical')],by=c("column"="transcript"))
  flatcorrmat<-mutate(flatcorrmat,
                      pair_classification=case_when( 
                        is_canonical.x=='canonical' &  is_canonical.y=='canonical' ~ 'canonical-canonical',
                        (is_canonical.x=='noncanonical' & is_canonical.y=='canonical') | (is_canonical.x=='canonical' & is_canonical.y=='noncanonical') ~ 'canonical-noncanonical',
                        is_canonical.x=='noncanonical' &  is_canonical.y=='noncanonical' ~ 'noncanonical-noncanonical',
                        TRUE ~'other')) 
  return(flatcorrmat)
}

#clr function from Quinn et al (2017)
clr<-function(counts)
{ 
  counts<-t(counts)
  logX <- log(counts)
  ref <- rowMeans(logX,na.rm = T)
  lr <- sweep(logX, 1, ref, "-")
  lr<-t(lr)
  return(lr)
}

getClrData<-function(raw_data,raw_cut,sample_cut,NA_cut){
  plot_df<-data.frame('transcript'=rownames(raw_data),
                      'sample_count'=apply(raw_data,1,function(x){length(which(x>raw_cut))}),stringsAsFactors = F)
  
  transcripts2keep<-filter(plot_df,sample_count>=sample_cut)$transcript
  raw_data<-raw_data[transcripts2keep,]
  #convert 0s to NA
  raw_data[raw_data<NA_cut]<-NA

  clr_data<-clr(raw_data)
  return(clr_data)
}

#get cor_analysis df:
getOrfPairsForAnalysis<-function(rho_matrix,clr_data,coexp_pairs='yes'){
  rho_flat<-flattenCorrMatrix(rho_matrix)
  whole_cutoff<-quantile(rho_flat$cor,0.99, na.rm=TRUE)
  rho_flat<-addAnnoInfo(rho_flat,coexpression_orf_list)
  rho_flat<-filter(rho_flat, pair_classification=='canonical-noncanonical' & row %in% rownames(clr_data) & column %in% rownames(clr_data))
  if(coexp_pairs=='yes'){
    pairs_df<-rho_flat  %>% filter(cor>whole_cutoff) 
  }else{
    #ie use random pairs
    pairs_df<-rho_flat %>% sample_n(10000) 
  }
  cor_analysis<-data.frame('orf1'=pairs_df$row, 'orf2'=pairs_df$column,'rho'=NA,'orf1_nsamples'=NA,'orf2_nsamples'=NA,
                           'nsamples_total'=NA, 'nsample_where_cor_stabalizes_005'=NA,'nsample_where_cor_stabalizes_01'=NA,
                           'rho_005'=NA,'rho_01'=NA)
  return(cor_analysis)
}


clr_5<-getClrData(raw,raw_cut=5,sample_cut=400,NA_cut=5)
orf_pairs_5_coexp<-getOrfPairsForAnalysis(rho,clr_5,coexp_pairs='yes')

saveRDS(orf_pairs_5_coexp,'orf_pairs_raw5_coexp_input.RDS')
saveRDS(clr_5,'clr_data_raw5_sample400.RDS')

clr_data<-clr_5
cor_analysis<-orf_pairs_5_coexp

#function to calculate coexpression
rho<-function(gene1,gene2)
{
  idx2keep<-intersect(which(!is.na(gene1)),which(!is.na(gene2)))
  gene1<-gene1[idx2keep]
  gene2<-gene2[idx2keep]
  return(1- (var(gene1-gene2,na.rm = T)/(var(gene1,na.rm = T)+var(gene2,na.rm = T) )))
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
roundDown <- function(x,m) m*floor(x / m)

for(i in 1:nrow(cor_analysis)){
  exp_df<-data.frame('orf1'=clr_data[cor_analysis$orf1[i],], 'orf2'=clr_data[cor_analysis$orf2[i],])
  cor_analysis$orf1_nsamples[i]<-nrow(filter(exp_df, !is.na(orf1)))
  cor_analysis$orf2_nsamples[i]<-nrow(filter(exp_df, !is.na(orf2)))
  
  exp_df<-filter(exp_df,!is.na(orf1) & !is.na(orf2))
  if(nrow(exp_df)>20){
    cor_analysis$rho[i]<-rho(exp_df$orf1,exp_df$orf2)
    cor_analysis$nsamples_total[i]<-nrow(exp_df)
    
    nrows_power_df<-roundDown(nrow(exp_df),10)
    power_analysis<-data.frame('sample_size'=1:nrows_power_df,'rho'=NA)
    for(j in 1:nrows_power_df){
      samples2use<-exp_df %>% sample_n(j) 
      power_analysis$rho[j]<-rho(samples2use$orf1, samples2use$orf2)
    }
    idxs<-chunk2(power_analysis$rho,nrow(power_analysis)/10)
    cor_diffs<-vector(length = length(idxs))
    for(k in 1:length(idxs)){
      cor_diffs[k]<-max(idxs[[k]],na.rm=TRUE)-min(idxs[[k]],na.rm=TRUE)
    }
    if(!(all(cor_diffs>=0.05))){
      cor_analysis$nsample_where_cor_stabalizes_005[i]<-min(which(cor_diffs<=0.05))*10
      cor_analysis$rho_005[i]<-power_analysis$rho[cor_analysis$nsample_where_cor_stabalizes_005[i]]
    }
    if(!(all(cor_diffs>=0.1))){
      cor_analysis$nsample_where_cor_stabalizes_01[i]<-min(which(cor_diffs<=0.1))*10
      cor_analysis$rho_01[i]<-power_analysis$rho[cor_analysis$nsample_where_cor_stabalizes_01[i]]
    }
  }
}

saveRDS(cor_analysis,'orf_pairs_raw5_coexp.RDS')


#plot 
#make histogram for distrbution of sample size for stabilization at each sample cutoff
library(dplyr)
library(RMariaDB)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(scales)

axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)



getNumSamplesHistogram<-function(pair_df,plot_filename){
  pair_df_melted<- melt(pair_df, id.vars=c('orf1','orf2','orf1_nsamples','orf2_nsamples','nsamples_total',
                                           'rho_005','rho_01','rho'))
  
  pair_df_melted<-mutate(pair_df_melted,fluctuations=case_when(variable == 'nsample_where_cor_stabalizes_005' ~ 'fluctuation <= 0.05',
                                                               variable == 'nsample_where_cor_stabalizes_01' ~ 'fluctuation <= 0.1' ))
  
  median_df<-data.frame('fluctuations'=c('fluctuation <= 0.05','fluctuation <= 0.1'),
                        'median'=c(median(filter(pair_df_melted, fluctuations=='fluctuation <= 0.05')$value,na.rm=TRUE),
                                   median(filter(pair_df_melted, fluctuations=='fluctuation <= 0.1')$value,na.rm=TRUE)))
  pair_df_melted<-mutate(pair_df_melted, fluctuations = 
                           case_when(fluctuations=='fluctuation <= 0.05' ~ sprintf('fluctuation <= 0.05\nmedian samples=%d', median_df$median[which(median_df$fluctuations=='fluctuation <= 0.05')]),
                                     fluctuations=='fluctuation <= 0.1' ~ sprintf('fluctuation <= 0.1\nmedian samples=%d', median_df$median[which(median_df$fluctuations=='fluctuation <= 0.1')])))
  median_df<-mutate(median_df, fluctuations = 
                      case_when(fluctuations=='fluctuation <= 0.05' ~ sprintf('fluctuation <= 0.05\nmedian samples=%d', median_df$median[which(median_df$fluctuations=='fluctuation <= 0.05')]),
                                fluctuations=='fluctuation <= 0.1' ~ sprintf('fluctuation <= 0.1\nmedian samples=%d', median_df$median[which(median_df$fluctuations=='fluctuation <= 0.1')])))
  
  
  p<-pair_df_melted %>% 
    ggplot(aes(x=value))+geom_histogram()+facet_grid(~fluctuations)+theme_minimal()+theme(axis.text=element_text(size=axis_text_size),
                                                                                          axis.title = element_text(size=axis_title_size)) +
    labs(x='number of samples need for rho to converge',y='number of ORF pairs')+
    geom_vline(data = median_df,aes(xintercept = median),lty='dashed',color='red', lwd=1.3)
  
  pdf(sprintf('%s',plot_filename),width=6,height=3)
  print(p)
  dev.off()
}

getMedianDf<-function(pair_df){
  pair_df_melted<- melt(pair_df, id.vars=c('orf1','orf2','orf1_nsamples','orf2_nsamples','nsamples_total',
                                           'rho_005','rho_01','rho'))
  
  pair_df_melted<-mutate(pair_df_melted,fluctuations=case_when(variable == 'nsample_where_cor_stabalizes_005' ~ 0.05,
                                                               variable == 'nsample_where_cor_stabalizes_01' ~ 0.1 ))
  
  median_df<-data.frame('fluctuations'=c(0.05,0.1),
                        'median'=c(median(filter(pair_df_melted, fluctuations==0.05)$value,na.rm=TRUE),
                                   median(filter(pair_df_melted, fluctuations==0.1)$value,na.rm=TRUE)))
  
  return(median_df)
  
}


getNumSamplesHistogram(cor_analysis,'raw5_coexpPairs_NumSamples2Converge_histogram.pdf')




