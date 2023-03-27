#script to plot distribution of coexpression values binned by expression level before and after applying spqn (from wang et al 2022)
library(dplyr)
library(ggplot2)
library(scales)
library(spqn)

quick_hist = function(values_df, breaks=30,num_bins=10,hex_color) {
  res = hist(filter(values_df,bin==1)$cor, plot=FALSE, breaks=breaks)
  
  dat = data.frame(xmin=head(res$breaks, -1L),
                   xmax=tail(res$breaks, -1L),
                   ymin=0.0,
                   ymax=res$counts,
                   bins=1)
  for(i in 2:num_bins)
  {
    res = hist(filter(values_df,bin==i)$cor, plot=FALSE, breaks=breaks)
    
    temp = data.frame(xmin=head(res$breaks, -1L),
                      xmax=tail(res$breaks, -1L),
                      ymin=0.0,
                      ymax=res$counts,
                      bins=i)
    dat<-rbind(dat,temp)
  }
  
  p<- ggplot(dat, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) +
    geom_rect(size=0.5, colour="white", fill=hex_color)+ facet_wrap(~bins,nrow = num_bins)+theme_minimal()+
    theme(text=element_text(size=10))+
    scale_x_continuous(limits = c(-0.5,1))+
    labs(x='rho',y='count')
  print(p)
}

flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}

get_hist_exp_plot<-function(expinfo,coexpinfo,numberofbins,output_file){
  plot_df<-data.frame('transcript'=rownames(expinfo),'mean_value'=rowMeans(expinfo,na.rm = T),stringsAsFactors = FALSE)
  plot_df<-plot_df[order(plot_df$mean_value),]
  binsize<- floor(nrow(plot_df)/numberofbins)
  x<-flattenCorrMatrix(coexpinfo[plot_df$transcript[1:binsize],plot_df$transcript[1:binsize]])
  x$bin<-1
  for(i in 2:numberofbins)
  {
    idx<-binsize*i
    temp<- flattenCorrMatrix(coexpinfo[plot_df$transcript[(idx-binsize+1):idx],plot_df$transcript[(idx-binsize+1):idx]])
    temp$bin<-i
    x<-rbind(x,temp)
  }
  
  pdf(output_file,width=2,height=5)
  quick_hist(x, breaks=30,num_bins = numberofbins, 'blue') 
  dev.off()
}

#import coexpression values
rho_matrix_5_400<-readRDS('rho_raw5_sample400.RDS')

#import clr normalized expression values
clr_data_5_400<-readRDS('clr_raw5_sample400.RDS')

#get mean clr value across all samples for each ORF 
exp_df_5_400<-data.frame(transcript=rownames(clr_data_5_400),mean_value=rowMeans(clr_data_5_400,na.rm = T),stringsAsFactors = F)
rownames(exp_df_5_400)<-exp_df_5_400$transcript
all(rownames(rho_matrix_5_400)==rownames(clr_data_5_400))

#normalize coexpression values based on expression level (mean clr) 
cor_m_spqn_5_400 <- normalize_correlation(as.matrix(rho_matrix_5_400), ave_exp=exp_df_5_400$mean_value, ngrp=20, size_grp=300, ref_grp=18)
rownames(cor_m_spqn_5_400)<-rownames(rho_matrix_5_400)
colnames(cor_m_spqn_5_400)<-colnames(rho_matrix_5_400)

#make plots
get_hist_exp_plot(clr_data_5_400,rho_matrix_5_400,5,'./20230106_figures/rho_dist_unnorm.pdf')
get_hist_exp_plot(clr_data_5_400,cor_m_spqn_5_400,5,'./20230106_figures/rho_dist_spqn_norm.pdf')

