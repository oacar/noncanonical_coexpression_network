#script to plot the association between being coexpressed and having a common TF bound 200 bp upsteam of TSS
#defining coexpressed as rho>99.8th percentile and TF bound using chip exo data from rossi et al
#need to run scripts in get_TSS.R, get_TF_binding_distance_matrix.R and get_TFs_in_common.R first
library(dplyr)
library(RMariaDB)

#load orf info
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

#import coexpression matrix
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA


#get rho cutoffs before sub setting the network

rf<-flattenCorrMatrix(rho)
cutoff_rf998<-quantile(rf$cor,0.998,na.rm=TRUE)

#subset to contain only orfs that have at least 1 TFBS
TF_distance<-readRDS('ORF_TF_distances_TSS_all.RDS')

TF_distance[TF_distance<=200]<-1
TF_distance[TF_distance>200]<-0
TF_distance[is.na(TF_distance)]<-0

TF_count<-rowSums(TF_distance,na.rm=TRUE)
orfsWithTFBS<-names(which(TF_count>0))

orfsWithTFBS<-intersect(orfsWithTFBS,rownames(rho))
# TF_distance<-TF_distance[orfsInCommon,]

rho<-rho[orfsWithTFBS,orfsWithTFBS]

#import TFs in common matrix 
#read in TFs in common matrix (orf x orf matrix where 0 corresponds to having no TFs in common and value of 1
# means that orf pairs share atleast 1 TF in common.) #using 200 bp upstream of TSS for cutoff 
TFsincommon<-readRDS('TFsincommonMatrix_TSS_200_all.rds')



TFsincommon<-TFsincommon[orfsWithTFBS,orfsWithTFBS]
# rho<-rho[orfsInCommon,orfsInCommon]

getContingencyMatrix<-function(TF_data,coexp_data,cutoff,orf_df,comparison_type){
  coexp_data<-flattenCorrMatrix(coexp_data)
  temp<-flattenCorrMatrix(TF_data)
  colnames(temp)[3]<-'TFsincommon'
  coexp_data<-left_join(temp,coexp_data)
  coexp_data<-filter(coexp_data,!is.na(cor))
  coexp_data<-left_join(coexp_data,orf_df[,c('transcript','is_canonical')],by=c('row'='transcript'))
  coexp_data<-left_join(coexp_data,orf_df[,c('transcript','is_canonical')],by=c('column'='transcript'))
  if(comparison_type=='cc'){
    coexp_data<-filter(coexp_data,is_canonical.x=='canonical' & is_canonical.y=='canonical')
  }else if(comparison_type=='cn'){
    coexp_data<-filter(coexp_data, (is_canonical.x=='canonical' & is_canonical.y=='noncanonical') |
                         (is_canonical.x=='noncanonical' & is_canonical.y=='canonical'))
  }else if(comparison_type=='nn'){
    coexp_data<-filter(coexp_data,is_canonical.x=='noncanonical' & is_canonical.y=='noncanonical')
  }else {}
  
  result<-mutate(coexp_data,coexpressed=case_when(cor> cutoff ~ 1,
                                                  TRUE ~ 0))
  cm=table(result$TFsincommon,result$coexpressed)
  return(cm)
}



#association between being coexpressed and sharing a TF for canonical-canonical ORF pairs
cm_rho_998_cc<-getContingencyMatrix(TFsincommon,rho,cutoff_rf998,orf_info,comparison_type='cc')

#association between being coexpressed and sharing a TF for noncanonical-noncanonical ORF pairs
cm_rho_998_nn<-getContingencyMatrix(TFsincommon,rho,cutoff_rf998,orf_info,comparison_type='nn')

#association between being coexpressed and sharing a TF for canonical-noncanonical ORF pairs
cm_rho_998_cn<-getContingencyMatrix(TFsincommon,rho,cutoff_rf998,orf_info,comparison_type='cn')


#make plots
getPlotdf<-function(cm){
  cp = (cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = (cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed','not coexpressed'),
                      'prop'=c(cp, ncp),'count'=c(cm['1','1'],cm['1','0']))
  plot_df$se<-sqrt(plot_df$prop * (1 - plot_df$prop) / plot_df$count)
  plot_df$prop<-plot_df$prop*100
  plot_df$se<-plot_df$se*100
  return(plot_df)
}

library(ggplot2)
library(dplyr)
library(scales)
library(rstatix)
library(tibble)
library(ggpubr)

file_path<-'./20230106_figures/'
significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

#convert contingency matrix into a data frame that contains the proportion of coexpressed pairs that have shared TF and 
#proportion of noncoexpressed pairs that share a TF
cc<-getPlotdf(cm_rho_998_cc)
cc$pair_type<-'canonical\ncanonical'

cn<-getPlotdf(cm_rho_998_cn)
cn$pair_type<-'canonical\nnoncanonical'

nn<-getPlotdf(cm_rho_998_nn)
nn$pair_type<-'noncanonical\nnoncanonical'

plot_df<-rbind(cc,cn)
plot_df<-rbind(plot_df,nn)

stats_cc<-pairwise_fisher_test(ccm, p.adjust.method = "holm",detailed=TRUE) %>% add_column('pair_type'='canonical\ncanonical')
#cc odds ratio = 4.2848 p-value <2.2e-16
stats_cn<-pairwise_fisher_test(cnm, p.adjust.method = "holm") %>% add_column('pair_type'='canonical\nnoncanonical')
#cn odds ratio = 2.9981 p-value < 2.2e-16
stats_nn<-pairwise_fisher_test(nnm, p.adjust.method = "holm") %>% add_column('pair_type'='noncanonical\nnoncanonical')
#nn odds ratio = 3.8575 p-value <2.2e-16
stats_tf<- bind_rows(stats_cc,stats_cn,stats_nn)


pdf(sprintf('%sTFBS_Enrichment_rho_cc_cn_nn_998.pdf',file_path),width=3.5, height = 2.5)
ggbarplot(plot_df, x="pair_type", y="prop", xlab="ORF pair type",ylab ="% of pairs with their promoter bound by a common TF",
          fill = "coexpressed", color=NA, palette = c('#3FC1C9','#E32726'),
          label = FALSE,   legend = "bottom",
          position = position_dodge(0.8))+  rremove("xlab")+
  geom_text(aes(label = format(plot_df$count, big.mark = ",")), size = 2.5, vjust = -0.5,hjust=0.5) +
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "black") +
  stat_pvalue_manual(stats_tf, x='pair_type',label="p.adj.signif",y.position=40,size = significance_size)+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) +
  font("legend.text", size=legend_text_size)+
  font("legend.title", size=legend_title_size)
dev.off()