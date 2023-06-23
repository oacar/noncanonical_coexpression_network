#supplementary figure 4
#script to plot the association between being coexpressed and physically interacting and forming a protein complex
#using many coexpression cutoffs to define coexpressed
#defined forming a protein complex using complex data from Pu et al (2009)

library(dplyr) 
library(RMariaDB)
library(reshape2)
library(ggplot2)
library(scales)
library(rstatix)
library(tibble)
library(ggpubr)

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


complexes_df<-read.delim('CYC2008_complex.tab')

#import rho matrix
#import coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA


rf<-flattenCorrMatrix(rho)
cutoff_rf998<-quantile(rf$cor,0.998,na.rm=TRUE)
cutoff_rf999<-quantile(rf$cor,0.999,na.rm=TRUE)
cutoff_rf99<-quantile(rf$cor,0.99,na.rm=TRUE)
cutoff_rf95<-quantile(rf$cor,0.95,na.rm=TRUE)
cutoff_rf90<-quantile(rf$cor,0.90,na.rm=TRUE)



orfsInRho<-filter(orf_info, transcript %in% colnames(rho) & !is.na(gene))
genesInCommon<-intersect(orfsInRho$gene,complexes_df$ORF)

# #remove genes from complex data not in coexpression data 
complexes_df<-inner_join(complexes_df,orfsInRho[,c('gene','transcript')],by=c('ORF'='gene'))

#remove orfs from coexpression network not in complex data for rho
idx<-which(rownames(rho) %in% complexes_df$transcript)
rho<-rho[idx,idx]
all(rownames(rho) %in% complexes_df$transcript)

getComplexEnrichmentRho<-function(cutoff){
  #initialize zero matrix to store protein complex information 
  complexes<-array(0L, dim(rho))
  colnames(complexes)<-colnames(rho)
  rownames(complexes)<-rownames(rho)
  
  #loop through each matrix to fill in complex association
  for(i in 1:nrow(complexes))
  {
    #for each ORF (ie row) in complexes matrix, find orf in complexes_df
    #get name of complex ORF is a part of 
    #get all unique ORFs also involved in that complex
    #find index associated with these other orfs (idx)
    #set [i,idx]<-1 to indicate complex association 
    orf<-rownames(complexes)[i]
    complex_i<-filter(complexes_df, transcript==orf)$Complex
    other_ORFS<-unique(filter(complexes_df,Complex %in% complex_i)$transcript)
    complexes[i,other_ORFS]<-1
  }
  
  #flatten coexpression matrix
  rho_flat<-flattenCorrMatrix(rho)
  #flatten complex matrix
  complexes_flat<-flattenCorrMatrix(complexes)
  colnames(complexes_flat)[3]<-'complex'
  #combine two lists, so that its orf1, orf2, rho, complex (1/0)
  df<-left_join(rho_flat,complexes_flat)
  df<-filter(df, !is.na(cor))
  
  df<-mutate(df, coexpressed= case_when(cor >cutoff ~1, TRUE ~ 0))
  
  #get contingency matrix
  contingency_table<-table(df$complex,df$coexpressed)
  return(contingency_table)
}


contingency_998rho<-getComplexEnrichmentRho(cutoff_rf998)
write.csv(contingency_998rho,'CYC2008_cocomplex_enrichment_rho_998th.csv')

contingency_999rho<-getComplexEnrichmentRho(cutoff_rf999)
write.csv(contingency_999rho,'CYC2008_cocomplex_enrichment_rho_999th.csv')


contingency_99rho<-getComplexEnrichmentRho(cutoff_rf99)
write.csv(contingency_99rho,'CYC2008_cocomplex_enrichment_rho_99th.csv')

contingency_95rho<-getComplexEnrichmentRho(cutoff_rf95)
write.csv(contingency_95rho,'CYC2008_cocomplex_enrichment_rho_95th.csv')

contingency_90rho<-getComplexEnrichmentRho(cutoff_rf90)
write.csv(contingency_90rho,'CYC2008_cocomplex_enrichment_rho_90th.csv')




file_path<-'./20230106_figures/'
significance_size <- 4 
axis_title_size = 12
legend_text_size = 10
legend_title_size = 12
axis_text_size <- 10

getPlotdf_csv<-function(cm){
  cm$X<-NULL
  colnames(cm)<-c('0','1')
  rownames(cm)<-c('0','1')
  cp = (cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = (cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed pairs','not coexpressed pairs'),'prop'=c(cp, ncp),'count'=c(cm['1','1'],cm['1','0']))
  plot_df$se<-sqrt(plot_df$prop * (1 - plot_df$prop) / plot_df$count)
  plot_df$prop<-plot_df$prop*100
  plot_df$se<-plot_df$se*100
  return(plot_df)
}
getPlotdf_rds<-function(cm){
  cp = (cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = (cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed','not coexpressed'),'prop'=c(cp, ncp),'count'=c(cm['1','1'],cm['1','0']))
  plot_df$se<-sqrt(plot_df$prop * (1 - plot_df$prop) / plot_df$count)
  plot_df$prop<-plot_df$prop*100
  plot_df$se<-plot_df$se*100
  return(plot_df)
}

contingency_matrix_95<-read.csv('CYC2008_cocomplex_enrichment_rho_95th.csv')
complex_df_95<-getPlotdf_csv(contingency_matrix_95)
complex_stats_95<-pairwise_fisher_test(contingency_matrix_95[,2:3]) %>% add_column('coexpressed'='coexpressed')
complex_df_95$coexp_cutoff<-0.95
complex_stats_95 <- complex_stats_95 %>% add_column ('coexp_cutoff'=0.95)

contingency_matrix_90<-read.csv('CYC2008_cocomplex_enrichment_rho_90th.csv')
complex_df_90<-getPlotdf_csv(contingency_matrix_90)
complex_stats_90<-pairwise_fisher_test(contingency_matrix_90[,2:3]) %>% add_column('coexpressed'='coexpressed')
complex_stats_90 <- complex_stats_90 %>% add_column ('coexp_cutoff'=0.90)
complex_df_90$coexp_cutoff<-0.90

contingency_matrix_99<-read.csv('CYC2008_cocomplex_enrichment_rho_99th.csv')
complex_df_99<-getPlotdf_csv(contingency_matrix_99)
complex_stats_99<-pairwise_fisher_test(contingency_matrix_99[,2:3]) %>% add_column('coexpressed'='coexpressed')
complex_stats_99 <- complex_stats_99 %>% add_column ('coexp_cutoff'=0.99)
complex_df_99$coexp_cutoff<-0.99

contingency_matrix_999<-read.csv('CYC2008_cocomplex_enrichment_rho_999th.csv')
complex_df_999<-getPlotdf_csv(contingency_matrix_999)
complex_stats_999<-pairwise_fisher_test(contingency_matrix_999[,2:3]) %>% add_column('coexpressed'='coexpressed')
complex_stats_999 <- complex_stats_999 %>% add_column ('coexp_cutoff'=0.999)
complex_df_999$coexp_cutoff<-0.999

plot_df<-rbind(complex_df_90,complex_df_95)
plot_df_stats<-rbind(complex_stats_90,complex_stats_95)

plot_df<-rbind(plot_df,complex_df_99)
plot_df_stats<-rbind(plot_df_stats, complex_stats_99)

plot_df<-rbind(plot_df,complex_df_999)
plot_df_stats<-rbind(plot_df_stats, complex_stats_999)


p<-plot_df %>% mutate(coexpressed =case_when(coexpressed=='coexpressed pairs' ~'coexpressed',
                                             coexpressed=='not coexpressed pairs' ~'not\ncoexpressed')) %>%
  ggbarplot( x="coexpressed", y="prop", xlab="cORF pair type",ylab ="% of pairs in a protein complex",
             fill = "coexpressed", color=NA, palette = c('#3FC1C9','#E32726'), 
             label = plot_df$count) + 
  # geom_text(aes(y = ..count.., group = interaction(coexpressed, coexp_cutoff)), 
  #           size = 2.5, vjust = -0.5, hjust = 0.5, position = position_dodge(width = .9)) +
  # geom_text(aes(label = format(plot_df_stats$count, big.mark = ",")), size = 2.5, vjust = -0.5,hjust=0.5) +
  stat_pvalue_manual(plot_df_stats, x='coexpressed', label="p.adj.signif",y.position=9.5, size = significance_size)+
  geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2, color = "black") +
  font("x.text", size = 9) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + grids(axis = c("xy", "x", "y"), color = "grey92", size = NULL, linetype = NULL)+
  rremove("legend")
pdf(sprintf('%sCYC2008_cocomplexEnrichment_supplementary_figure.pdf',file_path),width=7.5, height = 3.75)
facet(p,facet.by="coexp_cutoff" ,nrow = 1)
dev.off()

