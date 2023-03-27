#script to plot the association between being coexpressed and physically interacting and forming a protein complex
#defining coexpressed as rho>99.8th percentile and forming a protein complex from Pu et al (2009)

library(dplyr) 
library(RMariaDB)
library(reshape2)
library(ggplot2)
library(scales)

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
    cor  =(cormat)[ut]
  )
}

complexes_df<-read.delim('CYC2008_complex.tab')


#import coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA

#get cutoff for defining coexpressed (ie rho >99.8th percentile of all rho values)
rf<-flattenCorrMatrix(rho)
cutoff_rf98<-quantile(rf$cor,0.998,na.rm=TRUE)

#keep only annotated genes
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


contingency_998rho<-getComplexEnrichmentRho(cutoff_rf98)

library(ggplot2)
library(dplyr)
library(scales)
library(rstatix)
library(tibble)
library(ggpubr)
rho_hex<-'#DC0000B2'
file_path<-'./20230106_figures/'
significance_size <- 4 # default is 3.88 and does not scale with other size parameters, e.g., element_text(size=4) is different size
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

#convert contingency matrix into a data frame that contains the proportion of coexpressed pairs that form protein complex
#proportion of non-coexpressed pairs that form protein complex

getPlotdf<-function(cm){
  cp = 100*(cm['1','1']/(cm['1','1']+cm['0','1']))
  ncp = 100*(cm['1','0']/(cm['1','0']+cm['0','0']))
  plot_df<-data.frame('coexpressed'=c('coexpressed','not coexpressed'),'percentage'=c(cp, ncp))
  return(plot_df)
}


complex_df<-getPlotdf(contingency_998rho)
complex_stats<-pairwise_fisher_test(contingency_998rho[,1:2]) %>% add_column('coexpressed'='coexpressed')


pdf(sprintf('%sCYC2008_cocomplexEnrichment_rho998.pdf',file_path),width=2.5, height = 3)
complex_df %>% mutate(coexpressed =case_when(coexpressed=='coexpressed pairs' ~'coexpressed',
                                             coexpressed=='not coexpressed pairs' ~'not\ncoexpressed')) %>%
  ggbarplot( x="coexpressed", y="percentage", xlab="cORF pair type",ylab ="% of pairs in a protein complex",
             fill = "coexpressed", color=NA, palette = c('#BFBADA','#E32726'),
             label = FALSE)+
  stat_pvalue_manual(complex_stats, x='coexpressed', label="p.adj.signif",y.position=10, size = significance_size)+
  font("x.text", size = axis_text_size) +
  font("y.text", size = axis_text_size) +
  font("ylab", size = axis_title_size) +
  font("xlab", size = axis_title_size) + grids(axis = c("xy", "x", "y"), color = "grey92", size = NULL, linetype = NULL)+
  rremove("legend")
dev.off()



