#script that calculates the association between an ORF pair being coexpressed and having their promoters bound by a common TF. 
#coexpression network method used is RPKM + pearson correlation
#this script generates the contingency matrices at many different coexp cutoffs, and this will be used to compare to other methods for constructing coexpression 

#coexpression_pearson_rpkm generated in get_RPKM_pearson_coexpression.R
#TF_distance data generated in: get_TF_binding_distance_matrix.R
#TFsincommon data generated in: get_TFs_in_common.R

library(dplyr)
library(RMariaDB)

#load orf annotation info
#get ORF data

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = 'paris.csb.pitt.edu')
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


coexpression_pearson_rpkm<-readRDS('rpkm_raw5_sample400_coexp_pearson.RDS')

#do TF enrichment
flattenCorrMatrix <- function(cormat){#, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut]#,
    #p = pmat[ut]
  )
}

#get rho cutoffs before subsetting the network

rf<-flattenCorrMatrix(coexpression_pearson_rpkm)
cutoff_rf999<-quantile(rf$cor,0.999,na.rm=TRUE)
cutoff_rf998<-quantile(rf$cor,0.998,na.rm=TRUE)
cutoff_rf99<-quantile(rf$cor, 0.99, na.rm=TRUE)
cutoff_rf95<-quantile(rf$cor,0.95,na.rm=TRUE)
cutoff_rf90<-quantile(rf$cor, 0.90, na.rm=TRUE)

#subset to contain only orfs that have atleast 1 TFBS
TF_distance<-readRDS('ORF_TF_distances_TSS_all.RDS')

TF_distance[TF_distance<=200]<-1
TF_distance[TF_distance>200]<-0
TF_distance[is.na(TF_distance)]<-0

TF_count<-rowSums(TF_distance,na.rm=TRUE)
orfsWithTFBS<-names(which(TF_count>0))

orfsWithTFBS<-intersect(orfsWithTFBS,rownames(coexpression_pearson_rpkm))


#how many nORFs with TFBS? how many cORFs with TFBS?
length(orfsWithTFBS) #1909
x<-left_join(data.frame(transcript=orfsWithTFBS), orf_info)
table(x$is_canonical) 
#   canonical noncanonical
#       973          936

coexpression_pearson_rpkm<-coexpression_pearson_rpkm[orfsWithTFBS,orfsWithTFBS]

#import TFs in common matrix 
#read in TFs in common matrix (orf x orf matrix where 0 corresponds to having no TFs in common and value of 1
# means that orf pairs share atleast 1 TF in common.) #using 200 bp upstream of TSS for cutoff 
TFsincommon<-readRDS('TFsincommonMatrix_TSS_200_all.rds')



TFsincommon<-TFsincommon[orfsWithTFBS,orfsWithTFBS]


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

cm_pearson_90_cc<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf90,orf_info,comparison_type='cc')
saveRDS(cm_pearson_90_cc,'pearson_rpkm_cc_pairswithTFBS_200_90.RDS')
cm_pearson_95_cc<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf95,orf_info,comparison_type='cc')
saveRDS(cm_pearson_95_cc,'pearson_rpkm_cc_pairswithTFBS_200_95.RDS')
cm_pearson_99_cc<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf99,orf_info,comparison_type='cc')
saveRDS(cm_pearson_99_cc,'pearson_rpkm_cc_pairswithTFBS_200_99.RDS')
cm_pearson_998_cc<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf998,orf_info,comparison_type='cc')
saveRDS(cm_pearson_998_cc,'pearson_rpkm_cc_pairswithTFBS_200_998.RDS')
cm_pearson_999_cc<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf999,orf_info,comparison_type='cc')
saveRDS(cm_pearson_999_cc,'pearson_rpkm_cc_pairswithTFBS_200_999.RDS')


cm_pearson_90_nn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf90,orf_info,comparison_type='nn')
saveRDS(cm_pearson_90_nn,'pearson_rpkm_nn_pairswithTFBS_200_90.RDS')
cm_pearson_95_nn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf95,orf_info,comparison_type='nn')
saveRDS(cm_pearson_95_nn,'pearson_rpkm_nn_pairswithTFBS_200_95.RDS')
cm_pearson_99_nn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf99,orf_info,comparison_type='nn')
saveRDS(cm_pearson_99_nn,'pearson_rpkm_nn_pairswithTFBS_200_99.RDS')
cm_pearson_998_nn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf998,orf_info,comparison_type='nn')
saveRDS(cm_pearson_998_nn,'pearson_rpkm_nn_pairswithTFBS_200_998.RDS')
cm_pearson_999_nn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf999,orf_info,comparison_type='nn')
saveRDS(cm_pearson_999_nn,'pearson_rpkm_nn_pairswithTFBS_200_999.RDS')

cm_pearson_90_cn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf90,orf_info,comparison_type='cn')
saveRDS(cm_pearson_90_cn,'pearson_rpkm_cn_pairswithTFBS_200_90.RDS')
cm_pearson_95_cn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf95,orf_info,comparison_type='cn')
saveRDS(cm_pearson_95_cn,'pearson_rpkm_cn_pairswithTFBS_200_95.RDS')
cm_pearson_99_cn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf99,orf_info,comparison_type='cn')
saveRDS(cm_pearson_99_cn,'pearson_rpkm_cn_pairswithTFBS_200_99.RDS')
cm_pearson_998_cn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf998,orf_info,comparison_type='cn')
saveRDS(cm_pearson_998_cn,'pearson_rpkm_cn_pairswithTFBS_200_998.RDS')
cm_pearson_999_cn<-getContingencyMatrix(TFsincommon,coexpression_pearson_rpkm,cutoff_rf999,orf_info,comparison_type='cn')
saveRDS(cm_pearson_999_cn,'pearson_rpkm_cn_pairswithTFBS_200_999.RDS')

