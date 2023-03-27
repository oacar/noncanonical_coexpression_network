#create a matrix that is ORFxORF and contains a 1 if the two ORFs share at least 1 TF in common and 0 if not
#sharing a TF defined as having a TF bind within 200 bp of an ORFs TSS (binding determined by chip exo bed files from rossi et al)
library(dplyr)
require(Matrix)
require(igraph)
library(scales)
library(parallel)
library(preprocessCore)
library(doParallel)


#import coexpression matrix
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
#keep only orf pairs that have atleast 400 samples expressing both ORFs
rho[num_obs<400]<-NA


getTFsincommon<-function(orf1,orf2)
{
  idxi<-which(orf1==1)
  idxj<-which(orf2==1)
  if(length(intersect(idxi,idxj))>0)
  {
    return(1)
  }else{
    return(0)
  }
}

#matrix that is ORF x TF and contains the min distance of TF binidng to an ORFs TSS
TF_distance<- readRDS('ORF_TF_distances_TSS_all.RDS')

dist_threshold <- 200
#if a TF binds within 200 BP of an ORFs TSS set [ORF,TF] <-1 otherwise set to 0
TF_distance[TF_distance<=dist_threshold]<-1
TF_distance[TF_distance>dist_threshold]<-0
TF_distance[is.na(TF_distance)]<-0




# now create ORF x ORF matrix where the value is 1 if the two ORFs have the same TF binding within 200 bp of their TSS and 0 if not
numCores <- detectCores()
registerDoParallel(numCores-5)

opts <- list(chunkSize=500)
#create tf orf-orf matrix that contains a 1 if the two orfs share atleast 1 TF in common and 0 if not
TFsincommonMatrix <- foreach(i = 1:nrow(TF_distance), .combine = cbind,.options.nws=opts,.packages= c('base','foreach')) %dopar%{
  foreach(j = 1:nrow(TF_distance) , .combine = c, .inorder=TRUE,.packages = c('base','foreach')) %do%{
    
    return(getTFsincommon(TF_distance[i,],TF_distance[j,]))
    
  }
}
colnames(TFsincommonMatrix)<-rownames(TF_distance)
rownames(TFsincommonMatrix)<-rownames(TF_distance)

stopImplicitCluster()

TFsincommonMatrix<-TFsincommonMatrix[rownames(rho),rownames(rho)]
saveRDS(TFsincommonMatrix,'TFsincommonMatrix_TSS_200_all.rds')

