library(dplyr)
library(ggplot2)

#orf info including TSS info
translatome<-readRDS('./orf_info_TSS.RDS')

#import coexpression matrix
rho<-readRDS('/home/aar75/coexpression/20221110_rho/spqn_raw5_sample400.RDS')
#import matrix that contains number of samples an orf pair is detected in
num_obs<-readRDS('/home/aar75/coexpression/20221110_rho/numobs_raw5_sample400.RDS')
#set pairs where there are fewer than 400 samples expressing both orfs to NA
rho[num_obs<400]<-NA


#loop through all TF bed files from Rossi et al 
filenames <- list.files(path = 'yep-peaks', pattern="*.bed", full.names=TRUE)
TF_names<-gsub(x = filenames,pattern='^.*peaks/',replacement = "")
TF_names<-gsub(x = TF_names,pattern='\\.multi.*$',replacement = "")

#subset TF matrix to contain only actual TFs
TFs<-c('Yrr1','Zap1','Yap7','Yap5','Yap1','War1','Urc2','Tda9','Swi6','Swi5',
       'Sut1','Sum1','Stp2','Stp1','Ste12','Stb5','Stb4','Rlm1','Rox1','Rpn4','Rtg3',
       'Sfl1','Sfp1','Skn7','Sko1','Sok2','Spt3','Rfx1','Rds2','Rcs1','Pip2','Phd1',
       'Pdr3','Pdr1','Oaf3','Nrm1','Nrg1','Ndd1','Mss11','Mig2','Met4','Met32','Met31','Mcm1',
       'Mbp1','Mac1','Lys14','Leu3','Kar4','Ino4','Ino2','Hsf1','Hot1','Ace2','Aro80','Azf1',
       'Bas1','Cad1','Cha4','Cin5','Crz1','Ecm22','Fkh1','Fkh2','Fzf1','Gal4','Gcn4','Gcr1',
       'Gis1','Gln3','Gzf3','Hal9','Hap1')

idx<-which(TF_names %in% TFs)
TF_names<-TF_names[idx]
filenames<-filenames[idx]

#matrix that is ORFs x TFs and the elements are the min distance to that TFs binding site
TF_distances<-matrix(nrow = nrow(translatome),ncol=length(TF_names))
rownames(TF_distances)<-translatome$transcript
colnames(TF_distances)<-TF_names

for(i in 1:length(filenames))
{
  TF<-read.delim(filenames[i],col.names=c('chr','pos1','pos2','method','score'),header=F,stringsAsFactors=F)
  TF$chr<-as.numeric(gsub(TF$chr,pattern = 'chr',replacement = ''))
  for(j in 1:nrow(translatome))
  {
    
    if(!is.na(translatome$TSS[j])){
      if(translatome$strand[j]=='-')
      {
        x<-filter(TF,chr==translatome$chr_num[j] & pos1>=translatome$TSS[j]) 
        if(nrow(x)>0)
        {
          TF_distances[j,i]<-min(x$pos1)-translatome$TSS[j]
        }
      }else
      {
        x<-filter(TF,chr==translatome$chr_num[j] & pos1<=translatome$TSS[j])
        if(nrow(x)>0)
        {
          TF_distances[j,i]<-translatome$TSS[j]-max(x$pos1)
        }
      }
    }
  }
}


saveRDS(TF_distances,'ORF_TF_distances_TSS_all.RDS')


