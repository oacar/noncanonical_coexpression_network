#!/usr/bin/env Rscript

#get all samples names from the names of the directories in the counts folder
#then extract the TPM column from each sample's quant.sf file 
setwd("counts/SE/")
dirs<-list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
x<-read.delim(paste(dirs[1],'quant.sf',sep = '/'),stringsAsFactors = FALSE)

RAW<-x[,c(1,5)]

colnames(RAW)[2]<-gsub(dirs[1],pattern = "./",replacement = '')

for (i in 2:length(dirs))
  #for(i in 2:10)
{
  if(file.exists(paste(dirs[i],'quant.sf',sep = '/')))
  {
    temp<-read.delim(paste(dirs[i],'quant.sf',sep = '/'),stringsAsFactors = FALSE)
    if(all(temp$Name!=RAW$Name))
    {
      print(dirs[i])
      print('not the same transcript names')
    }
   
    RAW<-cbind(RAW,temp$NumReads)
    colnames(RAW)[i+1]<-gsub(dirs[i],pattern = "./",replacement = '')
   
    
  }
  else
  {
    print(dirs[i])
  }
}


#clean up transcript names

RAW$Name<-gsub(pattern = '::.*$',replacement = '',x = RAW$Name)
write.csv(RAW,'RAW_counts_SE.csv',row.names=FALSE)
