#see if terminators are present between conserved genes and neighboring de novo ORFs using chip exo
# binding info from nrd1 and pcf11 (from rossi et al 2021)

library(dplyr)
library(ggplot2)
library(RMariaDB)
library(data.table)


conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_pairs<-dbGetQuery(conn, "SELECT * FROM omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
orf_info_df<-dbGetQuery(conn,'SELECT * FROM omer.coexpressionOrfList_blaste4')
dbDisconnect(conn)

orf_pairs<-filter(orf_pairs,distance<=1000)

#import terminator chip exo data 
#chip exo bed files from rossi et al 2021
importChipExoData<-function(protein_name){
  filename<-sprintf('./yep-peaks/%s.multi_%s.filtered.bed',protein_name,protein_name)
  data<-read.delim(filename,
                   header = FALSE,col.names=c('chr','pos1','pos2','method','score'))
  data<-filter(data,score>=43.2)
  data$chr<-gsub(data$chr,pattern = 'chr',replacement = '')
  data$chr<-as.integer(data$chr)
  return(data)
}
nrd1_data<-importChipExoData('Nrd1')
pcf11_data<-importChipExoData('Pcf11')





findbindingUpstreamSame<-function(gene_info,orf_info,chipExo){
  if(orf_info$strand == '+'){
    x<-filter(chipExo, chr == orf_info$chr_num & pos1 >= orf_info$coor2 & pos1<= gene_info$coor1)
  }else{
    #on - strand 
    x<-filter(chipExo, chr== orf_info$chr_num & pos1 >= gene_info$coor2 & pos1 <= orf_info$coor1)
  }
  return(nrow(x))
}

findbindingUpstreamOpp<-function(gene_info,orf_info,chipExo){
  if(orf_info$strand == '+'){
    x<-filter(chipExo, chr==orf_info$chr_num & pos1 >= gene_info$coor2 & pos1 <= orf_info$coor1)
    
  }else{
    x<-filter(chipExo, chr==orf_info$chr_num & pos1 >= orf_info$coor2 & pos1 <= gene_info$coor1)
  }
  return(nrow(x))
}

findbindingDownstreamSame<-function(gene_info,orf_info,chipExo){
  if(orf_info$strand == '+'){
    x<-filter(chipExo,chr==orf_info$chr_num & pos1 >= gene_info$coor2 & pos1 <=orf_info$coor1)
  }else{
    x<-filter(chipExo,chr==orf_info$chr_num & pos1 >= orf_info$coor2 & pos1 <= gene_info$coor1)
  }
  return(nrow(x))
}

findbindingDownstreamOpp<-function(gene_info,orf_info,chipExo){
  if(orf_info$strand == '+'){
    x<-filter(chipExo,chr==orf_info$chr_num & pos1>= orf_info$coor2 & pos1 <=gene_info$coor1)
  }else{
    x<-filter(chipExo,chr==orf_info$chr_num & pos1 >= gene_info$coor2 & pos1 <=orf_info$coor1)
  }
  return(nrow(x))
}


orf_pairs$nrd1<-NA
orf_pairs$pcf11<-NA
for(i in 1:nrow(orf_pairs)){
  gene<-filter(orf_info_df, transcript== orf_pairs$neighbor_orf[i])
  orf<-filter(orf_info_df,transcript == orf_pairs$primary_orf[i])
  if(orf_pairs$orientation[i]=='upsame'){
    orf_pairs$nrd1[i]<-findbindingUpstreamSame(gene,orf,nrd1_data)
    orf_pairs$pcf11[i]<-findbindingUpstreamSame(gene,orf,pcf11_data)
  }else if(orf_pairs$orientation[i]=='upopposite'){
    orf_pairs$nrd1[i]<-findbindingUpstreamOpp(gene,orf,nrd1_data)
    orf_pairs$pcf11[i]<-findbindingUpstreamOpp(gene,orf,pcf11_data)
  }else if(orf_pairs$orientation[i]=='downsame'){
    orf_pairs$nrd1[i]<-findbindingDownstreamSame(gene,orf,nrd1_data)
    orf_pairs$pcf11[i]<-findbindingDownstreamSame(gene,orf,pcf11_data)
  }else if(orf_pairs$orientation[i]=='downopposite'){
    orf_pairs$nrd1[i]<-findbindingDownstreamOpp(gene,orf,nrd1_data)
    orf_pairs$pcf11[i]<-findbindingDownstreamOpp(gene,orf,pcf11_data)
  }else{}
}

saveRDS(orf_pairs,'terminator_binding.RDS')





