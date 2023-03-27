library(RMariaDB)
library(dplyr)
library("Biostrings")
library(stringr)

#get orf and gene list
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
coexpression_orf_list<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

coexpression_orf_list$length<-coexpression_orf_list$coor2-coexpression_orf_list$coor1

#filter to keep orfs that are canonical OR have evidence of translation 
coexpression_orf_list<- filter(coexpression_orf_list, orf_class %in% c('Verified','Uncharacterized','transposable_element_gene') | qval<0.05) #26869
coexpression_orf_list<-filter(coexpression_orf_list, length>=25) #25,724

#make bed file and then use bedtools to get list of orfs to remove overlapping orfs
NCBI_chrID<-c('ref|NC_001133|','ref|NC_001134|','ref|NC_001135|','ref|NC_001136|','ref|NC_001137|','ref|NC_001138|',
              'ref|NC_001139|','ref|NC_001140|','ref|NC_001141|','ref|NC_001142|','ref|NC_001143|','ref|NC_001144|',
              'ref|NC_001145|','ref|NC_001146|','ref|NC_001147|','ref|NC_001148|')
coexpression_orf_list$chr<-NA
for(i in 1:length(NCBI_chrID))
{
  idx<-which(coexpression_orf_list$chr_num==i)
  coexpression_orf_list$chr[idx]<-NCBI_chrID[i]
}


idx<-which(coexpression_orf_list$strand=='+')
coexpression_orf_list$coor1[idx]<-coexpression_orf_list$coor1[idx]-1
bed<-data.frame('chr'=coexpression_orf_list$chr,'start'=coexpression_orf_list$coor1,
                'end'=coexpression_orf_list$coor2,'name'=coexpression_orf_list$transcript,'dot'=1,
                'strand'=coexpression_orf_list$strand, stringsAsFactors = FALSE)
write.table(bed,file = 'all.bed',quote = FALSE,row.names = FALSE,sep="\t",col.names = FALSE)
