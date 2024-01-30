#neighbor voting like in Skinnider et al https://github.com/skinnider/SCT-MoA/blob/bfd455479d94db92d30153b763d06f5732879606/R/function/calculate-auroc.R#L9

library(dplyr)
library(RMariaDB)
library(ggpubr)
library(EGAD)
#load orf annotation info
#get ORF data

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)


#load coexpression network
rho<-readRDS('spqn_raw5_sample400.RDS')
num_obs<-readRDS('numobs_raw5_sample400.RDS')
rho[num_obs<400]<-NA


go_list<-read.delim('go_slim_mapping.tab',stringsAsFactors = F,header = F)
colnames(go_list)<-c('gene','common_name','sgd_id','go_type','go_term','go_id','orf_type')
go_list<-left_join(go_list, orf_info[,c('gene','transcript','orf_class')])
#filter to keep only biological process terms
go_list<-filter(go_list, orf_class=='Verified' & go_type=='P')

genes_in_common<-intersect(go_list$transcript, colnames(rho))

#filter to only include genes in common (5,133)
go_list<-filter(go_list,transcript %in% genes_in_common)
rho<-rho[genes_in_common, genes_in_common] 



#create go annotation matrix
go_annotations<-make_annotations(go_list[,c('transcript','go_id')],rownames(rho),unique(go_list$go_id))


#from skinnider's github "modified to allow missing values"
n = nrow(rho)
genes = rownames(rho)
net = matrix(rank(rho, na.last = "keep", ties.method = "average"),
             nrow = n, ncol = n)
rownames(net) = colnames(net) = genes
net = net / max(net, na.rm = T)
diag(net) = 1



gba = run_GBA(net, go_annotations) #using default EGAD parameters: min = 20, max = 1000, nfold = 3

aurocs = gba[[1]][, "auc"]
# get number of proteins to which each term was annotated
all_terms = colSums(go_annotations)
n_terms = all_terms[names(aurocs)]
# write result
result = data.frame(term = names(aurocs), auroc = aurocs, 
                    n_proteins = n_terms, 
                    pct_proteins = n_terms / length(genes))


significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


pdf("20231221_rho_canonical_BP_gba_auroc_hist.pdf", width=2,height=2)
gghistogram(result, x="auroc", y="count",fill = "#7570b3",
            xlab='average neighboring voting AUROC across 3-fold cross validation', 
            ylab='GO BP count') + geom_vline(xintercept = 0.5,linetype='dashed')+
  font("x.text", size = axis_text_size)+
  font("y.text",size=axis_text_size)+
  font("ylab",size=axis_title_size)+
  font("ylab",size=axis_title_size)
# font("legend.text", size=legend_text_size)+
# font("legend.title",size=legend_title_size)
dev.off()
