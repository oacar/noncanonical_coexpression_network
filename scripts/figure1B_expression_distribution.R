library(dplyr)
library(RMariaDB)
library(ggplot2)
library(scales)
#path to save plots
file_path<-'./20230106_figures/'

axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8

#read in raw rna seq counts 
raw_data<-readRDS('raw_counts.RDS')

#get length information about each ORF
conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)
orf_info$length<-orf_info$coor2-orf_info$coor1

#get number of samples an ORF is detected in (defining detected as raw >5)
plot_df<-data.frame('transcript'=rownames(raw_data),
                    'sample_count'=apply(raw_data,1,function(x){length(which(x>5))}),stringsAsFactors = F)

#set raw counts less than five to zero
raw_data[raw_data<=5]<-0

#calculate TPM 
length_df<-data.frame('transcript'=rownames(raw_data))
length_df<-left_join(length_df,orf_info[,c('transcript','length')])
lengths_matrix<-matrix(length_df$length, nrow=length(length_df$length), ncol=ncol(raw_data), byrow=FALSE)
lengths_matrix<-lengths_matrix/1000
rpk<-raw_data/lengths_matrix
scaling<-colSums(rpk)/1E6
scaling_matrix<-matrix(scaling, nrow=nrow(raw_data),ncol=length(scaling), byrow=TRUE)
tpm<-rpk/scaling_matrix
tpm[raw_data==0]<-NA

#get median tpm for each ORF across all samples 
median_tpm<-apply(tpm,1,median,na.rm=TRUE)
median_tpm<-data.frame('transcript'=names(median_tpm),'tpm'=median_tpm)
median_tpm<-left_join(median_tpm, plot_df[,c('transcript','sample_count')])
median_tpm<-left_join(median_tpm,orf_info[,c('transcript','is_canonical')])


rescale(c(0,50,100,500,1000), to =c(0,1))
pdf(sprintf("%sfigure1_expressionHexPlot_canonicalNoncanonical.pdf", file_path), width = 6, height = 3)
median_tpm %>%
  ggplot(aes(x=sample_count,y=tpm))+geom_bin_2d(bins = 15)+ #ylim(c(0,10))+
  scale_fill_distiller(palette = "Spectral",name='ORF count',#values = rescale(c(25,50,100,500,1000,3000), to =c(0,1)),
                       breaks=c(10,100,1000),
                       trans="log",
                       guide = guide_colorbar(barwidth = 0.8, barheight = 10))+
  facet_wrap(~is_canonical)+ scale_y_log10()+ 
  #scale_y_continuous(trans = pseudo_log_trans(base= 10),breaks=c(1,10,100,1000))+
  theme_minimal()+theme(legend.text=element_text(size=legend_text_size),
                        legend.title=element_text(size=legend_title_size),
                        axis.title = element_text(size=axis_title_size),axis.text=element_text(size=axis_text_size))+ 
  labs(x='sample count', y='median expression (tpm)')
dev.off()


