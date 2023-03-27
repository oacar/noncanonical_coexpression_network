#script to plot figure 5H: an example locus where there is a terminator present 
#terminator_info data generated in get_terminator_info.R
#transcript information from TIF seq from Pelechano et al 2013.
library(Rsamtools)
library(ggplot2)
library(dplyr)
library(scales)
library(RMariaDB)
library(gridExtra)
library(cowplot)
library(forcats)


significance_size <- 4 
axis_title_size = 12
legend_text_size = 8
legend_title_size = 10
axis_text_size <- 8


#define genomic coordinate range for plotting
pc1<-639000
pc2<-641800

#get gff info:

conn <- dbConnect(MariaDB(),
                  usr = sql.usr,
                  password = sql.pwd,
                  host = host)
terminator_info<-dbGetQuery(conn, "SELECT * FROM april.Terminator_info")
orf_info<-dbGetQuery(conn, "SELECT * FROM omer.coexpressionOrfList_blaste4")
dbDisconnect(conn)

gff_denovo<-filter(orf_info, transcript == 'chr16_641385')
gff_denovo$title<-'denovo'

gff_conserved<-filter(orf_info, transcript== 'chr16_639525')
gff_conserved$title<-'conserved'

gff<-rbind(gff_conserved[,c('strand','coor1','coor2','title')],gff_denovo[,c('strand','coor1','coor2','title')])
#get tif seq info:

library(data.table)
tif<-fread('./tif_seq/S1_TIFs.txt',fill=TRUE)
tif<-as.data.frame(tif)

idx<-which(tif$strand=='-')
temp<-tif$t5[idx]
tif$t5[idx]<-tif$t3[idx]
tif$t3[idx]<-temp
tif$tot_count<-tif$ypd+tif$gal

#get all transcripts in this genomic window
tif<-filter(tif,chr==gff_denovo$chr_num[1] & ((t3>pc1 & t3<pc2)|(t5>pc1 & t5<pc2)|(t5<pc1 & t3>pc2)) & strand==gff_denovo$strand[1])


colnames(tif)[c(1,2,3,4)]<-c('chr_num','strand','coor1','coor2')
tif$length<-tif$coor2-tif$coor1
idx<-which(tif$coor1 < pc1)
tif$coor1[idx]<- pc1
idx<-which(tif$coor2 > pc2)
tif$coor2[idx]<- pc2
tif$title<-'TIF-seq'

tif<-tif[with(tif, order(-coor2,-length)),]
tif$position<-rev((1:nrow(tif)))/10

#get terminator info
#import insulator chip exo data
importChipExoData<-function(protein_name){
  filename<-sprintf('./yep-peaks/%s.multi_%s.filtered.bed',protein_name,protein_name)
  data<-read.delim(filename,
                   header = FALSE,col.names=c('chr','pos1','pos2','method','score'))
  data<-filter(data,score>=43.2)
  data$chr<-gsub(data$chr,pattern = 'chr',replacement = '')
  data$chr<-as.integer(data$chr)
  return(data)
}
pcf11_data<-importChipExoData('Pcf11')
x<-filter(pcf11_data, chr == gff_denovo$chr_num[1] & pos1 < pc2 & pos1 >pc1)

preads1<-ggplot(tif)+ geom_segment(mapping=aes(x=coor1, xend=coor2, y=position,yend=position),color="#404142")+
  geom_vline(xintercept = x$pos1, color="red",linetype="dotted",linewidth=1)+ theme_minimal()+ 
  scale_x_continuous(limits=c(pc1-1, pc2+1),label=comma)+coord_cartesian(xlim=c(pc1-1, pc2+1), expand=F)+
  labs(y="TIF-seq reads",x=NULL)+theme(axis.text.y = element_blank(), axis.text.x=element_blank(), axis.title.y = element_text(size=axis_title_size))#+theme(axis.text = element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm"))

pgff<-ggplot(gff, aes(fill=title))+ geom_rect(mapping=aes(xmin=coor1, xmax=coor2, ymin=0.2,ymax=0.35))+theme_minimal()+ 
  geom_vline(xintercept = x$pos1, color="red",linetype="dotted",linewidth=1)+
  scale_x_continuous(limits=c(pc1-1, pc2+1),label=comma)+coord_cartesian(xlim=c(pc1-1, pc2+1), expand=F)+
  scale_y_continuous(breaks = NULL)+
  labs(y="",x="genomic coordinates")+theme(axis.text.y = element_blank(),
                                           axis.text.x= element_text(size=axis_text_size) ,
                                           axis.title.x = element_text(size=axis_title_size),
                                           legend.position = "none")+ scale_fill_manual(values=c('#F8766D','#619CFF')) 


pdf('./20230106_figures/loci_chr16_641385.pdf',width=4,height=3.5)
print(plot_grid(preads1,pgff,align = "v", axis = "lr",nrow=2, rel_heights = c(7/10,3/10)))
dev.off()
