#get transitivity and diameter for 1000 randomized networks 
library(igraph)
library(dplyr)
#folder containing 1000 randomized networks
networks<-list.files(path = "./randomized_networks",full.names = TRUE)
for(f in networks){
  network<-read.csv(f, header=FALSE)
  colnames(network)[2:3]<-c('row','column')
  g<-graph_from_data_frame(network[,c('row','column')],directed = FALSE) #keep only edges before doing this
  
  d<-diameter(g)
  tt<-transitivity(g)
  dd<-mean_distance(g)
  
  filepath<-"diameter_random.txt"
  write(d,file=filepath,append=TRUE)
  
  filepath<-"transitivity_random.txt"
  write(tt,file=filepath,append=TRUE)
  
}

