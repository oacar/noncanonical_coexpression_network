reorder_columns<-function(rho_matrix){
  rho_matrix<-rho_matrix[,rownames(rho_matrix)]
  return(rho_matrix)
}

rho_5_400<-readRDS('rho_raw5_sample400.RDS')
rho_5_400<-reorder_columns(rho_5_400)
num_5_400<-readRDS('numobs_raw5_sample400.RDS')
num_5_400<-reorder_columns(num_5_400)

saveRDS(rho_5_400,'rho_raw5_sample400.RDS')
saveRDS(num_5_400,'numobs_raw5_sample400.RDS')
