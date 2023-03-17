library(tidyverse)
library(GOSemSim)
# loop over all ORFs to find GSEA similarity with their neighbors
# Rscript scripts/gsea_jaccard.r gsea_results_file_name.csv neighbor_dat_file_name.csv output_file_name.csv
args <- commandArgs(trailingOnly = TRUE)
gsea_results <- args[1]
# neighbor_data_file <- args[2]
output_file <- args[2]
gsea_rho <- data.table::fread("data/interim/fgsea_slim_bp_noexponent.csv") # "data/interim/gsea/gsea_output_rho_spqn.csv")

library(RMariaDB)

# import list of genes
con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
# allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList")
allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")
neighbor_df <- DBI::dbGetQuery(con, "select * from omer.pairwise_orf_distances_denovo_conserved_closest_blaste4")
annotatedORFs <- allORFs %>% filter(is_canonical=='canonical')#drop_na(gene)
orf_names <- gsea_rho$transcript[gsea_rho$transcript %in% annotatedORFs$transcript] %>% unique()
res_df <- neighbor_df %>% filter(primary_orf %in% gsea_rho$transcript, neighbor_orf %in% gsea_rho$transcript)
n_iter <- nrow(res_df) # Number of iterations of the loop
# Initializes the progress bar
pb <- txtProgressBar(
    min = 1, # Minimum value of the progress bar
    max = n_iter, # Maximum value of the progress bar
    style = 3, # Progress bar style (also available style = 1 and style = 2)
    width = 50, # Progress bar width. Defaults to getOption("width")
    char = "="
) # Character used to create the bar
d <- godata("org.Sc.sgd.db", ont = "BP", computeIC = TRUE)
for (i in 1:n_iter) {
    orf_name <- res_df$primary_orf[i]
    go_ids <- gsea_rho$pathway[gsea_rho$transcript == orf_name]
    orf_name_j <- res_df$neighbor_orf[i]
    go_ids_j <- gsea_rho$pathway[gsea_rho$transcript == orf_name_j]
    # semsim_wang <- mgoSim(go_ids, go_ids_j, d, measure = "Wang")
    # semsim_resnik <- mgoSim(go_ids, go_ids_j, d, measure = "Resnik")
    # semsim_lin <- mgoSim(go_ids, go_ids_j, d, measure = "Lin")
    semsim_rel <- mgoSim(go_ids, go_ids_j, d, measure = "Rel")
    # semsim_jiang <- mgoSim(go_ids, go_ids_j, d, measure = "Jiang")
    # common <- length(intersect(go_ids, go_ids_j))
    # total <- length(union(go_ids, go_ids_j))
    # jc <- common / total
    results <- c(semsim_rel)
    res_df[i, c("rel")] <- results
    setTxtProgressBar(pb, i)
}


close(pb) # Close the connection

res_df %>% write_csv("data/interim/fgsea_slim_similarity.csv")
res_df %>% DBI::dbWriteTable(con, "fgsea_slim_similarity", ., row.names = FALSE, overwrite = TRUE)


DBI::dbSendQuery(con,"rename table fgsea_slim_bp_noexponent_similarity to fgsea_slim_bp_noexponent_similarity_03052023")
# gsea_rho %>%
#     group_by(transcript) %>%
#     count() %>%
#     # mutate(n = n()) %>%
#     filter(n == 5)
