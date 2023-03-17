suppressMessages(library(tidyverse))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Sc.sgd.db))
# example usage
set.seed(11111)
# Rscript scripts/run_gsea.R coexpression_matrix_file_name.csv pairwise_observation_file_name.csv output_file_name.csv
# args <- commandArgs(trailingOnly = TRUE)

coexp_mat_file <- "/home/oma21/coexpression/data/interim/coexpression_matrix_rho_spqn.csv" #
obs_file <- "/home/oma21/coexpression/data/interim/pairwise_obs.csv" #
# output_file <- args[3]

# print(output_file)


coexpression_matrix <- data.table::fread(coexp_mat_file)
obs <- data.table::fread(obs_file)

mat_matrix <- as.matrix(coexpression_matrix[, 2:ncol(coexpression_matrix)])
obs_matrix <- as.matrix(obs[, 2:ncol(coexpression_matrix)])

mat_matrix[obs_matrix < 400] <- NA

rownames(mat_matrix) <- colnames(mat_matrix)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
# translation <- DBI::dbGetQuery(con, "select * from aaron.scer_orfs_translation_info")
coexpression_orf_list <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")


# df <- data.table::fread(cluster_df_file) #' ../data/interim/wgcna_2/df_6.0_50.csv')
annotated_orfs <- coexpression_orf_list %>%
  filter(is_canonical == "canonical")
nrow(annotated_orfs)
mat_matrix_annotated <- mat_matrix[annotated_orfs$transcript[annotated_orfs$transcript %in% rownames(mat_matrix)], ]

rownames(mat_matrix_annotated) <- annotated_orfs$gene[annotated_orfs$transcript %in% rownames(mat_matrix)]

orf_names <- coexpression_orf_list %>%
  filter(transcript %in% rownames(mat_matrix)) %>%
  dplyr::select(transcript) %>%
  pull()
length(orf_names)
prepGeneList <- function(orf_name, matrix, df) {
  gene_name <- ifelse(df$gene[df$transcript == orf_name] == "X", NA, df$gene[df$transcript == orf_name])
  geneList <- matrix[, orf_name]

  geneList <- geneList[is.na(geneList) == FALSE]
  if (is.na(gene_name) == F) {
    geneList <- geneList[names(geneList) != gene_name]
  }
  geneList <- sort(geneList, decreasing = TRUE)
  geneList
}

n_iter <- length(orf_names) # Number of iterations of the loop

# Initializes the progress bar
pb <- txtProgressBar(
  min = 1, # Minimum value of the progress bar
  max = n_iter, # Maximum value of the progress bar
  style = 3, # Progress bar style (also available style = 1 and style = 2)
  width = 50, # Progress bar width. Defaults to getOption("width")
  char = "="
) # Character used to create the bar

library(fgsea)
go_list <- read.delim("/home/aar75/coexpression/12_20_21_go/go_slim_mapping.tab", stringsAsFactors = F, header = F)
colnames(go_list) <- c("transcript", "common_name", "sgd_id", "go_type", "go_term", "go_id", "orf_type")
# go_list %>% head()
go_list_bp <- go_list[go_list$go_type == "P", ]
go.pathways.from.table <- list()
unique.go.ids <- unique(go_list_bp$go_id)
for (i in 1:length(unique.go.ids)) {
  go.pathways.from.table[[unique.go.ids[i]]] <- go_list$transcript[go_list$go_id == unique.go.ids[i]]
}
GO.pathways <- go.pathways.from.table
results <- list()
results_neg <- list()
exponent <- 1
for (i in seq(1:n_iter)) {
  orf_name <- orf_names[i]
  geneList <- prepGeneList(
    orf_name,
    mat_matrix_annotated, coexpression_orf_list
  )
  fgsea_res <- fgseaMultilevel(GO.pathways, geneList, minSize = 10, nPermSimple = 1000, gseaParam = exponent)
  fgsea_res_neg <- fgsea_res %>%
    filter(pval < 0.05 & NES < 0) %>%
    mutate(transcript = orf_name) %>%
    drop_na()
  fgsea_res_pos <- fgsea_res %>%
    filter(pval < 0.05 & NES > 0) %>%
    mutate(transcript = orf_name) %>%
    drop_na()
  results <- append(results, list(fgsea_res_pos))
  results_neg <- append(results_neg, list(fgsea_res_neg))
  setTxtProgressBar(pb, i)
  # print(i)
}


res_df <- data.table::rbindlist(results)
res_df %>% write_csv("/home/oma21/coexpression/data/interim/fgsea_slim_bp_noexponent.csv")
res_df_neg <- data.table::rbindlist(results_neg)
res_df_neg %>% write_csv("/home/oma21/coexpression/data/interim/fgsea_slim_bp_noexponent_purified.csv")

close(pb) # Close the connection

print("done")


res_df %>% filter(
  transcript == "chr2_614024",
  pathway %in% c(
    "GO:0055085",
    "GO:0006811",
    "GO:0048193",
    "GO:0016192",
    "GO:0016197"
  )
)
coexpression_orf_list %>%
  filter(classification == "denovo", transcript %in% rownames(mat_matrix)) %>%
  nrow()
res_df %>%
  filter(pathway == "GO:0006810") %>%
  nrow()
coexpression_orf_list %>%
  filter(classification == "denovo", transcript %in% rownames(mat_matrix)) %>%
  inner_join(res_df %>% filter(pathway == "GO:0006810"), by = "transcript") %>%
  nrow()
