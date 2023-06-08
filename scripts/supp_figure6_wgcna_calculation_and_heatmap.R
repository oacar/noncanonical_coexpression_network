# wcgna_try.R
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(tidygraph))
suppressMessages(library(ggraph))
suppressMessages(library(WGCNA))
suppressMessages(library(ComplexHeatmap))

args <- commandArgs(trailingOnly = TRUE)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
respectD <- TRUE
dynamic <- F
if (length(args) == 0) {
  mat_file <- "/home/oma21/coexpression/data/interim/coexpression_matrix_rho_spqn.csv"
  obs_file <- "/home/oma21/coexpression/data/interim/pairwise_obs.csv"
  threshold <- 12 # as.numeric(args[3])
  tree_filename <- NA # args[4]
  heatmap_filename <- NA # args[5]
  size_distribution_plot <- NA # args[6] depreceated
  adjacency_file <- NA # args[7]
  df_filename <- NA # args[8]
  min_module_size <- 50 # as.numeric(args[9])
} else {
  mat_file <- args[1]
  obs_file <- args[2]
  threshold <- as.numeric(args[3])
  tree_filename <- args[4]
  heatmap_filename <- args[5]
  # size_distribution_plot <- args[6] # depreceated
  adjacency_file <- args[6]
  df_filename <- args[7]

  min_module_size <- as.numeric(args[8])
}
print(args)
allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")
mat <- data.table::fread(mat_file)
obs <- data.table::fread(obs_file)

mat_matrix <- as.matrix(mat[, 2:ncol(mat)])
obs_matrix <- as.matrix(obs[, 2:ncol(mat)])

mat_matrix[obs_matrix < 400] <- 0
rownames(obs_matrix) <- colnames(obs_matrix)
rownames(mat_matrix) <- colnames(mat_matrix)
translated_in_network <- allORFs$transcript[allORFs$transcript %in% colnames(mat_matrix)]

allORFs %>%
  filter(transcript %in% colnames(mat_matrix)) %>%
  select(classification) %>%
  table()
mat.matrix.translated <- mat_matrix[translated_in_network, translated_in_network]

print("Matrix created")
print(dim(mat.matrix.translated))
# Enable parallel backend
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
# Call the network topology analysis function

soft_thr <- threshold
adj <- adjacency.fromSimilarity(mat.matrix.translated, power = soft_thr)
if (is.na(adjacency_file) == F) {
  write.table(adj, adjacency_file, quote = FALSE, sep = ",")
}
print("Adj created")
print(dim(adj))
TOM <- TOMsimilarity(adj)
dissTOM <- 1 - TOM

print("TOM calculated")

# Call the hierarchical clustering function
geneTree <- hclust(as.dist(dissTOM), method = "ward.D2")

# We like large modules, so we set the minimum module size relatively high:
# Module identification using dynamic tree cut:
if (dynamic) {
  dynamicMods <- cutreeDynamic(
    dendro = geneTree, distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = respectD, minClusterSize = min_module_size
  )
} else {
  dynamicMods <- cutreeDynamic(
    dendro = geneTree, distM = dissTOM, method = "tree",
    deepSplit = TRUE,
    minClusterSize = min_module_size
  )
  #                            pamRespectsDendro = respectD, minClusterSize = minModuleSize);
}


print("Tree was cut")
# Convert numeric lables into colors
dynamicColors <- labels2colors(dynamicMods)
# Plot the dendrogram and colors underneath
if (is.na(tree_filename) == F) {
  png(tree_filename, width = 8, height = 8, units = "in", res = 150)
  plotDendroAndColors(geneTree, dynamicColors, "Clusters", dendroLabels = FALSE, hang = 0.03, addGuide = T, guideHang = 0.05, main = "")
  dev.off()
}
print(paste0("There are ", length(unique(dynamicMods)), " modules"))
orf_cluster <- c(dynamicMods)
# names(orf_cluster) <- geneTree$order


df_w_cls <- data.frame(orf_cluster) %>%
  rownames_to_column("idx") %>%
  mutate(idx = as.numeric(idx)) %>%
  data.frame(transcript = colnames(mat.matrix.translated), idx = 1:nrow(mat.matrix.translated)) %>%
  left_join(allORFs, by = c("transcript" = "transcript"))
print(dim(df_w_cls))
if (is.na(df_filename) == F) {
  df_w_cls %>%
    mutate(min_size = min_module_size) %>%
    readr::write_csv(df_filename)
}
# annot_color <- c("other" = "green", "denovo" = "blue", "conserved" = "orange")
annot_color <- c("canonical" = "#7570b3", "noncanonical" = "#1b9e77")
names(dynamicColors) <- dynamicMods
row_ha <- rowAnnotation(
  df = data.frame(
    classification = df_w_cls$is_canonical,
    cluster = df_w_cls$orf_cluster
  ),
  col = list(classification = annot_color, cluster = dynamicColors)
)
library(circlize)
col_fun <- colorRamp2(c(0, 0.3), c("white", "red"))
hm <- Heatmap(adj,
  col = col_fun,
  cluster_rows = as.dendrogram(geneTree),
  cluster_columns = as.dendrogram(geneTree),
  right_annotation = row_ha,
  # row_dend_width = unit(1, "in"),
  # column_dend_height = unit(1, "in"),
  # show_column_dend = F,
  # show_row_dend = F,
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = T,
  raster_resize_mat = max,
  raster_device = "png"
)

print("Plotting heatmap")
if (is.na(heatmap_filename) == F) {
  pdf(heatmap_filename, width = 12, height = 12, units = "in", res = 150) # , units = "in", res = 150)
  pdf("reports/figures/paper_figures_01062023/heatmap_with_dend.pdf", width = 12, height = 12) # , units = "in", res = 150)
  hm
  dev.off()
}


allORFs %>%
  filter(transcript %in% colnames(mat_matrix) & classification == "denovo") %>%
  count(is_canonical) %>%
  ggplot(aes(x = is_canonical, fill = is_canonical, y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = annot_color) +
  theme_bw() +
  labs(x = "", y = "Number of de novo ORFs") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 8), axis.title = element_text(size = 8), axis.text = element_text(size = 8), plot.title = element_text(size = 8))
ggsave("reports/figures/paper_figures_01062023/denovo_canonical_barplot.pdf", width = 1.5, height = 3)
