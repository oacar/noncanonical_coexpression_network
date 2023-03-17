library(tidyverse)
library(RMariaDB)

library(GenomicRanges)
# use zsh for loop to get all the files
# for i in {1..16}; do Rscript scripts/all_by_all_distance.R output_folder/ $i; done
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Usage: Rscript scripts/all_by_all_distance.R <output_dir> <chr_num>")
}

output_dir <- args[1]
chromosome <- args[2]

# import list of genes
con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
allORFs <- DBI::dbGetQuery(con, "select * from omer.coexpressionOrfList_blaste4")

grs <- GenomicRanges::GRanges(seqnames = allORFs$chr_num, strand = allORFs$strand, ranges = IRanges::IRanges(allORFs$coor1, allORFs$coor2,
    names = allORFs$transcript
))
# distance(grs["chr2_89147"],
#     grs["chr2_88523"],
#     ignore.strand = T
# )
# res_list <- list()
i <- chromosome
# i <- 2
# for (i in 1:16) {
print(i)
# chr = i
allORFs_sub <- allORFs %>% filter(chr_num == i)
allORFs_sub_distance <- allORFs_sub %>%
    inner_join(allORFs_sub %>% select(chr_num, transcript, strand, coor1), by = c("chr_num")) %>%
    filter(transcript.x != transcript.y) %>%
    mutate(dist = distance(grs[transcript.x], grs[transcript.y], ignore.strand = T)) %>%
    mutate(
        strands = ifelse(strand.x == strand.y, "same", "different"),
        up_down = case_when(
            strands == "same" & dist == 0 ~ "sense-overlap",
            strands == "different" & dist == 0 ~ "antisense-overlap",
            strand.x == "+" & strand.y == "+" ~ ifelse(coor1.x < coor1.y, "upsame", "downsame"),
            strand.x == "-" & strand.y == "-" ~ ifelse(coor1.x < coor1.y, "downsame", "upsame"),
            strand.x == "+" & strand.y == "-" ~ ifelse(coor1.x < coor1.y, "downopposite", "upopposite"),
            strand.x == "-" & strand.y == "+" ~ ifelse(coor1.x < coor1.y, "upopposite", "downopposite"),
            TRUE ~ "overlap"
        )
    )

#    res_list[[i]] <- allORFs_sub_distance
# }
# allORFs_distances <- data.table::rbindlist(res_list)
# allORFs_distances %>% nrow()
# allORFs_distances %>% colnames()

pairwise_distances <- allORFs_sub_distance %>% select(
    primary_orf_gene_id = gene,
    chr = chr_num, primary_orf = transcript.x, neighbor_orf = transcript.y, primary_orf_strand = strand.x,
    neighbor_orf_strand = strand.y, primary_orf_first_coord = coor1.x,
    primary_orf_second_coord = coor2,
    neighbor_orf_first_coord = coor1.y,
    orf_class, distance = dist, strands, orientation = up_down
)

pairwise_distances %>% write_csv(paste0(output_dir, "/pairwise_orf_distances_", chromosome, ".csv"))

# DBI::dbWriteTable(con, "pairwise_orf_distances", pairwise_distances, overwrite = T)
DBI::dbDisconnect(con)

# gc()
