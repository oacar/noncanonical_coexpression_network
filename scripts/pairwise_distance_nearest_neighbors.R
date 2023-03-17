library(tidyverse)

con <- DBI::dbConnect(RMariaDB::MariaDB(), groups = "mariaDB")
allORFs <- DBI::dbGetQuery(con, "select * from coexpressionOrfList_blaste4")

prev <- DBI::dbGetQuery(con, "select * from coexpressionOrfList_blaste4_02282023")

foo <- data.table::fread("/home/aar75/coexpression/20221110_rho/orfsNotInRNAseq.csv")

prev %>% filter(transcript %in% foo$transcript)
denovo <- allORFs %>% filter(classification == "denovo")
conserved <- allORFs %>% filter(classification == "conserved")

for (i in 1:16) {
    chr <- i
    df <- data.table::fread(glue::glue("/home/oma21/coexpression/data/interim/pairwise_orf_distances/pairwise_orf_distances_{chr}.csv"))
    df_sub <- df %>% filter(primary_orf %in% denovo$transcript & neighbor_orf %in% conserved$transcript)

    df_sub %>%
        group_by(primary_orf, orientation) %>%
        filter(distance == min(distance)) %>%
        ungroup() %>%
        select(primary_orf, neighbor_orf, distance, orientation) %>%
        arrange(primary_orf, distance) %>%
        write_csv(glue::glue("/home/oma21/coexpression/data/interim/pairwise_orf_distances/pairwise_orf_distances_{chr}_denovo_conserved_closest_blaste4.csv"))

    rm(df)
    rm(df_sub)
    gc()
}


pairwise_closest_denovo_conserved <- data.table::rbindlist(lapply(1:16, function(i) {
    data.table::fread(glue::glue("/home/oma21/coexpression/data/interim/pairwise_orf_distances/pairwise_orf_distances_{i}_denovo_conserved_closest_blaste4.csv"))
})) %>%
    left_join(
        allORFs[, c("transcript", "coor1", "strand")],
        by = c("primary_orf" = "transcript")
    ) %>%
    left_join(
        allORFs[, c("transcript", "coor1")],
        by = c("neighbor_orf" = "transcript")
    ) %>%
    mutate(
        distance = ifelse(orientation == "sense-overlap", 1, distance),
        orientation = case_when(
            orientation == "sense-overlap" & strand == "+" & coor1.x < coor1.y ~ "upsame",
            orientation == "sense-overlap" & strand == "+" & coor1.x > coor1.y ~ "downsame",
            orientation == "sense-overlap" & strand == "-" & coor1.x < coor1.y ~ "downsame",
            orientation == "sense-overlap" & strand == "-" & coor1.x > coor1.y ~ "upsame",
            TRUE ~ orientation
        )
    ) %>%
    select(primary_orf, neighbor_orf, distance, orientation)

pairwise_closest_denovo_conserved %>%
    distinct() %>%
    write_csv("/home/oma21/coexpression/data/interim/pairwise_orf_distances/pairwise_orf_distances_denovo_conserved_closest_all_blaste4.csv")
DBI::dbWriteTable(con, "pairwise_orf_distances_denovo_conserved_closest_blaste4", pairwise_closest_denovo_conserved %>% distinct(), overwrite = T)
