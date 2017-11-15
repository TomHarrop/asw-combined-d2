#!/usr/bin/env Rscript

library(data.table)

minimap_results_file <- "output/minimap/k_99/trim_decon/results.paf"

# read results file
minimap_results <- fread(minimap_results_file, fill = TRUE)
col_names <- c("query_name",
               "query_length",
               "query_start",
               "query_end",
               "strand",
               "target_name",
               "target_length",
               "target_start",
               "target_end",
               "no_residue_matches",
               "alignment_block_length",
               "mapping_quality")
setnames(minimap_results,
         names(minimap_results)[c(1:12)],
         col_names)

# parse for 1:1 matches
no_self_hits <- minimap_results[query_name != target_name]
no_self_hits[, query_block_cov := no_residue_matches / query_length]
no_self_hits[, target_block_cov := no_residue_matches / target_length]

hit_cutoff <- 0.75
putative_hits <- no_self_hits[target_block_cov > hit_cutoff & 
                 query_block_cov > hit_cutoff]

putative_hits[, sum(query_length)]

setorder(putative_hits, target_length)

