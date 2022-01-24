#!/usr/bin/env Rscript --vanilla
# .libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths()))
library(tidyverse)
library(linkagemapping)
library(qtl)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - empty cross object
# 2 - scan2 object
# 3 - permutations


#insert cross data
crossobj <- get(linkagemapping::load_cross_obj(args[1]))

# load scan2
load(args[2])

# load perms
perms <- readr::read_tsv(args[3])

# summarize scantwo and add GWER thresholds defined by permutations
scan2_summary <- summary(scan) %>%
    dplyr::mutate(fv1_thresh = quantile(perms$fv1, probs = 0.95),
                  full_thresh = quantile(perms$full, probs = 0.95),
                  add_thresh = quantile(perms$add, probs = 0.95),
                  av1_thresh = quantile(perms$av1, probs = 0.95),
                  int_thresh = quantile(perms$int, probs = 0.95)) %>%
    dplyr::mutate(pos1f = as.character(pos1f),
                  pos2f = as.character(pos2f),
                  pos1a = as.character(pos1a),
                  pos2a = as.character(pos2a))

# how to get marker positions for both N2xCB cross objects
if(args[1] == "N2xCB4856cross_full") {
    # riail marker conversion
    mappos <- qtl::pull.map(crossobj, as.table = TRUE) %>%
        dplyr::mutate(marker = rownames(.),
                      cM = as.character(pos)) %>%
        dplyr::select(-pos, -chr) %>%
        dplyr::distinct(cM, .keep_all = T) 
    
    # convert genetic pos to genomic pos
    test <- scan2_summary %>%
        # pos1f
        dplyr::left_join(mappos, by = c("pos1f" = "cM")) %>%
        dplyr::mutate(pos1f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos2f
        dplyr::left_join(mappos, by = c("pos2f" = "cM")) %>%
        dplyr::mutate(pos2f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos1a
        dplyr::left_join(mappos, by = c("pos1a" = "cM")) %>%
        dplyr::mutate(pos1a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos2a
        dplyr::left_join(mappos, by = c("pos2a" = "cM")) %>%
        dplyr::mutate(pos2a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker)
} else if(args[1] == "N2xCB4856cross") {
    data("N2xCB4856markers")
    
    # riail marker conversion
    mappos <- qtl::pull.map(crossobj, as.table = TRUE) %>%
        dplyr::mutate(marker = rownames(.),
                      cM = as.character(pos)) %>%
        dplyr::select(-pos, -chr) %>%
        dplyr::distinct(cM, .keep_all = T) %>%
        dplyr::left_join(N2xCB4856markers) %>%
        dplyr::mutate(marker = paste0(chr.roman, "_", position)) %>%
        dplyr::select(marker, cM)
    
    # convert genetic pos to genomic pos
    test <- scan2_summary %>%
        # pos1f
        dplyr::left_join(mappos, by = c("pos1f" = "cM")) %>%
        dplyr::mutate(pos1f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos2f
        dplyr::left_join(mappos, by = c("pos2f" = "cM")) %>%
        dplyr::mutate(pos2f = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos1a
        dplyr::left_join(mappos, by = c("pos1a" = "cM")) %>%
        dplyr::mutate(pos1a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker) %>%
        # pos2a
        dplyr::left_join(mappos, by = c("pos2a" = "cM")) %>%
        dplyr::mutate(pos2a = as.numeric(stringr::str_split_fixed(marker, "_", 2)[,2])) %>%
        dplyr::select(-marker)
} else {
    # return scan2summary in cM instead of bp
    test <- scan2_summary
}

test <- test %>%
        dplyr::mutate(trait = perms$trait[1])

readr::write_tsv(test, glue::glue("{test$trait[1]}_scan2summary.tsv"))
