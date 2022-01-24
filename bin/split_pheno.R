#!/usr/bin/env Rscript --vanilla
# .libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths()))
library(tidyverse)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - traitfile
# 2 - threshold

df <- readr::read_tsv(args[1])

# if theshold is FDR, do not split by trait
if(args[2] == "FDR") {
	readr::write_tsv(df, path = "phenotype.tsv")
} else if(args[2] == "GWER") {
	# split phenotype by trait
    for(i in unique(df["trait"])[[1]]) {
    	traitdf <- df %>%
        	dplyr::filter(trait == i)
    	readr::write_tsv(traitdf, path = paste0(i, "-phenotype.tsv"))
	}
} else {
	stop("Error: Please select 'GWER' or 'FDR' threshold")
}