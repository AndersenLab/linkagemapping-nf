#!/usr/bin/env Rscript --vanilla
.libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths()))
library(tidyverse)
library(qtl)
library(linkagemapping)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - cross object
# 2 - permutation

# load cross object
load(args[1])

# run scan2 with 5 permutations
scan2thousand <- qtl::scantwo(drugcross, n.perm=1, pheno.col=3, method="mr")

# make output into dataframe
df <- data.frame(full = scan2thousand$full[[1]], fv1 = scan2thousand$fv1[[1]], int = scan2thousand$int[[1]], 
					add = scan2thousand$add[[1]], av1 = scan2thousand$av1[[1]], one = scan2thousand$one[[1]],
					trait = names(drugcross$pheno)[3])

# save dataframe
readr::write_tsv(df, paste0(names(drugcross$pheno)[3], "_scan2thousand_", args[2], ".tsv"))