#!/usr/bin/env Rscript --vanilla
.libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths()))
library(tidyverse)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - mapping file

# load file
df <- readr::read_tsv(args[1])

# convert to numerics
df["var_exp"] <- as.numeric(df["var_exp"][[1]])
df["eff_size"] <- as.numeric(df["eff_size"][[1]])
df["ci_l_pos"] <- as.numeric(df["ci_l_pos"][[1]])
df["ci_r_pos"] <- as.numeric(df["ci_r_pos"][[1]])

# save dataframe as Rda
annotatedmap <- df
save(annotatedmap, file = paste0("annotatedmap.Rda"))