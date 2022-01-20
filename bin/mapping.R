#!/usr/bin/env Rscript --vanilla
.libPaths(c("/projects/b1059/software/R_lib_3.6.0", .libPaths()))
library(tidyverse)
library(linkagemapping)
library(qtl)

# load arguments
args <- commandArgs(trailingOnly = TRUE)

# arguments
# 1 - phenotype file
# 2 - threshold
# 3 - cross object
# 4 - set
# 5 - perms
# 6 - ci
# 7 - scan2 boolean

# load phenotype data
df <- readr::read_tsv(args[1])

# get trait name
threshold <- args[2]
if(threshold == "FDR") {
	phenotype_name <- "phenotype"
} else {
	phenotype_name <- glue::glue("{df$condition[1]}.{df$trait[1]}")
}


# load cross object
assign('cross', get(linkagemapping::load_cross_obj(args[3])))

# choose the markerset
if(args[3] == "N2xCB4856cross_full") {
	markers <- NA
} else {
	# markers <- stringr::str_split_fixed("${params.cross}", "cross", 2)[,1]
	markers <- strsplit("N2xCB4856cross", "cross")[[1]]
}


# make trait cross object
drugcross <- linkagemapping::mergepheno(cross, df, set = args[4])
save(drugcross, file = glue::glue("{phenotype_name}-mapcross.Rda"))


# perform the mapping
if(args[8] %in% c("TRUE", "true", TRUE)) {
  map <- linkagemapping::fsearch(drugcross, permutations = args[5], thresh = threshold, markerset = markers)

  # annotate map
  cilod <- args[6]
  
  # check to make sure there is something to annotate
  peaks <- map %>%
    na.omit()
  if(nrow(peaks) > 0) {
    annotatedmap <- linkagemapping::annotate_lods(map, drugcross, cutoff = cilod)
  } else {
    annotatedmap <- map %>%
      dplyr::mutate(var_exp = NA, eff_size = NA, ci_l_marker = NA, ci_l_pos = NA, ci_r_marker = NA, ci_r_pos = NA)
  }
  
  # save annotated map
  readr::write_tsv(annotatedmap, paste0(phenotype_name, "-", threshold, ".", cilod, ".annotated.tsv"))
  
  # plot LOD plot
  lod <- linkagemapping::maxlodplot(annotatedmap)
  ggsave(lod, filename = paste0(phenotype_name, "-", threshold, ".", cilod, ".lod.png"), width = 8, height = 4)
  
  # plot PXG plot
  pxg <- linkagemapping::pxgplot(drugcross, annotatedmap)
  ggsave(pxg, filename = paste0(phenotype_name, "-", threshold, ".", cilod, ".pxg.png"))
  
  
}

if(args[7] %in% c("TRUE", "true", TRUE)) {
	# run scan2
	scan <- qtl::scantwo(drugcross, pheno.col=3, method="mr")

	# make output into dataframe
	save(scan, file = glue::glue("{phenotype_name}_scan2.Rda"))

	# plot scan2
	png(glue::glue("{phenotype_name}_scan2plot.png"))
	plot(scan)
	dev.off()
}
