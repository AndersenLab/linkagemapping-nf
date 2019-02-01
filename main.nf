#!/usr/bin/env nextflow

// params:
// [1] riail phenotype dataframe - required (trait is in condition.trait format)
// [2] cross object - optional. default = marker (because it can be loaded with linkagemapping)
// [3] GWER/FDR - optional. default = GWER
// [4] CI 1.5-LOD - optional. default = chromosomal
// [5] condtrt flag - do I make a new column for condition.trait?

riails = params.in.replace(".tsv","")
params.cross = 'marker'
params.thresh = 'GWER'
params.ci = "chromosomal"
params.out = riails + "-mapping"
params.nperm = 1000
params.set = 2

// // perform principal component analysis on 24 main traits (exclude iqr, f., fluorescence)
// process perform_pca {
// 	publishDir "analysis-${params.out}/", mode: 'copy'
// 	input:
// 		file 'infile' from Channel.fromPath(params.in)

// 	output:
// 		file 'allRIAILsPCregressed.tsv' into riail_pheno

// 	"""
// 	#!/usr/bin/env Rscript --vanilla
// 	library(tidyverse)
// 	library(linkagemapping)

// 	# read in phenotypes and cross object
// 	# assign("cross", get(linkagemapping::load_cross_obj("${params.cross}")))

// 	if("${params.cross}" == "marker") {
// 		data("N2xCB4856cross")
// 		cross <- N2xCB4856cross
// 		markers <- "N2xCB4856"
// 	} else {
// 		load("${params.cross}")
// 		cross <- get(grep("cross", ls(), value = T))
// 		markers <- "full"
// 	}

// 	allRIAILregressed <- readr::read_tsv("${infile}")

// 	# keep only the strains in set2
// 	phenocross <- linkagemapping::mergepheno(cross, allRIAILregressed, set = 2)
// 	drugregressed <- phenocross\$pheno %>%
// 	    dplyr::filter(set == 2)

// 	allRIAILregressed <- allRIAILregressed %>%
// 	    dplyr::filter(strain %in% drugregressed\$strain) %>%
//     	dplyr::select(strain, condition, trait, phenotype)

// 	#for each drug, find the PCs and summarize them

// 	allPCs <- NULL
// 	for(d in unique(allRIAILregressed\$condition)){
	    
// 	    pc_traits <- dplyr::filter(allRIAILregressed, condition == d)%>%
// 	        dplyr::mutate(drugtrait = paste0(condition, ".", trait)) %>%
// 	        dplyr::select(drugtrait, strain, phenotype)%>%
// 	        dplyr::ungroup()%>%
// 	        dplyr::filter(!strain %in% c("N2","CB4856"), !grepl("red|green|yellow|f.ad|f.L1|f.L2L3|f.L4|iqr", drugtrait)) %>%
// 	        unique() %>%
// 	        tidyr::spread(drugtrait, phenotype) %>%
// 	        na.omit()
	    
// 	    row.names(pc_traits) <- pc_traits\$strain
	    
// 	    pc_traits <- pc_traits %>%
// 	        dplyr::select(-strain)
	    
// 	    scales_pc_traits <- as.data.frame(scale(pc_traits))
	    
// 	    pca_obj <- princomp(scales_pc_traits)
	    
// 	    # figure out how many PCs to keep that will explain > 90% of the variance
// 	    # pull the total variance explained with each PC and call it "drug.cumsum"
// 	    cumsums <- t(as.matrix(cumsum(pca_obj\$sdev^2/sum(pca_obj\$sdev^2))))
// 	    rownames(cumsums) <- paste(d, "cumsum", sep = ".")
// 	    cumsums <- as.data.frame(cumsums) %>%
// 	        tidyr::gather(comp, var) %>%
// 	        dplyr::filter(var > 0.9) 
	    
// 	    # the first component in the dataframe is the first one that goes over 90%, so we want to keep it and everything else below it
// 	    keep <- as.numeric(stringr::str_split_fixed(cumsums\$comp[1], "Comp.", 2)[,2])
	    
// 	    # keep only those PCs phenotypes for all RIAIL strains
// 	    colnames(pca_obj\$scores) <- paste(d, colnames(pca_obj\$scores), sep = "_")
	    
// 	    pcaoutput <- data.frame(pca_obj\$scores[,1:keep]) %>%
// 	        dplyr::mutate(strain = rownames(.)) %>%
// 	        tidyr::gather(trait, phenotype, -strain)
	    
// 	    allPCs <- rbind(allPCs, pcaoutput)
// 	}

// 	# tidy
// 	allPCs\$trait <- gsub("_Comp.", "_PC", allPCs\$trait)
// 	allRIAILsPC <- allPCs %>%
// 	    tidyr::separate(trait, c("condition", "trait"), sep = "_")

// 	# combine allRIAILsregressed and allRIAILsPC
// 	combined_phenotype <- rbind(allRIAILregressed, allRIAILsPC) %>%
// 		dplyr::mutate(condtrt = paste0(condition, ".", trait))

// 	readr::write_tsv(combined_phenotype, "allRIAILsPCregressed.tsv")

// 	"""
// }

// split phenotypes into separate files
process split_pheno {
	input:
		file 'infile' from Channel.fromPath(params.in)

	output:
        file '*.tsv' into phenotypes

    """
    #!/usr/bin/env Rscript --vanilla
    library(dplyr)
    library(readr)

    df <- readr::read_tsv("${infile}")

    # if theshold is FDR, do not split by trait
    if("${params.thresh}" == "FDR") {
    	readr::write_tsv(df, paste0("${riails}", "-phenotype.tsv"))
    } else if("${params.thresh}" == "GWER") {
    	# split phenotype by trait
	    for(i in unique(df["condtrt"])[[1]]) {
	    	traitdf <- df %>%
	        	dplyr::filter(condtrt == i)
	    	readr::write_tsv(traitdf, path = paste0(i, "-phenotype.tsv"))
		}
    } else {
    	stop("Error: Please select 'GWER' or 'FDR' threshold")
    }

    """
}

phenotypes_split = phenotypes.flatten()


// Perform the mapping and annotate
process mapping {
	cpus 4
	tag { input_tsv }
	publishDir "analysis-${params.out}/mappings", mode: 'copy'

	input:
		file input_tsv from phenotypes_split

	output:
		file("*.unannotated.tsv") into unannotated_maps
		file("*.annotated.tsv") into annotated_maps

	"""

	#!/usr/bin/env Rscript --vanilla
	library(linkagemapping)
	library(dplyr)
	library(readr)

	# load phenotype data
	df <- readr::read_tsv("${input_tsv}")

	# get trait name
	threshold <- "${params.thresh}"
	if(threshold == "FDR") {
		phenotype_name <- "${riails}"
	} else {
		phenotype_name <- unique(df["condtrt"])[[1]]
	}
	

	# load cross object
	# assign('cross', get(linkagemapping::load_cross_obj("${params.cross}")))

	if("${params.cross}" == "marker") {
		data("N2xCB4856cross")
		cross <- N2xCB4856cross
		markers <- "N2xCB4856"
	} else {
		load("${params.cross}")
		cross <- get(grep("cross", ls(), value = T))
		markers <- "full"
	}

	# make trait cross object
	drugcross <- linkagemapping::mergepheno(cross, df, set = "${params.set}")

	# choose the markerset
	#if(deparse(substitute(${params.cross})) == "N2xCB4856cross") {
	#	markers <- "N2xCB4856"
	#} else if(deparse(substitute(${params.cross})) == "N2xCB4856cross_full2") {
	#	markers <- 'full'
	#}

	# perform the mapping
	map <- linkagemapping::fsearch(drugcross, permutations = $params.nperm, thresh = threshold, markerset = markers)

	# save unannotated map
	readr::write_tsv(map, paste0(phenotype_name, "-", threshold, ".unannotated.tsv"))

	# annotate map
	cilod <- "${params.ci}"

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

	"""
}

// combine the annotated maps for each drug and output as r dataframe
process concatenate_drug {
	publishDir "analysis-${params.out}/", mode: 'copy'

	input:
		val(input_mapping) from annotated_maps.toSortedList()

	output:
		file("*.Rda")

	"""
	#!/usr/bin/env Rscript --vanilla
	library(tidyverse)

	# get the drugs
	# pheno <- readr::read_tsv(${params.in})

	drugs <- c('carmustine', 'chlorothanil', 'daunorubicin', 'docetaxel', 'etoposide', 'fluoxetine.125', 'fluoxetine.250', 'irinotecan', 
           'methotrexate.625', 'methotrexate.3125', 'thiabendazole.625', 'thiabendazole.125', 'tunicamycin', 'abamectin', 'albendazole', 
           'amsacrine', 'bortezomib', 'chlorpyrifos', 'dactinomycin', 'fenbendazole.15', 'fenbendazole.30', 'mebendazole', 'topotecan', 
           'mianserin', 'monepantel', 'arsenicdibasic', 'arsenictrioxide', 'bleomycin', 'cadmium', 'copper', 'deiquat', 'FUdR', 'mechlorethamine', 
           'nickel', 'paraquat', 'puromycin', 'silver', 'vincristine', 'zinc', 'cisplatin.250', 'cisplatin.500', 'lysate.175', 'OP50', 'DA837', 
           'JUb68', 'HT115')

    threshold <- "${params.thresh}"
    cilod <- "${params.ci}"

    for(num in 1:length(drugs)) {
    	# assign drug
    	d <- drugs[num]

    	# load all data associated with that drug
    	files <- grep("chrom", grep(d, list.files(), value = T), value = T)
    	annotatedmap <- NULL
    	for(f in files) {
    		traitdf <- readr::read_tsv(f)
    		traitdf$ci_l_pos <- as.numeric(traitdf$ci_l_pos)
			traitdf$ci_r_pos <- as.numeric(traitdf$ci_r_pos)
			traitdf$pos <- as.numeric(traitdf$pos)
			traitdf$lod <- as.numeric(traitdf$lod)
    		annotatedmap <- rbind(annotatedmap, traitdf)
    	}

    	save(annotatedmap, file = paste0("/projects/b1059/projects/Katie/20180622_alldrugGWERmap/analysis-20180827_allAnnotatedLods/", d, "-", threshold, ".", cilod, ".annotated.Rda"))
    }


    for(num in 1:length(drugs)) {
    	# assign drug
    	d <- drugs[num]

    	# load all data associated with that drug
    	load(grep(d, list.files(), value = T))
		annotatedmap$ci_l_pos <- as.numeric(annotatedmap$ci_l_pos)
		annotatedmap$ci_r_pos <- as.numeric(annotatedmap$ci_r_pos)
		annotatedmap$pos <- as.numeric(annotatedmap$pos)
		annotatedmap$lod <- as.numeric(annotatedmap$lod)
		annotatedmap$var_exp <- as.numeric(annotatedmap$var_exp)
		annotatedmap$eff_size <- as.numeric(annotatedmap$eff_size)


    	save(annotatedmap, file = paste0("/projects/b1059/projects/Katie/20180622_alldrugGWERmap/analysis-20180827_allAnnotatedLods/", d, "-", threshold, ".", cilod, ".annotated.Rda"))
    }

	"""
}


// combine the annotated mappings
// process concatenate_mappings {

//     publishDir "analysis-${params.out}/", mode: 'copy'
    
//     input:
//         val(input_mapping) from annotated_maps.toSortedList()

//     output:
//         file("${params.out}.tsv") into output

//     """
//     # use this to only print the header of the first line
// 	awk 'FNR>1 || NR==1' ${input_mapping.join(" ")} > ${params.out}.tsv
//     """

// }

// convert to Rda format
// process convertR {

// 	publishDir "analysis-${params.out}/", mode: 'copy'

// 	input:
// 		file("mappingfile") from output

// 	output:
// 		file("${params.out}.Rda")

// 	"""
// 	#!/usr/bin/env Rscript --vanilla
// 	library(readr)

// 	# load file
// 	df <- readr::read_tsv("$mappingfile")

// 	# convert to numerics
// 	df["var_exp"] <- as.numeric(df["var_exp"][[1]])
// 	df["eff_size"] <- as.numeric(df["eff_size"][[1]])
// 	df["ci_l_pos"] <- as.numeric(df["ci_l_pos"][[1]])
// 	df["ci_r_pos"] <- as.numeric(df["ci_r_pos"][[1]])

// 	# save dataframe as Rda
// 	${params.out} <- df
// 	save(${params.out}, file = paste0("${params.out}", ".Rda"))

// 	"""

// }