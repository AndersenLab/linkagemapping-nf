#!/usr/bin/env nextflow

// params:
// [1] riail phenotype dataframe - required (trait is in condition.trait format)
// [2] cross object - optional. default = marker (because it can be loaded with linkagemapping)
// [3] GWER/FDR - optional. default = GWER
// [4] CI 1.5-LOD - optional. default = chromosomal

riails = params.in.replace(".tsv","")
params.cross = 'marker'
params.thresh = 'GWER'
params.ci = "chromosomal"
params.out = riails + "-mapping"
params.nperm = 1000

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
    df <- readr::read_tsv("${infile}") %>%
    	dplyr::mutate(condtrt = paste0(condition, ".", trait))

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
	drugcross <- linkagemapping::mergepheno(cross, df, set = 2)

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


// combine the annotated mappings
process concatenate_mappings {

    publishDir "analysis-${params.out}/", mode: 'copy'
    
    input:
        val(input_mapping) from annotated_maps.toSortedList()

    output:
        file("${params.out}.tsv") into output

    """
    # use this to only print the header of the first line
	awk 'FNR>1 || NR==1' ${input_mapping.join(" ")} > ${params.out}.tsv
    """

}

// convert to Rda format
process convertR {

	publishDir "analysis-${params.out}/", mode: 'copy'

	input:
		file("mappingfile") from output

	output:
		file("${params.out}.Rda")

	"""
	#!/usr/bin/env Rscript --vanilla
	library(readr)

	# load file
	df <- readr::read_tsv("$mappingfile")

	# convert to numerics
	df["var_exp"] <- as.numeric(df["var_exp"][[1]])
	df["eff_size"] <- as.numeric(df["eff_size"][[1]])
	df["ci_l_pos"] <- as.numeric(df["ci_l_pos"][[1]])
	df["ci_r_pos"] <- as.numeric(df["ci_r_pos"][[1]])

	# save dataframe as Rda
	${params.out} <- df
	save(${params.out}, file = paste0("${params.out}", ".Rda"))

	"""

}