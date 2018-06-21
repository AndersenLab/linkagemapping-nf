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
params.out = riails + "mapping"
params.nperm = 1000

// for use on personal computer
// params.out = "/Users/katieevans/Dropbox/AndersenLab/QTLpaper/scripts/GWER_mapping/"
// params.in = "/Users/katieevans/Dropbox/AndersenLab/QTLpaper/scripts/GWER_mapping/"
// riails = "/Users/katieevans/Dropbox/AndersenLab/RCode/Linkage\\ mapping/RIAILsMappings/allRIAILsregressed.Rda"


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

    # split phenotype by trait
    for(i in unique(df["trait"])[[1]]) {
    	traitdf <- df %>%
        	dplyr::filter(trait == i)
    	readr::write_tsv(traitdf, path = paste0(i, "-phenotype.tsv"))
	}

    """
}

phenotypes_split = phenotypes.flatten()


// Perform the mapping and annotate
process mapping {
	cpus 4
	tag { input_tsv }
	publishDir "analysis-${riails}/mappings", mode: 'copy'

	input:
		file input_tsv from phenotypes_split

	output:
		file("*.annotated.tsv") into annotated_maps

	"""

	#!/usr/bin/env Rscript --vanilla
	library(linkagemapping)
	library(dplyr)
	library(readr)

	# load phenotype data
	df <- readr::read_tsv("${input_tsv}")

	# get trait name
	phenotype_name <- unique(df["trait"])[[1]]

	# load cross object
	if("${params.cross}" == "marker") {
		data("N2xCB4856cross")
		cross <- N2xCB4856cross
	} else {
		load("${params.cross}")
		cross <- get(grep("cross", ls(), value = T))
	}

	# make trait cross object
	drugcross <- linkagemapping::mergepheno(cross, df, set = 2)

	# perform the mapping
	threshold <- "${params.thresh}"
	map <- linkagemapping::fsearch(drugcross, permutations = $params.nperm, thresh = threshold, markerset="N2xCB4856")

	# annotate map
	cilod <- "${params.ci}"
	annotatedmap <- linkagemapping::annotate_lods(map, drugcross, cutoff = cilod)

	# save annotated map
	readr::write_tsv(annotatedmap, paste0(phenotype_name, "-", threshold, ".", cilod, ".annotated.tsv"))

	"""
}


// combine the annotated mappings
process concatenate_mappings {

    publishDir "analysis-${riails}/", mode: 'copy'
    
    input:
        val(input_mapping) from annotated_maps.toSortedList()

    output:
        file("${params.out}.tsv") into output

    """
	cat ${input_mapping.join(" ")} > ${params.out}.tsv
    """

}
