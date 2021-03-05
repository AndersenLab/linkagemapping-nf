#!/usr/bin/env nextflow

nextflow.preview.dsl=2

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
date = new Date().format( 'yyyyMMdd' )
params.nperm = 10
nperms = params.nperm
params.set = 2
params.scan = "FALSE"
params.out = "Analysis-${date}"


def log_summary() {

    out =  '''

---------------------------------
---------------------------------
     * Linkagemapping-nf *
---------------------------------
---------------------------------
                                              
'''

out += """

To run the pipeline:

nextflow main.nf --in=/path/phenotype.tsv

    parameters                 	description                           	Set/Default
    ==========                 	===========                           	========================
    --crossobj					Cross object file name 					${params.cross}
    --thresh 					Threshold (GWER or FDR)  				${params.thresh}
    --ci 						Confidence interval method 				${params.ci}
    --nperm 					Number of permutations 					${params.nperm}
    --set 						RIAIL set   							${params.set}
    --scan 						To perform scan2 or not  				${params.scan}

---
"""
out
}

log.info(log_summary())

workflow {

	// linkage mapping
	Channel.fromPath(params.in) | split_pheno
	split_pheno.out.flatten() | mapping
	mapping.out.annotated_map.collectFile(keepHeader: true, name: 'annotatedmap.tsv', storeDir: "Analysis-${date}") | convertR

	// Scan2000
	seeds = Channel
			.from(1..100000)
			.randomSample(nperms)
			.combine(mapping.out.crossobj) | scan2000

	scanmap = mapping.out.scan2_object
			.map { file -> tuple(file.baseName.split('_')[0], file) } // grab the trait name from file and create tuple/set

	// summarize scan2
	scan2000.out
			.map { file -> tuple(file.baseName.split('_')[0], file) } // grab the trait name from file and create tuple/set
			.groupTuple() // group by trait
			.collectFile(keepHeader: true, sort: {it[0]}) {it[1]} // collect all the perms by trait and feed into one file
			.map { file -> tuple(file.baseName.split('_')[0], file) } // grab the trait name from file and create tuple/set
			.join(scanmap).view() | summarize_scan2 // join with scan2.Rda object and feed into summarize


}


// split phenotypes into separate files
process split_pheno {
	input:
		file('infile')

	output:
        file('*.tsv')

    """
    Rscript --vanilla ${workflow.projectDir}/bin/split_pheno.R ${infile} ${params.thresh}

    """
}



// Perform the mapping and annotate
process mapping {
	cpus 4
	memory '20 GB'

	publishDir "${params.out}/mappings", mode: 'copy', pattern: "*.tsv"
	publishDir "${params.out}/scan2", mode: 'copy', pattern: "*.Rda"
	publishDir "${params.out}/scan2", mode: 'copy', pattern: "*.png"

	input:
		file("input_tsv")

	output:
		path "*.annotated.tsv", emit: annotated_map
		path "*scan2.Rda", emit: scan2_object
		path "*mapcross.Rda", emit: crossobj
		file "*scan2plot.png"

	"""
	Rscript --vanilla ${workflow.projectDir}/bin/mapping.R ${input_tsv} ${params.thresh} ${params.cross} ${params.set} ${params.nperm} ${params.ci} ${params.scan}

	"""
}


// convert to Rda format
process convertR {

	publishDir "${params.out}", mode: 'copy'

	input:
		file("annotatedmap")

	output:
		file("annotatedmap.Rda")

	"""
	Rscript --vanilla ${workflow.projectDir}/bin/convert_to_rda.R ${annotatedmap}

	"""

}


// do the permutations
process scan2000 {
	cpus 4
	memory '20 GB'
	tag { s }

	input:
		tuple val("s"), file("mapcross")

	output:
		file("*scan2thousand*.tsv")

	//when:
	//	"${params.scan}" == TRUE

	"""
	Rscript --vanilla ${workflow.projectDir}/bin/scan_two_thousand.R ${mapcross} ${s}
	trait=`echo *scan2thousand*.tsv | cut -d '_' -f1`

	"""
}


// summarize scan2 and add thresholds for significance from permutations
process summarize_scan2 {

	publishDir "${params.out}/scan2", mode: "copy"

	input:
		tuple val("trait"), file("scantwothousand"), file("scan2")

	output:
		file("*scan2summary.tsv")

	//when:
	//	"${params.scan}" == TRUE


	"""
	Rscript --vanilla ${workflow.projectDir}/bin/scan2_summary.R ${params.cross} ${scan2} ${scantwothousand}

	"""


}



