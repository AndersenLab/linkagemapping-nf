#!/usr/bin/env nextflow
// test
nextflow.preview.dsl=2

date = new Date().format( 'yyyyMMdd' )

// Debug parameters
if(params.debug) {
	println """

        *** Using debug mode ***

    """
	params.nperm = 10
	params.out = "debug-${date}"
	params.in = "${workflow.projectDir}/test_data/testdf.tsv"
} else {
	params.nperm = 1000
	params.out = "Analysis-${date}"
}

// rest of the parameters
riails = params.in.replace(".tsv","")
params.cross = 'N2xCB4856cross'
params.thresh = 'GWER'
params.ci = "chromosomal"
nperms = params.nperm
params.set = 2
params.scan = false
params.map = true


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

nextflow main.nf --debug
nextflow main.nf --in=/path/phenotype.tsv

    parameters                  description                           	Set/Default
    ==========                  ===========                           	========================
    --debug                     Use --debug to indicate debug mode      ${params.debug}
    --in                        Trait file                              ${params.in}
    --crossobj                  Cross object file name                  ${params.cross}
    --thresh                    Threshold (GWER or FDR)                 ${params.thresh}
    --ci                        Confidence interval method              ${params.ci}
    --nperm                     Number of permutations                  ${params.nperm}
    --set                       RIAIL set                               ${params.set}
    --scan                      To perform scan2 or not                 ${params.scan}
    --map                       To perform linkage mapping or not       ${params.map}

    username                                                            ${"whoami".execute().in.text}


---
"""
out
}

log.info(log_summary())

workflow {

	// linkage mapping
	Channel.fromPath(params.in) | split_pheno
	split_pheno.out.flatten() | mapping

	if(params.map == true) {
		mapping.out.annotated_map.collectFile(keepHeader: true, name: 'annotatedmap.tsv', storeDir: "Analysis-${date}") | convertR
	}

	if(params.scan == true) {
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
				.join(scanmap) | summarize_scan2 // join with scan2.Rda object and feed into summarize
	 }

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
	publishDir "${params.out}/scan2", mode: 'copy', pattern: "*scan2.Rda"
	publishDir "${params.out}/plots", mode: 'copy', pattern: "*.png"
	publishDir "${params.out}/plots", mode: 'copy', pattern: "*.pdf"
	publishDir "${params.out}/mappings", mode: 'copy', pattern: "*mapcross.Rda"

	input:
		file("input_tsv")

	output:
		path "*.annotated.tsv", emit: annotated_map, optional: true
		path "*scan2.Rda", emit: scan2_object, optional: true
		path "*mapcross.Rda", emit: crossobj
		path "*.png", optional: true
		path "*.pdf", optional: true


	"""
	Rscript --vanilla ${workflow.projectDir}/bin/mapping.R ${input_tsv} ${params.thresh} ${params.cross} ${params.set} ${params.nperm} ${params.ci} ${params.scan} ${params.map}

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

	"""
	Rscript --vanilla ${workflow.projectDir}/bin/scan_two_thousand.R ${mapcross} ${s}
	trait=`echo *scan2thousand*.tsv | cut -d '_' -f1`

	"""
}


// summarize scan2 and add thresholds for significance from permutations
process summarize_scan2 {

	publishDir "${params.out}/scan2", mode: "copy"

	cpus 4
	memory 16.GB

	input:
		tuple val("trait"), file("scantwothousand"), file("scan2")

	output:
		file("*scan2summary.tsv")


	"""
	Rscript --vanilla ${workflow.projectDir}/bin/scan2_summary.R ${params.cross} ${scan2} ${scantwothousand}

	"""


}

// Make HTML report and plots

// process html_report {

// 	input:
// 		tuple path("pheno"), path("map")

// 	output:
// 		file("*")
	
// 	"""


// 	"""

// }

