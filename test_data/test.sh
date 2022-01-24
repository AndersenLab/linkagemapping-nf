# test runs for pipeline

# first, debug mapping with low permutations, no scan2
nextflow run ../main.nf --debug

# next, debug mapping with low permutations, with scan2
nextflow run ../main.nf --debug --scan true

# next, debug scan2 with no mapping
nextflow run ../main.nf --debug --scan true --map false

# next, run mapping debug with full cross object
nextflow run ../main.nf --debug --cross N2xCB4856cross_full --scan true

# next, run a full mapping with scan2
nextflow run ../main.nf --in ../test_data/testdf.tsv --scan true --cross N2xCB4856cross_full