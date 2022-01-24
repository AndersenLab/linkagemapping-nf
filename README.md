# linkagemapping-nf

Perform linkage mapping for one or more traits using nextflow.

### Requirements
The linkagemapping-nf pipeline requires installation of several R packages including [linkage mapping]("https://github.com/AndersenLab/linkagemapping"), `qtl` and `tidyverse`. These packages are all available on Quest with no additional work for the user.

This pipeline also requires Nextflow version 20.0+ and Singularity (or docker). On Quest, run with:

```
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env

module load singularity
```

### Usage
```
nextflow run AndersenLab/linkagemapping-nf --in=<input.tsv>
```

### Input
For the Andersen Lab, the input file is generally going to be the output from the `easysorter` pipeline. Multiple traits can be represented, and you want the trait column to be in the format of condition_trait. In general, you need a tsv dataframe of the following format:

| condition | trait | strain | phenotype |
| --- | --- | --- | --- |
| docetaxel | mean.TOF | QX318 | -31.25 |
| docetaxel | mean.TOF | QX299 | 5.154 |
| ... | ... | ... | ... |
| zinc | norm.n | QX198 | 0.735 |
| zinc | norm.n |QX269 | 25.68 | 

#### Optional parameters
| param | default | optional | explanation |
| --- | --- | --- | --- |
| `--cross` | 'N2xCB4856cross' | 'file_to_cross' | default 'N2xCB4856cross' uses cross object generated from genotypes of 1400 markers. Other option to use WGS cross object (N2xCB4856cross_full) |
| `--thresh` | 'GWER' | 'FDR' | defines method for setting significance threshold in linkage mapping. Options are GWER (default) and FDR. GWER selects the threshold for significant QTL for each trait independently. FDR selects the threshold for significant QTL for all traits/drugs at the same time. |
| `--ci` | 'chromosomal' | 'proximal' | method for defining 1.5-LOD drop confidence intervals for QTL. Default is 'chromosomal' meaning the 1.5-LOD drop applies to the entire chromosome. The other option is 'proximal' which means the confidence interval will span 1.5-LOD units directly to the left/right of the peak marker. |
| `--nperm` | 1000 | n | Number of permutations to run for mapping and scan2. Minimum recommended is 100. |
| `--set` | 2 | 1 | N2xCB4856 RIAIL set to include for mapping. |
| `--scan` | FALSE | TRUE | Option to run scan2 and permutations for each trait mapped. |
| `--map` | TRUE | FALSE | Option to not run mapping (only use if scan = TRUE). |

#### Example with options
```
nextflow run AndersenLab/linkagemapping-nf --in='input.tsv' --ci='proximal' --cross='N2xCB4856cross_full' --scan TRUE
```
