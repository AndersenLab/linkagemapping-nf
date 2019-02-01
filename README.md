# linkagemapping-nf

Perform linkage mapping for one or more traits using nextflow.

### Requirements
The linkagemapping-nf pipeline requires installation of several R packages including [linkage mapping]("https://github.com/AndersenLab/linkagemapping"), `dplyr`, and `readr`.

```
# Install linkagemapping package
devtools::install_github("AndersenLab/linkagemapping")
```

### Usage
```
nextflow run AndersenLab/linkagemapping-nf --in=<input.tsv>
```

### Input
For the Andersen Lab, the input file is generally going to be the output from the `easysorter` pipeline. Multiple traits can be represented, and you want the trait column to be in the format of condition_trait. In general, you need a tsv dataframe of the following format:

| trait | strain | phenotype |
| --- | --- | --- |
| docetaxel.mean.TOF | QX318 | -31.25 |
| docetaxel.mean.TOF | QX299 | 5.154 |
| ... | ... | ... |
| zinc.norm.n | QX198 | 0.735 |
| zinc.norm.n | QX269 | 25.68 | 

#### Optional parameters
| param | default | optional | explanation |
| --- | --- | --- | --- |
| `--cross` | 'marker' | 'file_to_cross' | default 'marker' uses cross object generated from genotypes of 1400 markers. Otherwise, include full path to whole genome cross object |
| `--thresh` | 'GWER' | 'FDR' | defines method for setting significance threshold in linkage mapping. Options are GWER (default) and FDR. GWER selects the threshold for significant QTL for each trait independently. FDR selects the threshold for significant QTL for all traits/drugs at the same time. |
| `--ci` | 'chromosomal' | 'proximal' | method for defining 1.5-LOD drop confidence intervals for QTL. Default is 'chromosomal' meaning the 1.5-LOD drop applies to the entire chromosome. The other option is 'proximal' which means the confidence interval will span 1.5-LOD units directly to the left/right of the peak marker. |
| `--out` | <name_of_input-mapping.tsv> | <output_file.tsv> | define name of output .tsv file | 

#### Example with options
```
nextflow run AndersenLab/linkagemapping-nf --in='input.tsv' --ci='proximal' --out='output.tsv'
```

Test
