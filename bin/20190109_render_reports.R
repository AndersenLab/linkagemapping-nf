library(rmarkdown)
library(linkagemapping)
library(dplyr)

# Generate HTML report for each of the drugs mapped in linkage mapping
setwd("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/reports")

# all 46 drugs
drugs <- c('carmustine', 'chlorothanil', 'daunorubicin', 'docetaxel', 'etoposide', 'fluoxetine.125', 'fluoxetine.250', 'irinotecan', 
           'methotrexate.625', 'methotrexate.3125', 'thiabendazole.625', 'thiabendazole.125', 'tunicamycin', 'abamectin', 'albendazole', 
           'amsacrine', 'bortezomib', 'chlorpyrifos', 'dactinomycin', 'fenbendazole.15', 'fenbendazole.30', 'mebendazole', 'topotecan', 
           'mianserin', 'monepantel', 'arsenicdibasic', 'arsenictrioxide', 'bleomycin', 'cadmium', 'copper', 'deiquat', 'FUdR', 'mechlorethamine', 
           'nickel', 'paraquat', 'puromycin', 'silver', 'vincristine', 'zinc', 'cisplatin.250', 'cisplatin.500', 'lysate.175', 'OP50', 'DA837', 
           'JUb68', 'HT115')

# load phenotype data and genotype data
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsPCregressed.Rda")
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/RIAIL_data/allRIAILsregressed.Rda")
load_cross_obj("N2xCB4856cross_full")

for(num in 1:length(drugs)) {
    # assign drug
    d <- drugs[num]
    
    # load data
    load(paste0("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/data/", d, "-GWER.chromosomal.annotated.Rda"))
    
    # change concentrations to drug-conc
    if(d == "fluoxetine.125") annotatedmap$trait <- gsub("fluoxetine.125", "fluoxetine-125", annotatedmap$trait) else
    if(d == "fluoxetine.250") annotatedmap$trait <- gsub("fluoxetine.250", "fluoxetine-250", annotatedmap$trait) else
    if(d == "methotrexate.625") annotatedmap$trait <- gsub("methotrexate.625", "methotrexate-625", annotatedmap$trait) else
    if(d == "methotrexate.3125") annotatedmap$trait <- gsub("methotrexate.3125", "methotrexate-3125", annotatedmap$trait) else
    if(d == "thiabendazole.625") annotatedmap$trait <- gsub("thiabendazole.625", "thiabendazole-625", annotatedmap$trait) else
    if(d == "fenbendazole.15") annotatedmap$trait <- gsub("fenbendazole.15", "fenbendazole-15", annotatedmap$trait) else
    if(d == "fenbendazole.30") annotatedmap$trait <- gsub("fenbendazole.30", "fenbendazole-30", annotatedmap$trait) else
    if(d == "cisplatin.250") annotatedmap$trait <- gsub("cisplatin.250", "cisplatin-250", annotatedmap$trait) else
    if(d == "cisplatin.500") annotatedmap$trait <- gsub("cisplatin.500", "cisplatin-500", annotatedmap$trait) else
    if(d == "lysate.175") annotatedmap$trait <- gsub("lysate.175", "lysate-175", annotatedmap$trait) else
    if(d == "thiabendazole.125") annotatedmap$trait <- gsub("thiabendazole.125", "thiabendazole-125", annotatedmap$trait)
    
    d <- gsub("\\.", "-", d)
    
    # phenotype data
    pheno <- allRIAILsPCregressed %>%
        dplyr::filter(condition == d)
    
    # make cross object
    drugcross <- linkagemapping::mergepheno(N2xCB4856cross_full2, pheno, set = 2)
    
    # run markdown
    rmarkdown::render("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/scripts/20190109_fullgenomemapping_reports.Rmd", 
                      output_file = paste0("~/Dropbox/AndersenLab/RCode/Linkage mapping/RIAILsMappings/20180829/reports/Report_", d, ".html"))
}
