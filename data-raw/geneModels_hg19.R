library(usethis)
geneModels_hg19 = GenomicDistributionsData::buildGeneModels("hg19")
usethis::use_data(geneModels, overwrite=TRUE)
