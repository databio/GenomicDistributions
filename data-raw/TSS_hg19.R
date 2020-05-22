library(usethis)
TSS_hg19 = GenomicDistributionsData::buildTSS("hg19")
usethis::use_data(TSS_hg19, overwrite=TRUE)
