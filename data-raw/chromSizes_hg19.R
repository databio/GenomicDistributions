library(usethis)
chromSizes_hg19 = GenomicDistributionsData::buildChromSizes("hg19")
usethis::use_data(chromSizes_hg19, overwrite=TRUE)
