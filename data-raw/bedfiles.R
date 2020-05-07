fileNameList = c("vistaEnhancers.bed.gz", "setB_100.bed.gz")

for (fileName in fileNameList) {
    storedObjectName = strsplit(fileName, "\\.")[[1]][1]
    x = rtracklayer::import(system.file("extdata", fileName, package = "GenomicDistributions"))
    assign(storedObjectName, x)
    do.call("use_data", list(as.name(storedObjectName), overwrite = TRUE))
    rm(feats, storedObjectName)
}