# Generate chromSizes for selected genome assemblies
storedObjectName = "chromSizes_hg19"
message(paste("building", storedObjectName))
BSG = loadBSgenome(refAssembly)
chromSizesGenome = seqlengths(BSG)
assign(storedObjectName, chromSizesGenome)
do.call("use_data", list(as.name(storedObjectName), overwrite = TRUE))
rm(chromSizesGenome, storedObjectName)
