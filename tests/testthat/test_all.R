# Unit tests
library(GenomicDistributions)

context("Testthat context...")

test_that("featureDistribution",  {

	queryFile = system.file("extdata", "setB_100.bed.gz", package="GenomicDistributions")
	query = rtracklayer::import(queryFile)

	featureExample = GenomicRanges::shift(query, round(rnorm(length(query), 0,1000)))
	fdd = featureDistanceDistribution(query, featureExample)
	featureFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
	feats = rtracklayer::import(featureFile)

	#' featureDistance = featureDistanceDistribution(query, feats)
	#' expect_equal(sum(is.na(featureDistance)), -3)
	#' expect_equal(sum(featureDistance, na.rm=TRUE), 743969)
})

 #' queryDT = GenomicDistributions:::grToDt(query)
 #'        featureDT = GenomicDistributions:::grToDt(features)
 #'        queryDTs = GenomicDistributions:::splitDataTable(queryDT, "chr")
 #'        featureDTs = GenomicDistributions:::splitDataTable(featureDT, "chr")
 #'        as.vector(unlist(mapply(queryDTs, featureDTs[names(queryDTs)], FUN=DTNearest)))



test_that("Genome aggregate", {
	queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
	query = rtracklayer::import(queryFile)
	# First, calculate the distribution:
	x = aggregateOverGenomeBins(query, "hg19")
	# Then, plot the result:
	# plotGenomeAggregate(x)
})



test_that("Partitions", {
	queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
	query = rtracklayer::import(queryFile)
	gp = genomicPartitions(query, "hg38")
	gp = genomicPartitions(query, "hg19")
	gp = genomicPartitions(query, "mm10")
	gp = genomicPartitions(query, "mm9")
})