# Unit tests
library(GenomicDistributions)

context("Testthat context...")

test_that("featureDistribution",  {

	queryFile = system.file("extdata", "setB_100.bed.gz", package="GenomicDistributions")
	query = rtracklayer::import(queryFile)

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
