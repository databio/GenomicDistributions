# Unit tests
library(GenomicDistributions)

context("Testthat context...")

#############################################################################
# Test data should be with toy examples you can work out by hand
# that way you can calculate by hand and compare to the output of the function

# toy data for testing functions
# if altered, tests relying on these objects will be disrupted
start1 = c(seq(from=1, to = 2001, by = 1000), 800)
start2 = c(seq(from=126, to = 2126, by = 1000), 100, 2500)
chrString1 = c(rep("chr1", 3), "chr2")
chrString2 = c(chrString1, "chr3")

origCoordDT1 <- data.table(chr=chrString1,
                           start = start1,
                           end = start1 + 250)
origCoordDT2 = data.table(chr=chrString2,
                          start=start2,
                          end=start2+150)
coordDT1 = copy(origCoordDT1)
coordDT2 = copy(origCoordDT2)

testGR1 = dtToGr(coordDT1)
testGR2 = dtToGr(coordDT2)
###############################################################################


# "featureDistanceDistribution" function is now named "calcFeatureDist"
# reset test data in case it was changed by another unit test section
coordDT1 = copy(origCoordDT1)
coordDT2 = copy(origCoordDT2)
testGR1 = dtToGr(coordDT1)
testGR2 = dtToGr(coordDT2)
test_that("featureDistribution",  {
    
    ############# old
    queryFile = system.file("extdata", "setB_100.bed.gz", package="GenomicDistributions")
    query = rtracklayer::import(queryFile)
    
    featureExample = GenomicRanges::shift(query, round(rnorm(length(query), 0,1000)))
    fdd = featureDistanceDistribution(query, featureExample)
    featureFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
    feats = rtracklayer::import(featureFile)
    
    #' featureDistance = featureDistanceDistribution(query, feats)
    #' expect_equal(sum(is.na(featureDistance)), -3)
    #' expect_equal(sum(featureDistance, na.rm=TRUE), 743969)
    ############# old
    
    coordDT1$end[1] = 100
    coordDT1$start[2] = 200
    coordDT1$end[2] = 400
    testGR1 = dtToGr(coordDT1)
    # DTNearest
    # @param DT1 data.table Has start and end column
    # @param DT2 
    # @return numeric vector. Distance from region set to closest other region set.
    # Distance from the midpointof each region to the midpoint.
    nearestVec = DTNearest(coordDT1, coordDT2)
    nearestVec
    expect_equal(nearestVec, c(150, -99, 75, -750))
    nearestVec2 = DTNearest(coordDT2, coordDT1)
    nearestVec2 # actual: c( 99, -901, -1276, 750, -1650)
    coordDT1
    coordDT2 #3 is matching to chr2 from coordDT1
    expect_equal(nearestVec2, c( 99, -901, -75, 750, NA))
    
    featureDistance = calcFeatureDist(testGR1, testGR2)
    featureDistance
    expect_equal(featureDistance, c(150, -99, 75, -750))
    featureDistance2 = calcFeatureDist(testGR2, testGR1)
    featureDistance2
    expect_equal(featureDistance2, c( 99, -901, -75, 750, NA))
    
    coordDT1$chr = "chr2"
    testGR1 = dtToGr(coordDT1)
    featureDistance = calcFeatureDist(testGR1, testGR2)
    featureDistance
    featureDistance2 = calcFeatureDist(testGR2, testGR1)
    featureDistance2
    coordDT1 = rbind(coordDT1, data.frame(chr="chr1", 
                                          start=1000000, 
                                          end=1000006, 
                                          mid=1000003))
    testGR1 = dtToGr(coordDT1)
    featureDistance = calcFeatureDist(testGR1, testGR2)
    featureDistance
    coordDT1
    coordDT2
    #expect_equal(featureDistance, )
    
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


# "genomicPartitions" function changed to "calcPartitionsRef"

test_that("Partitions", {
    
    ################### old
    queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
    query = rtracklayer::import(queryFile)
    gp = genomicPartitions(query, "hg38")
    gp = genomicPartitions(query, "hg19")
    gp = genomicPartitions(query, "mm10")
    gp = genomicPartitions(query, "mm9")
    plotPartitions(gp)
    ################### old
    
    # test calcPartitions()
    calcPartitions()
})