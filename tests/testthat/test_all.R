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

origCoordDT1 = data.table(chr=chrString1,
                           start = start1,
                           end = start1 + 250)
origCoordDT2 = data.table(chr=chrString2,
                          start=start2,
                          end=start2+150)
coordDT1 = copy(origCoordDT1)
coordDT2 = copy(origCoordDT2)

testGR1 = dtToGr(coordDT1)
testGR2 = dtToGr(coordDT2)
testGR3 = GenomicRanges::shift(testGR2, 1000)
testGR4 = GenomicRanges::shift(testGR2, 2500)
testGR5 = GenomicRanges::shift(testGR2, 4000)
###############################################################################

# test for calcOLCount
# reset test data in case it was changed by another unit test section
coordDT1 = copy(origCoordDT1)
coordDT2 = copy(origCoordDT2)
testGR1 = dtToGr(coordDT1)
testGR2 = dtToGr(coordDT2)
test_that("calcOLCount", {
    
    # uses midpoint coordinate of queryRegionDT
    testGRList = GRangesList(dtToGr(data.table(chr=c("chr1", "chr1"),
                                   start = c(1, 2001),
                                   end = c(2000, 4000))), 
                             dtToGr(data.table(chr=c("chr2", "chr2"),
                                               start = c(1, 2001),
                                               end = c(2000, 4000))),
                             dtToGr(data.table(chr=c("chr3", "chr3"),
                                               start = c(1, 2001),
                                               end = c(2000, 4000))))
    olCount1 = calcOLCount(queryRegionDT = coordDT2, regionsGRL = testGRList)
    expect_equal(olCount1$N, c(2, 1, 1, 1))
    expect_equal(olCount1$regionGroupID, c(1, 1, 2, 3))
    
    # only expect one overlap: chr2
    olCount2 = calcOLCount(coordDT2, dtToGr(data.table(chr=c("chr1", "chr1", "chr2"),
                                           start = c(1, 250, 170),
                                           end = c(150, 300, 180))))
    olCount2=as.data.frame(olCount2)
    expectedOut = data.frame(regionID=3, chr="chr2", start=170, end=180, withinGroupID=3, regionGroupID=1, N=1, stringsAsFactors = FALSE)
    expect_equal(olCount2, expectedOut)
})




# "featureDistanceDistribution" function is now named "calcFeatureDist"
# reset test data in case it was changed by another unit test section
# and select just one chromosome - since DTNearest is help function calculating
# distances within one chromosome
coordDT1 = copy(origCoordDT1)
coordDT2 = copy(origCoordDT2)
testGR1 = dtToGr(coordDT1)
testGR2 = dtToGr(coordDT2)
test_that("featureDistribution",  {
    
    ############# old
    # queryFile = system.file("extdata", "setB_100.bed.gz", package="GenomicDistributions")
    # query = rtracklayer::import(queryFile)
    # 
    # featureExample = GenomicRanges::shift(query, round(rnorm(length(query), 0,1000)))
    # fdd = featureDistanceDistribution(query, featureExample)
    # featureFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
    # feats = rtracklayer::import(featureFile)
    
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
    expect_equal(nearestVec, c(124, -99, 276, 75))
  

    # DTNearest ignores chromosome completely. By design.
    # DTNearest shouldn't be used with data from different chromosomes.
    # Suggested to split by chromosome when such case presents (e.g chrom1).
    DT1chrom1 = coordDT1[coordDT1$chr == "chr1"]
    DT2chrom1 = coordDT2[coordDT2$chr == "chr1"]
    nearestVec2C1 = DTNearest(DT2chrom1, DT1chrom1)
    expect_equal(nearestVec2C1, c(99, -901, -75))
    
    featureDistance = calcFeatureDist(testGR1, testGR2)
    featureDistance
    expect_equal(featureDistance, c(150, -99, 75, -750))
    featureDistance2 = calcFeatureDist(testGR2, testGR1)
    featureDistance2
    
    expect_equal(featureDistance2, c( 99, -901, -75, 750, NA))
    
    # coordDT1$chr = "chr2"
    # testGR1 = dtToGr(coordDT1)
    # featureDistance = calcFeatureDist(testGR1, testGR2)
    # featureDistance
    # featureDistance2 = calcFeatureDist(testGR2, testGR1)
    # featureDistance2

    
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
    #queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
    #query = rtracklayer::import(queryFile)
    #gp = genomicPartitions(query, "hg38")
    #gp = genomicPartitions(query, "hg19")
    #gp = genomicPartitions(query, "mm10")
    #gp = genomicPartitions(query, "mm9")
    #plotPartitions(gp)
    ################### old
    
    # test calcPartitions()
    # GenomePartitionList 
    promCore = GenomicRanges::reduce(trim(promoters(testGR2, upstream=100, downstream=0)))
    promProx = GenomicRanges::reduce(trim(promoters(testGR2, upstream=2000, downstream=0)))
    promoterProx = GenomicRanges::setdiff(promProx, promCore)
    
    # remove any possible overlaps between classes
    testGR5 = GenomicRanges::setdiff(testGR5, testGR4)
    testGR3 = GenomicRanges::setdiff(testGR3, testGR4)
    testGR3 = GenomicRanges::setdiff(testGR3, testGR5)
    
    nonThree = GenomicRanges::setdiff(testGR2, testGR4)
    nonThreeFive = GenomicRanges::setdiff(nonThree, testGR5)
    intronGR = GenomicRanges::setdiff(nonThreeFive, testGR3)
    
    partList = list(promoterCore=GenomicRanges::reduce(trim(promoters(testGR2, upstream=100, downstream=0))),
                    promoterProx=promoterProx, 
                    threeUTR=testGR4, 
                    fiveUTR=testGR5,
                    exon=testGR3,
                    intron=intronGR)
  
    gp = genomePartitionList(testGR2, testGR3, testGR4, testGR5)
    expect_equal(gp, partList)
    
    # calcPartitions
    partition = rep(0, length(testGR1))
    for (i in seq_along(partList)) {
      ols = countOverlaps(testGR1[partition==0], partList[[i]])
      partition[partition==0][ols > 0] = names(partList)[[i]]
    }
    partition[partition=="0"] = "intergenic"
    testPartitions = data.frame(table(partition))
    
    testPartitionNames = c("promoterCore", "promoterProx", "threeUTR", "fiveUTR",
                           "exon", "intron", "intergenic")
    if (!all(testPartitionNames %in% testPartitions$partition)){
      notIncluded = testPartitionNames[!(testPartitionNames %in% 
                                           testPartitions$partition)]
      addRows = data.frame(partition = notIncluded, 
                           Freq = rep(0, length(notIncluded)))
      testPartitions = rbind(testPartitions, addRows)
    }
    
    Partitions = calcPartitions(testGR1, partList)
    expect_equal(Partitions, testPartitions)
    
})

test_that("Neighbor distances", {
  
  testGRdt = grToDt(sort(testGR1))
  splitdt = splitDataTable(testGRdt, "chr")
  chromTest = splitdt[[1]]
  # Compare bp distance generated by neighbordt
  distancesExp = neighbordt(chromTest)
  # Calculated by hand c(749, 749)
  expect_equal(distancesExp, c(749, 749))
  
  # Compare  distances from calcNeighborDist
  distances = calcNeighborDist(testGR1)
  expect_equal(distances, c(749, 749))
  
})

test_that("Nearest Neighbor distances", {
  
  testGR2dt = grToDt(sort(testGR2))
  splitdt2 = splitDataTable(testGR2dt, "chr")
  chromTest2 = splitdt2[[1]]
  # Compare bp distance generated by neighbordt
  nearestDistancesExp = neighbordt(chromTest2)
  up = nearestDistancesExp[-length(dist)]
  down = nearestDistancesExp[-1]
  dt = data.table(i=up, j=down)
  pairmins = dt[, pmin(i, j)]
  nNeighbors = c(nearestDistancesExp[1], pairmins, 
                 nearestDistancesExp[length(dist)])
  
  # Calculated by hand c(849, 849, 849)
  expect_equal(nNeighbors, rep(849, 3))
  
  # Compare  distances from calcNeighborDist
  nearestNeighborsTest = calcNearestNeighbors(testGR2)
  expect_equal(nearestNeighborsTest, rep(849, 3))
  
})

