# library(GenomicDistributions)
library(testthat)
library(data.table)

# data
query = vistaEnhancers
querySftd = GenomicRanges::shift(query, 100000)
queryList = GRangesList(q1=query, q2=querySftd)

context("general")
test_that("binRegion works with binSize and binCount", {
    for(s in seq(1, 100, by=50)){
        for(e in seq(1000, 10000, by=5000)){
            expect_visible(binRegion(start=s, end=e, binSize=10))
            expect_visible(binRegion(start=s, end=e, binCount=10))
        }
    }
})

test_that("calcChromBinsRef works with list input", {
    expect_visible(calcChromBinsRef(queryList, "hg19"))
})

context("result")
test_that("binRegion returns result of correct length", {
    expect_equal(
        binRegion(start=1, end=100, binSize=10),
        binRegion(start=1, end=100, binCount=10),
    )
    expect_length(binRegion(start=1, end=100, binSize=10), 5)
    expect_equal(NROW(binRegion(start=1, end=100, binSize=10)), 10)
})

test_that("calcChromBinsRef returns a proper object type, length ad includes all the regions", {
    result = calcChromBinsRef(query, "hg19")
    expect_is(result, "data.table")
    expect_length(result, 6)
    expect_equal(sum(result$N), length(query))
})
