# library(GenomicDistributions)
library(testthat)
library(data.table)

context("general")
test_that("binRegion works with binSize and binCount", {
    for(s in seq(1,100, by=10)){
        for(e in seq(1000, 10000, by=100)){
            expect_visible(binRegion(start=s, end=e, binSize=10))
            expect_visible(binRegion(start=s, end=e, binCount=10))
        }
    }
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



