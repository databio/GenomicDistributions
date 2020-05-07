# lib
library(data.table)
library(testthat)
library(GenomicDistributions)

# data
query = vistaEnhancers
querySftd = GenomicRanges::shift(query, 100000)
queryList = GRangesList(q1=query, q2=querySftd)

# tests
context("general")
test_that("calcPartitionsRef works", {
    lapply(queryList, function(x) expect_visible(calcPartitionsRef(x, "hg19")))
})

context("result")
test_that("calcPartitionsRef returns a result of a proper class", {
    expect_true(is(calcPartitionsRef(query, "hg19"), "data.frame"))
})

test_that("calcPartitionsRef returns a result of a proper length", {
    expect_length(calcPartitionsRef(query, "hg19"), 2)
    expect_equal(NROW(calcPartitionsRef(query, "hg19")), 5)
})

test_that("calcPartitionsRef returns different results for different queries", {
    expect_false(all(calcPartitionsRef(query, "hg19")$Freq == 
                     calcPartitionsRef(querySftd, "hg19")$Freq))
})