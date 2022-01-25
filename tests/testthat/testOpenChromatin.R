# lib
library(data.table)
library(testthat)
library(GenomicDistributions)
# data
cellMatrix = exampleOpenSignalMatrix_hg19
query = vistaEnhancers
querySftd = GenomicRanges::shift(query, 100)
queryList = GRangesList(q1=query, q2=querySftd)

# tests
context("general")
test_that("calcSummarySignal works", {
    expect_visible(calcSummarySignal(query, cellMatrix))
    expect_visible(calcSummarySignal(querySftd, cellMatrix))
})

test_that("ccalcSummarySignal works with multiple queries", {
    expect_visible(calcSummarySignal(queryList, cellMatrix))
})

context("result")
test_that("calcSummarySignal returns a result of a proper class", {
    expect_true(is(calcSummarySignal(query, cellMatrix), "list"))
    expect_true(is(calcSummarySignal(query, cellMatrix)[[1]], "data.table"))
    expect_true(is(calcSummarySignal(query, cellMatrix)[[2]], "data.frame"))
})

test_that("calcSummarySignal returns different results for different queries", {
    expect_false(identical(calcSummarySignal(query, cellMatrix)[[1]], 
                           calcSummarySignal(querySftd, cellMatrix)[[1]]))
})

test_that("calcSummarySignal combines results from multi-query runs", {
    ql = GRangesList(q1=query, q2=query)
    expect_true(NROW(calcSummarySignal(query, cellMatrix)[[1]])*2 == 
                    NROW(calcSummarySignal(ql, cellMatrix)[[1]]))
    expect_true(NROW(calcSummarySignal(query, cellMatrix)[[2]])*2 == 
                  NROW(calcSummarySignal(ql, cellMatrix)[[2]]))
})

