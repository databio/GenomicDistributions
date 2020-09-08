# lib
library(data.table)
library(testthat)
library(GenomicDistributions)
# data
cellMatrix = exampleOpenSignalMatrix_hg19
query = vistaEnhancers
querySftd = GenomicRanges::shift(query, 100000)
queryList = GRangesList(q1=query, q2=querySftd)

# tests
context("general")
test_that("calcOpenSignal works", {
    expect_visible(calcOpenSignal(query, cellMatrix))
    expect_visible(calcOpenSignal(querySftd, cellMatrix))
})

test_that("calcOpenSignal works with multiple queries", {
    expect_visible(calcOpenSignal(queryList, cellMatrix))
})

context("result")
test_that("calcOpenSignal returns a result of a proper class", {
    expect_true(is(calcOpenSignal(query, cellMatrix), "list"))
    expect_true(is(calcOpenSignal(query, cellMatrix)[[1]], "data.table"))
    expect_true(is(calcOpenSignal(query, cellMatrix)[[2]], "data.frame"))
})

test_that("calcOpenSignal returns different results for different queries", {
    expect_false(identical(calcOpenSignal(query, cellMatrix)[[1]], 
                           calcOpenSignal(querySftd, cellMatrix)[[1]]))
})

test_that("calcOpenSignal combines results from multi-query runs", {
    ql = GRangesList(q1=query, q2=query)
    expect_true(NROW(calcOpenSignal(query, cellMatrix)[[1]])*2 == 
                    NROW(calcOpenSignal(ql, cellMatrix)[[1]]))
    expect_true(NROW(calcOpenSignal(query, cellMatrix)[[2]])*2 == 
                  NROW(calcOpenSignal(ql, cellMatrix)[[2]]))
})

