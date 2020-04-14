# lib
library(data.table)
library(testthat)
library(GenomicDistributions)
# data
exampleCellMatrixFile = system.file("extdata", "example_cell_matrix.txt", 
                                    package="GenomicDistributions")
cellMatrix = data.table::fread(exampleCellMatrixFile)
f = system.file("extdata", "vistaEnhancers.bed.gz", 
                package="GenomicDistributions")
query = rtracklayer::import(f)
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
    expect_true(is(calcOpenSignal(query, cellMatrix), "data.table"))
})

test_that("calcOpenSignal returns different results for different queries", {
    expect_false(NROW(calcOpenSignal(query, cellMatrix)) == 
                     NROW(calcOpenSignal(querySftd, cellMatrix)))
})

test_that("calcOpenSignal combines results from multi-query runs", {
    ql = GRangesList(q1=query, q2=query)
    expect_true(NROW(calcOpenSignal(query, cellMatrix))*2 == 
                    NROW(calcOpenSignal(ql, cellMatrix)))
})

