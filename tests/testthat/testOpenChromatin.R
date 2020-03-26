# lib
library(data.table)
library(testthat)
library(GenomicDistributions)
# data
exampleCellMatrixFile = system.file("extdata", "example_cell_matrix.txt.gz", package="GenomicDistributions")
cellMatrix = data.table::fread(exampleCellMatrixFile)
f = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
query = rtracklayer::import(f)
querySftd = GenomicRanges::shift(query, 1000000)

# tests
context("general")
test_that("calcOpenSignal works", {
    expect_visible(calcOpenSignal(query, cellMatrix))
    expect_visible(calcOpenSignal(querySftd, cellMatrix))
})

context("result")
test_that("calcOpenSignal returns a proper result of a proper class", {
    expect_true(is(calcOpenSignal(query, cellMatrix), "data.table"))
})

test_that("calcOpenSignal different results for different queries", {
    expect_false(NROW(calcOpenSignal(query, cellMatrix))==NROW(calcOpenSignal(querySftd, cellMatrix)))
})