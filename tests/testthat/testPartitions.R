# lib
library(data.table)
library(testthat)
library(GenomicDistributions)

# data
refs = c("hg38", "hg19")
f = system.file("extdata", "vistaEnhancers.bed.gz", 
                package="GenomicDistributions")
query = rtracklayer::import(f)
querySftd = GenomicRanges::shift(query, 100000)
queryList = GRangesList(q1=query, q2=querySftd)

# tests
context("general")
test_that("calcPartitionsRef works", {
    for(r in refs){
        lapply(queryList, function(x) expect_visible(calcPartitionsRef(x, r)))
    }
})

context("result")
test_that("calcPartitionsRef returns a result of a proper class", {
    expect_true(is(calcPartitionsRef(query, refs[1]), "data.frame"))
})

test_that("calcPartitionsRef returns a result of a proper length", {
    expect_length(calcPartitionsRef(query, refs[1]), 2)
    expect_equal(NROW(calcPartitionsRef(query, refs[1])), 5)
})

test_that("calcPartitionsRef returns different results for different queries", {
    expect_false(all(calcPartitionsRef(query, refs[1])$Freq == 
                     calcPartitionsRef(querySftd, refs[1])$Freq))
})