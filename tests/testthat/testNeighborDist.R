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
test_that("calcNeighborDist works", {
  lapply(queryList, function(x) expect_visible(calcNeighborDist(x)))
})

context("result")
test_that("calcNeighborDist returns a result of a proper class", {
  expect_true(is(calcNeighborDist(query), "numeric"))
  expect_true(is(calcNeighborDist(queryList), "list" ))
})

test_that("calcNeighborDist returns the same result for a shifted region set", {
  expect_equal(calcNeighborDist(query), calcNeighborDist(querySftd))
})

# test_that("calcNeighborDist yields a numeric in range 0-10", {
#     x = calcNeighborDist(query)
#     for(i in x){
#         expect_gt(i, 0)    
#         expect_lt(i, 10)
#     }
#})
