# COMMENTED OUT DUE TO BSgenome.Hsapiens.UCSC.<genome>.masked PACKAGES 
# DEPENDANCIES, WHICH ARE NOT INCLUDED IN REQUIREMENTS DUE TO SIZE

# # lib
# library(testthat)
# 
# # data
# featureFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
# feats = rtracklayer::import(featureFile)
# refs = c("hg38", "hg19")
# 
# # tests
# context("general")
# test_that("calcGCContent works", {
#     for(r in refs){
#         expect_visible(calcGCContentRef(feats, r))
#     }
# })
# 
# context("result")
# test_that("calcGCContent yields results of proper length", {
#     expect_equal(length(calcGCContentRef(feats, "hg19")), length(feats))
# })
# 
# test_that("calcGCContent yields a numeric result", {
#     expect_true(is.numeric(calcGCContentRef(feats, "hg19")))
# })
# 
# test_that("calcGCContent yields a numeric in range 0-1", {
#     x = calcGCContentRef(feats, "hg19")
#     for(i in x){
#         expect_gt(i, 0)    
#         expect_lt(i, 1)
#     }
# })
