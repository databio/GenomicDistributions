# PACKAGE DOCUMENTATION
#' Produces summaries and plots of features distributed across genomes
#' 
#' If you have a set of genomic ranges, the GenomicDistributions R package can
#' help you with some simple visualizations. Currently, it can produce two kinds
#' of plots: First, the chromosome distribution plot, which visualizes how your
#' regions are distributed over chromosomes; and second, the feature
#' distribution plot, which visualizes how your regions are distributed relative
#' to a feature of interest, like Transcription Start Sites (TSSs).
#'
#'
#' @docType package
#' @name GenomicDistributions
#' @author Nathan C. Sheffield
#'
#' @references \url{http://github.com/nsheff}
#' @import GenomicRanges
#' @import ggplot2
#' @import data.table
#' @importFrom methods is
NULL


# Because of some issues with NOTEs on R CMD CHeck and CRAN submission,
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1") {
	utils::globalVariables(c(
	"cuts", "mid", "J", "chr", "N",
	"regionID", "x", "name", "BSFilter"))
}
