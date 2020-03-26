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
#' @references \url{http://github.com/databio/GenomicDistributions}
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#'             seqnames granges
#' @importFrom ggplot2 ggplot aes facet_grid geom_jitter geom_line
#'             geom_bar theme_classic xlab ylab geom_hline ylim scale_color_discrete
#'             scale_x_discrete scale_y_discrete scale_fill_brewer scale_color_manual
#'             scale_x_continuous ggtitle geom_vline scale_fill_discrete xlim
#'             scale_color_brewer theme element_blank unit element_text geom_density
#'             geom_point
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder copy is.data.table
#' @importFrom gridExtra grid.arrange
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
