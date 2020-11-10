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
#' @import ggplot2
#' @importFrom GenomicRanges GRanges GRangesList elementMetadata strand
#'             seqnames granges
#' @importFrom data.table ":=" setDT data.table setkey fread setnames 
#'             setcolorder rbindlist setattr setorder copy is.data.table
#'             tstrsplit as.data.table foverlaps
#' @importFrom reshape2 melt
#' @importFrom IRanges IRanges Views
#' @importFrom Biostrings alphabetFrequency
#' @importFrom methods is
#' @importFrom utils installed.packages getAnywhere data globalVariables

NULL

# You can either use 'import X' or 'importFrom X abcdefg'. importFrom 
# is better practice,
# but for ggplot2 we were simply importing so many functions that it makes 
# sense to just import the whole package
# @importFrom ggplot2 ggplot aes facet_grid geom_jitter geom_line
#             geom_bar theme_classic xlab ylab geom_hline ylim 
#             scale_color_discrete scale_x_discrete scale_y_discrete 
#             scale_fill_brewer scale_color_manual scale_x_continuous
#             ggtitle geom_vline scale_fill_discrete xlim
#             scale_color_brewer theme element_blank unit 
#             element_text geom_density geom_point guides geom_col 
#             theme_bw scale_fill_manual


# Because of some issues with NOTEs on R CMD check and CRAN submission,
# (see here: http://stackoverflow.com/questions/9439256/)
# I have to register stuff used in data.table as non-standard evaluation,
# in order to pass some R CMD check NOTES.
if(getRversion() >= "2.15.1") {
    utils::globalVariables(c(
    "cuts", "mid", "J", "chr", "N", "regionID", "x", "name", "BSFilter", 
    "start", "end", "findOverlaps", "queryHits", "subjectHits", "buildJ",
    "seqlengths", "IRanges", "seqlengths", "reduce", "seqlevels", "follow",
    "trim", "error", "nlist", "aggregate", "median",  "bgDists", "Freq", "bgX",
    "bgFreq", "value", "regionSet", "Group.1", "cellType", "spaceLabel", 
    "signal", "group", "medianBar", "partition", "Freq", "Freq", "cumsize", 
    "frif", "aggregate", "withinGroupID", "lowerCaseTissue", "boxplot.stats", 
    "median", "barplot", "legend", "promoters", "seqlevels", "width", 
    "precede", "elementMetadata", ".N", ".SD", "colorRampPalette", "count", 
    "countOverlaps", "distance", "elementMetadata<-", "elementNROWS", 
    "expected", "log10OE", "pintersect", "plot_labels", "query", 
    "regionGroupID", "seqlevels<-", "size", "tableCount", "V1", "queryPeak", 
    "xid", "yid", "na.omit", "peakName", "mixedVar",
    "cellTypeMetadata", "tissueType", "boxStats",
    "tissue", ".", "Percent", "Var1"))
}


