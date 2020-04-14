#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome as a BSgenome object and some ranges on that
#' reference, this function will return a vector of the same length as the
#' granges object, with percent of Cs and Gs.
#' 
#' @param query  A GenomicRanges or GenomicRangesList object with query regions.
#' @param ref Reference genome BSgenome object.
#' @return A numeric vector of list of vectors with the GC percentage of the query regions.
#' @export
#' @examples
#' \dontrun{
#'query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#'GRquery = rtracklayer::import(query)
#'bsg = loadBSgenome('hg19')
#'gcvec = calcGCContent(GRquery, bsg)
#' }
calcGCContent = function(query, ref) {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         ref="BSgenome"))
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcGCContent, ref)
        return(x)
    }
    v = IRanges::Views(ref, query)
    gcvec = apply(Biostrings::alphabetFrequency(v)[,c("C","G")],1, sum)/width(v)
    return(gcvec)
}


#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome as a BSgenome object and some ranges on that
#' reference, this function will return a vector of the same length as the
#' granges object, with percent of Cs and Gs.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes with
#'     \code{getTSSs}.
#' @return A numeric vector or list of vectors with the GC percentage of the query regions.
#' @export
#' @examples
#' \dontrun{
#' query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' GRquery = rtracklayer::import(query)
#' refAssembly = 'hg19'
#' GCcontent = calcGCContentRef(GRquery, refAssembly)
#' } 
calcGCContentRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         refAssembly="character"))
    ref = loadBSgenome(refAssembly)
    return(calcGCContent(query, ref))
}

#' Plots a density distribution of GC vectors

#' Give results from the \code{calcGCContent} function, this will produce a
#' density plot
#' @param gcvectors A numeric vector or list of numeric vectors of GC contents.
#' @return A ggplot object plotting distribution of GC content in query regions.
#' @export
#' @examples
#' numVector = rnorm(400, mean=0.5, sd=0.1)
#' GCplot = plotGCContent(numVector)
#' vecs = list(example1 = rnorm(400, mean=0.5, sd=0.1), example2 = rnorm(600, mean=0.5, sd=0.1))
#' GCplot = plotGCContent(numVector)
#' 
plotGCContent = function(gcvectors) {
    .validateInputs(list(gcvectors=c("numeric","list")))
    gcdf = lapply(gcvectors, as.data.frame)
    # reshape2 is deprecated, but there's no other way to do this easily...
    gcdfReshaped = reshape2::melt(gcdf)
    colnames(gcdfReshaped)[colnames(gcdfReshaped) == "L1"] = "regionSet"
    # plot multiple regionsets if gcvectors is a list
    if (is(gcvectors, "list")) {
        meansdf = aggregate(gcdfReshaped$value, list(gcdfReshaped$regionSet), mean)
        g = ggplot2::ggplot(gcdfReshaped, aes(x=value, colour=regionSet)) +
        geom_density() +
        geom_vline(data=meansdf, aes(xintercept=x, colour=Group.1),
                   linetype="dashed", size=0.5) +
        theme_classic() +
        theme(legend.position = "bottom")
    } else {
        # plot a single regionset
        g = ggplot2::ggplot(gcdfReshaped, aes(x=value)) + 
        geom_density() + 
        geom_vline(aes(xintercept=mean(value)),
                   color="red", linetype="dashed", size=0.5) + 
        theme_classic()
    }    
    g = g + xlim(0,1) + xlab("gc %")
    return(g)
}
