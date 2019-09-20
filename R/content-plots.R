#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome as a BSgenome object and some ranges on that
#' reference, this function will return a vector of the same length as the
#' granges object, with percent of Cs and Gs.
#' 
#' @param query  A GenomicRanges or GenomicRangesList object with query regions
#' @param ref BSgenome object
#' @export
#' @examples
#' \dontrun{
#' bsg = loadBSgenome('hg19')
#' gcvec = calcGCContent(query, bsg)
#' }
calcGCContent = function(query, ref) {
	# if (!requireNamespace(ref, quietly=TRUE)) {
	# 	message(ref, " package is not installed.")
	# }
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
#' @export
calcGCContentRef = function(query, refAssembly) {
	# if (!requireNamespace(ref, quietly=TRUE)) {
	# 	message(ref, " package is not installed.")
	# }
	ref = loadBSgenome(refAssembly)
	return(calcGCContent(query, ref))
}

#' Plots a density distribution of GC vectors

#' Give results from the \code{calcGCContent} function, this will produce a
#' density plot
#' @param gcvector A vector of GC contents
#' @export
plotGCContent = function(gcvectors) {

	gcdf = as.data.frame(list(gc=gcvec))
	g = ggplot2::ggplot(gcdf, aes(x=gc)) + 
		geom_density() + 
		theme_classic() + 
		xlim(0,1) +
		geom_vline(aes(xintercept=mean(gc)),
            color="blue", linetype="solid", size=0.5)

	return(g)
}