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
#'bsg = loadBSgenome('hg19')
#'gcvec = calcGCContent(query, bsg)
#' 
calcGCContent = function(query, ref) {
	# if (!requireNamespace(ref, quietly=TRUE)) {
	# 	message(ref, " package is not installed.")
	# }
	if (!(is(query, "GRanges") || is(query, "GRangesList" ))) {
	  stop("query should be a GRanges object or GRanges list. Check object class.")
	}
  if (!(is(ref, "BSgenome"))) {
    stop("ref should be a BSgenome object.")
  }
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
calcGCContentRef = function(query, refAssembly) {
	# if (!requireNamespace(ref, quietly=TRUE)) {
	# 	message(ref, " package is not installed.")
	# }
  if (!(is(query, "GRanges") || is(query, "GRangesList" ))) {
    stop("query should be a GRanges object or GRanges list. Check object class.")
  }  
  if (!(is(refAssembly, "character"))) {
    stop("refAssembly should be a character vector specifying the reference genome.")
  }
	ref = loadBSgenome(refAssembly)
	return(calcGCContent(query, ref))
}

#' Plots a density distribution of GC vectors

#' Give results from the \code{calcGCContent} function, this will produce a
#' density plot
#' @param gcvector A numeric vector or list of GC contents.
#' @return A ggplot object plotting distribution of GC content in query regions.
#' @export
plotGCContent = function(gcvectors) {
  if (!(is(gcvectors, "list") || is(gcvectors, "numeric"))) {
    stop("gcvectors should be a numeric vector or list of vectors. Check object class")
  }
	gcdf = as.data.frame(list(gc=gcvectors))
	gcdfReshaped = reshape2::melt(gcdf) 
	# plot multiple regionsets if gcvectors is a list
	if (is(gcvectors, "list")) {
	  meansdf = aggregate(gcdfReshaped$value, list(gcdfReshaped$variable), mean)
	  g = ggplot2::ggplot(gcdfReshaped, aes(x=value, colour=variable)) +
	    geom_density() +
	    geom_vline(data=meansdf, aes(xintercept=x, colour=Group.1),
	               linetype="solid", size=0.5) +
	    theme_classic() +
	    theme(legend.position = "bottom")
	} else {
	  # plot a single regionset
	  g = ggplot2::ggplot(gcdfReshaped, aes(x=value)) + 
	    geom_density() + 
	    geom_vline(aes(xintercept=mean(value)),
	               color="red", linetype="solid", size=0.5) + 
	    theme_classic()
	}    
	g = g +   
	 xlim(0,1) 

	return(g)
}
