#' Plot distance width of the genomic ranges in query vs. their distance
#' to TSS.
#
#' @param dists Numeric vector or list of distances of features to TSS.
#' @param widths Numeric vector or list of feature widths.
#' @return  A ggplot object
#' 
#' @export
#' @examples
#' \dontrun{
#' query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' GRquery = rtracklayer::import(query)
#' dists = calcFeatureDistRefTSS(query, "hg19")
#' widths = calcWidth(query)
#' plotDistWidth(dists, widths)
#' }
plotDistWidth = function(dists, widths){
  if(is(dists, "list") && is(widths, "list")){
    
    # check that names of the lists are matching
    if(is.null(names(dists)) && is.null(names(widths))) {
      nameList = seq_along(dists) # Fallback to sequential numbers
    } else if (identical(names(dists), names(widths))){
      nameList = names(dists)
    } else {
      stop("Dists and widths lists must have matching list names.")
    }
    
    # Make the table for ggplot
    unlistDist = unlist(dists)
    unlistWidths = unlist(widths)
    listName = rep(nameList, vapply(dists, length, integer(1)))
    plotTable = data.frame(TSSdist = unlistDist, 
                           widths = unlistWidths, 
                           name = listName)
    
  } else if (is(dists, "numeric") && is(widths, "numeric")) {
    plotTable = data.frame(TSSdist = dists, widths = widths)
  } else {
    stop("Objects dists and widths must be both numeric or both lists.")
  }
  
  # make scatter plot
  p = ggplot(plotTable, aes(TSSdist, widths))
  p = p + geom_point(alpha = 0.2) +
    geom_vline(xintercept = 0) +
    theme_classic()
  
  # if lists - make facets
  if ("name" %in% names(plotTable)){
    p = p + facet_grid(. ~name)
  }
  return(p)
}
