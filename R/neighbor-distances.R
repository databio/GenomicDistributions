


#' Group regions from the same chromosome together and
#' calculate the distances between neighboring regions. 
#' Distances are then lumped into a numeric vector. 
#'
#' @param query A GRanges or GRangesList object.
#'
#' @return A numeric vector or list with different vectors containing the
#'  distances within neighboring regions.
#' @export
#' @examples 
#' dist = calcNeighborDist(vistaEnhancers)
calcNeighborDist =  function(query) {
    .validateInputs(list(query=c("GRanges","GRangesList")))
    # lapply if a GRangeslist is provided
    if (is(query, "GRangesList")) {
        dist = lapply(query, calcNeighborDist)
        namelist = names(query)
        if (is.null(namelist)) {
            newnames = seq_along(query)
            namelist = newnames
            # Append names
            names(dist) = namelist
        }
        return(dist)
    }
    querydt = grToDt(sort(query))
    querydts = splitDataTable(querydt, "chr")
    distanceVectors = lapply(querydts, neighbordt)
    d = as.vector(unlist(distanceVectors))
    # remove overlaps for log10 transformation
    dcvec = d[!(d == "0")] 
    dcvec = log10(dcvec)
    return(dcvec)
}
  
#' Internal helper function to calculate distance between neighboring regions.
#'
#' @param querydt A data table with regions grouped according to
#' chromosome.
#' @return A numeric vector with the distances in bp 
neighbordt = function(querydt)  {
    # there should be at least 2 regions for each chr
    if (length(rownames((querydt))) > 1) {
        endVect = abs(querydt[, diff(end)])
        regionWidth = querydt[, (end-start)]
        widthVect = regionWidth[-1]
        distancesVector = endVect - widthVect
        # neg values represent overlaps between neighbor regions, set those to 0
        distancesVector[which(distancesVector < 0)] = 0 
        return(distancesVector)
  }
}
  

#' Plot the distances between neighboring regions.The distance in the 
#' x axis is log10 transformed for ease of comparison between 
#' different regionsets and to account for outliers. 
#' 
#' @param dcvec A numeric vector or list with vectors containing distances 
#' between neighbor regions. Produced by \code{calcNeighborDist}
#'
#' @return A ggplot density object showing the distribution of
#' log10 transformed distances.  
#' @export
#' @examples
#' numVector = rnorm(400, mean=5, sd=0.1)
#' d = plotNeighborDist(numVector)
plotNeighborDist = function(dcvec) {
    .validateInputs(list(dcvec=c("numeric","list")))
    distdf = lapply(dcvec, as.data.frame)
    distReshaped = reshape2::melt(distdf, id.vars=NULL)
    colnames(distReshaped)[colnames(distReshaped) == "L1"] = "regionSet"
    if (is(dcvec, "list")) {
        g = ggplot2::ggplot(distReshaped, aes(x=value, fill=regionSet, colour=regionSet)) +
          geom_density(alpha=0.5) +
          theme_classic() +
          theme(legend.position = "bottom")
    } else {
        g = ggplot2::ggplot(distReshaped, aes(x=value)) +
          geom_density(alpha=0.4) + 
          theme_classic()
    }
    g = g + 
        xlab(expression(log[10]*("bp distance"))) +
        xlim(0, 10) +
        theme(aspect.ratio=1) +
        theme_blank_facet_label() +
        ggtitle("Neighboring regions distance distribution") +
        theme(plot.title = element_text(hjust=0.5)) 
    return(g)
}

