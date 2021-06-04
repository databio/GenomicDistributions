


#' Group regions from the same chromosome together and
#' calculate the distances of a region to its upstream and
#' downstream neighboring regions.  
#' Distances are then lumped into a numeric vector. 
#'
#' @param query A GRanges or GRangesList object.
#' @param correctedRef A string indicating the reference genome
#' to use if distances are corrected for the number of 
#' regions in a regionSet. 
#' 
#' @return A numeric vector or list with different vectors containing the
#'  distances of regions to their upstream/downstream neighbors.
#' @export
#' @examples 
#' dist = calcNeighborDist(vistaEnhancers)
calcNeighborDist =  function(query, correctRef="None") {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         correctRef=c("character")))
    # lapply if a GRangeslist is provided
    if (is(query, "GRangesList")) {
        dist = lapply(query,
                      function(x){calcNeighborDist(x, correctRef = correctRef)})
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
    # remove overlaps
    dcvec = d[!(d == "0")]
    # Correct for number of regions
    if (!correctRef=="None") {
        chromSizes = getChromSizes(correctRef)
        genomelen = sum(chromSizes)
        meanWidth = mean(calcWidth(query))
        expectedDist = genomelen/nrow(querydt) - meanWidth
        correctedDist = log10(dcvec/expectedDist)
        return(correctedDist)
    # If we just want to look at the raw neighbor distances
    } else {
        return(dcvec)
    }
}
  
#' Internal helper function to calculate distance 
#' between neighboring regions.
#'
#' @param querydt A data table with regions grouped according to
#' chromosome.
#' @return A numeric vector with the distances in bp 
neighbordt = function(querydt)  {
    # there should be at least 2 regions for each chr
    if (nrow(querydt) > 1) {
        endVect = abs(querydt[, diff(end)])
        regionWidth = querydt[, (end-start+1)]
        distancesVector = endVect - regionWidth[-1]
        # neg values represent overlaps between neighbor regions, set those to 0
        distancesVector[which(distancesVector < 0)] = 0 
        return(distancesVector)
  }
}


#' Group regions from the same chromosome together and
#' compute the distance of a region to its nearest neighbor. 
#' Distances are then lumped into a numeric vector. 
#'
#' @param query A GRanges or GRangesList object.
#' @param correctedRef A string indicating the reference genome
#' to use if Nearest neighbor distances are corrected for the 
#' number of regions in a regionSet. 
#'
#' @return A numeric vector or list of vectors containing the
#'  distance of regions to their nearest neighbors.
#' @export
#' @examples 
#' Nneighbors = calcNearestNeighbors(vistaEnhancers)
calcNearestNeighbors = function(query, correctRef="None") {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         correctRef=c("character")))
    # lapply if a GRangeslist is provided
    if (is(query, "GRangesList")) {
        dist = lapply(query,
                      function(x){calcNearestNeighbors(x, correctRef = correctRef)})
        namelist = names(query)
        if (is.null(namelist)) {
            newnames = seq_along(query)
            namelist = newnames
            # Append names
            names(dist) = namelist 
        }
        return(dist)
    }
    # Calculate nearest neighbors in a vectorized manner
    dist = calcNeighborDist(query)
    upstream = dist[-length(dist)]
    downstream = dist[-1]
    dt = data.table(i=upstream, j=downstream)
    pairmins = dt[, pmin(i, j)]
    # First and last distances are default nearest neighbors
    nNeighbors = c(dist[1], pairmins, dist[length(dist)])
    # Correct for number of regions
    if (!correctRef=="None") {
        chromSizes = getChromSizes(correctRef)
        genomelen = sum(chromSizes)
        meanWidth = mean(calcWidth(query))
        expectedDist = genomelen/length(query) - meanWidth
        correctedDist = log10(nNeighbors/expectedDist)
        return(correctedDist)
    } else {
        return(nNeighbors)
    }
}

#' Plot the distances from regions to their upstream/downstream neighbors
#' or nearest neighbors. Distances can be passed as either raw bp or
#' corrected for the number of regions (log10(obs/exp)), but this has
#' to be specified in the function parameters. 
#' 
#' @param dcvec A numeric vector or list of vectors containing distances 
#' to upstream/downstream neighboring regions or to nearest neighbors. 
#' Produced by \code{calcNeighborDist} or \code{calcNearestNeighbors}
#' @param correctedDist A logical indicating if the plot axis should
#' be adjusted to show distances corrected for the number of regions
#' in a regionset.
#' @param Nneighbors A logical indicating whether legend should be adjusted
#' if Nearest neighbors are being plotted. Default legend shows distances
#' to upstream/downstream neighbors.  
#'
#' @return A ggplot density object showing the distribution of
#' raw or corrected distances.  
#' @export
#' @examples
#' numVector = rnorm(400, mean=5, sd=0.1)
#' d = plotNeighborDist(numVector)
plotNeighborDist = function(dcvec, correctedDist=FALSE,
                            Nneighbors=FALSE) {
    .validateInputs(list(dcvec=c("numeric","list")))
    # if input is list, convert it to a data frame with 
    # value and region set name, if input is vector - make a single
    # columns data.frame
    if (is(dcvec, "list")) {
        nameList = names(dcvec)
        vectorLengths = unlist(lapply(dcvec, length))
        distReshaped = data.frame(value = unlist(dcvec),
                                  regionSet = rep(nameList, vectorLengths))
        g = ggplot2::ggplot(distReshaped, aes(x=value,
                                              fill=regionSet,
                                              colour=regionSet)) +
          geom_density(alpha=0.4)
    } else {
        distReshaped = data.frame(value = dcvec)
        g = ggplot2::ggplot(distReshaped, aes(x=value)) +
          geom_density() 
    }
    if (correctedDist==TRUE) {
        g = g + 
          xlab(expression(log[10](over(Obs, Exp)))) +
          geom_vline(xintercept = 0, linetype="dashed") +
          ggtitle("Corrected neighboring regions distance distribution") 
    } else {
        g = g + 
          xlab(expression("bp distance")) +
          scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                        labels = scales::trans_format("log10", 
                                                  scales::math_format(10^.x))) +
          ggtitle("Neighboring regions distance distribution") 
    }
    g = g + 
      theme_classic() +
      theme(aspect.ratio=1,
            plot.title = element_text(hjust=0.5),
            legend.position = "bottom") +
      theme_blank_facet_label()
  
  # Adjust legend if plotting nearest neighbors 
    if (Nneighbors==TRUE){
        g = g + 
          labs(fill="regionSet Nneighbors", 
               colour="regionSet Nneighbors")
    }
    return(g)
}

