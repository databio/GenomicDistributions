


#' Internal helper function to calculate distance between neighboring regions
#'
#' @param querydt 
#' @return
neighbordt = function(querydt)  {
    if (length(rownames((querydt))) > 1) {
        endVect = abs(querydt[, diff(end)])
        regionWidth = querydt[, (end-start)]
        widthVect = regionWidth[-1]
        distancesVector = endVect - widthVect
        # neg values represent overlaps between neighboring regions, should set those to 0
        distancesVector[which(distancesVector < 0)] = 0 
        return(distancesVector)
    }
}
  


#' Group regions from the same chromosome together and
#' calculate the distances between neighboring regions. 
#' Distances are then lumped into a numeric vector. 
#'
#' @param query A GRanges object or a GRanges list
#'
#' @return A numeric vector with the distances within neighboring regions
#' @export
calcNeighborDist =  function(query) {
    .validateInputs(list(query=c("GRanges","GRangesList")))
    # lapply if a GRangeslist is provided
    if (is(query, "GRangesList")) {
        dist = lapply(query, calcNeighborDist)
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
  
  

#' Plot the distances between neighboring regions. 
#' 
#'
#' @param dcvec 
#'
#' @return
#' @export
#'
#' @examples
plotNeighborDist = function(dcvec) {
    .validateInputs(list(dcvec=c("numeric","list")))
    distdf = lapply(dcvec, as.data.frame)
    distReshaped = reshape2::melt(distdf, id.vars=NULL)
    colnames(distReshaped)[colnames(distReshaped) == "L1"] = "regionSet"
    if (is(dcvec, "list")) {
        g = ggplot2::ggplot(distReshaped, aes(x=value, fill=regionSet)) +
          geom_density(alpha=0.2) +
          theme_classic() +
          theme(legend.position = "bottom")
    } else {
        g = ggplot2::ggplot(distReshaped, aes(x=value)) +
          geom_density(alpha=0.4) + 
          theme_classic()
    }
    g = g + 
        xlab("Log10 (bp distance)") +
        xlim(0, 10) +
        theme_blank_facet_label() +
        ggtitle("Neighboring regions distance distribution") +
        theme(plot.title = element_text(hjust=0.5)) 
    return(g)
}

