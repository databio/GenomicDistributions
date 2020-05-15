#' Calculate the widths of regions
#' 
#' The length of a genomic region (the distance between the start and end) is called the width
#' When given a query set of genomic regions, this function returns the width
#' @param query A GRanges or GRangesList object with query sets
#' @return A vector of the widths (end-start coordinates) of GRanges objects.
#' @export
#' @examples
#' TSSdist = calcFeatureDistRefTSS(vistaEnhancers, "hg19")
#' plotFeatureDist(TSSdist, featureName="TSS")
calcWidth = function(query) { 
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcWidth)
        return(x) } 
    width(query)
}


#' Plot quantile-trimmed histogram
#' 
#' Given the results from \code{calcWidth}, plots a histogram with outliers trimmed.
#' 
#' x-axis breaks for the frequency calculations are based on the "divisions" 
#' results from helper function \code{calcDivisions}.
#' 
#' @param x Data values to plot
#' @param EndBarColor Color for the quantile bars on both ends of the graph
#'     (optional)
#' @param MiddleBarColor Color for the bars in the middle of the graph
#'     (optional)
#' @param quantile Quantile of data to be contained in each end bar (optional)
#' Quantiles must be under .2, optimal size is under .1
#' @param bins The number of bins for the histogram to allocate data to.
#'     (optional)
#' @param indep logical value which returns a list of plots that have had their
#'     bins calculated independently; the normal version will plot them on the 
#'     same x and y axis.
#' @param numbers a logical indicating whether the raw numbers should be 
#'     displayed, rather than percentages (optional).
#' @return A ggplot2 plot object
#' @export
#' @examples
#' plotQTHist(runif(500)*1000)
#' plotQTHist(list(q1=runif(500)*1000, q2=runif(500)*1000))
plotQTHist = function(x, EndBarColor = "gray57", MiddleBarColor = "gray27",
    quantile=NULL, bins=NULL, indep=FALSE, numbers=FALSE) {
    if (indep) {
        if (is(x, "list") | is(x, "List")) {
            x = lapply(x, plotQTHist)
            namesx = names(x)
            for (i in seq_along(x)){
                x[[i]] = x[[i]] + ggtitle(namesx[i])
            }
        return(x)
        # you can use grid.arrange like this to plot these           
        # do.call("grid.arrange", x)
        }
    }
    output = calcDivisions(x, quantile=quantile, bins=bins)
    # if all x are the same - recalculate divisions
    divisionCheck = output[["divisions"]]
    if (length(divisionCheck) > length(unique(divisionCheck))){
      if (length(unique(divisionCheck)) == 3){
        output[["divisions"]] = c(-Inf, divisionCheck[2], divisionCheck[2]+1, Inf)
        output[["bins"]] = 1
      } else {
        output[["divisions"]] = unique(divisionCheck)
        output[["bins"]] = (length(unique(divisionCheck)) - 3)
      }
    }
    if(is(x, "List")){
        x = as.list(x)
    }
    if(is.list(x)){
        nplots = length(x)
    } else {
        nplots = 1
    }

    df = cutDists(x, divisions=output[["divisions"]])
    if ("name" %in% names(df)){
        if (!numbers)
            df$Freq = df[, .(Freq.Per = (Freq / sum(Freq)) * 100), by = name]$"Freq.Per"

        g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
            facet_wrap(. ~name)
    } else {
        g = ggplot(df, aes(x=cuts, y=Freq))
    }
    # Create a vector for the colors
    colors_vect = c(EndBarColor ,
        rep(MiddleBarColor, (length(output[["divisions"]])-3)), EndBarColor)
    colors_vect = rep(colors_vect, nplots)

    nbars = output[["bins"]]+2
    qbaridx = sort(c(seq(1, nbars*nplots, by=nbars),
            seq(nbars, nbars*nplots, by=nbars)))
  
    g = g +
        geom_bar(stat="identity", fill = colors_vect) + 
        theme_classic() + 
        theme(aspect.ratio=1) + 
        theme_blank_facet_label() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) +
        theme(plot.title = element_text(hjust = 0.5)) + # Center title
        ggtitle("Quantile Trimmed Histogram") +
        theme(legend.position="bottom") +
        geom_text(aes(label= paste((output[["quantile"]]*100),"%", sep='')),
            data=df[qbaridx,], hjust=-1, angle=90, size=2.5)
    return(g)
}


# Internal helper function for \code{plotQTHist}
# 
# If the bins or quantiles for the hist. are specified by the user, those are 
# used otherwise, calculates bins based on size of the dataset, and quantile 
# based on bins.
#
# @param x A vector of GRanges x.
# @return A list of the divisions that will be used in plotting the histogram. 
# @examples
# calcDivisions(runif(500)*1000)
calcDivisions = function(x, bins=NULL, quantile = NULL){
    if(is.list(x))
        x=unlist(x)
        
    # calculating bins
    if(!is.null(bins)){
        b = bins
    } 
    else {
        n = length(x)
        if (n > 10000) {n = 10000}
        if (n < 500) {n = 500}
        # finding number of bins based on the size of dataset
        b = round(n^.15 + (n/200) )
    }
    # calculating quantiles
    if(!is.null(quantile)){
        if(quantile >= .2){
            stop("Quantile must be less than .2, Optimal size is under .1") }
        q = quantile
    }
    else{
        # finding the quantile on each side based on number of bins
        q = round(25/(b))/100
        # minimum on each side is 1%
        q = max(.01,q)
    }
    quant = unname(quantile(x, probs = c((q), (1-(q)))))
    seq_10 = seq(quant[1], quant[2], length = b+1)
    div = c(-Inf, round(seq_10), Inf)
    listOutput <- list("bins"= b,"quantile"= q, "divisions" = div)
    return(listOutput)
}
