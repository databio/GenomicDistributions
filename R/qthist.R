
library(gridExtra)

#' Calculate the widths of regions
#' 
#' The length of a genomic region (the distance between the start and end) is called the width
#' When given a query set of genomic regions, this function returns the width
#' @param query A GRanges or GRangesList object with query sets
#' 
#' @export
calcWidth = function(query) { 
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    x = lapply(query, calcWidth)
    return(x) }	
  width(query)
}

#' Plot quantile-trimmed histogram
#' 
#' given the results from \code{calcWidth}, plots a histogram of widths
#' 
#' x-axis breaks for the frequency calculations are based on the "divisions" results from 
#' helper function \code{calcDivisions}
#' 
#' @param widths Results from \code{calcWidths}
#' @param EndBarColor Color for the quantile bars on both ends of the graph (optional)
#' @param MiddleBarColor Color for the bars in the middle of the graph (optional)
#' @param quantile Quantile of data to be contained in each end bar (optional)
#' Must be between 1 and 20 percent, optimal size is under 10.
#' @param bins The number of bins for the histogram to allocate data to. (optional)
#' 
#' @export

plotQTHist = function(widths, EndBarColor = "gray57", MiddleBarColor = "gray27", quantile=NULL, bins=NULL) {
    if (is(widths, "list")) {
        x = lapply(widths, plotQTHist)
        nameswidths = names(widths)
        for (i in seq_along(x)){
            x[[i]] = x[[i]] + ggtitle(nameswidths[i]) }
        return(do.call("grid.arrange", x))
      }
    
    output = calcDivisions(widths, quantile=quantile, bins=bins)
    bins = output[["bins"]]
    quantile = output[["quantile"]]
    div = output[["divisions"]]
    
    colors_vect = c( EndBarColor , rep(MiddleBarColor, (length(div)-3)), EndBarColor) # creates a vector for the colors
    df = cutDists(widths, divisions= div) # calculating a frequency table with the specified divisions
    g = ggplot(df, aes(x=cuts, y=Freq))
    
    g = g +
        geom_bar(stat="identity", fill = colors_vect) + 
        theme_classic() + 
        theme(aspect.ratio=1) + 
        theme_blank_facet_label() +
        ylab("Frequency") +
        xlab("") +
        theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
        theme(plot.title = element_text(hjust = 0.5)) + # Center title
        ggtitle("Quantile Trimmed Histogram") +
        theme(legend.position="bottom") +
        geom_text(aes(label= paste(quantile,"%", sep='')), data=df[c(1,length(df$Freq)),], vjust=-1)
    return(g)
}

#' Internal helper function for \code{plotQTHist}
#' 
#' If the bins or quantiles for the hist. are specified by the user, those are used
#' otherwise, calculates bins based on size of the dataset, and quantile based on bins
#' Returns the divisions that will be used in plotting the histogram

calcDivisions = function(widths, bins=NULL, quantile = NULL){
    # calculating bins
    if(!is.null(bins)){
        b = bins
    } 
    else {
        n = length(widths)
        if (n > 10000) {n = 10000}
        if (n < 500) {n = 500}
        b = round(n^.15 + (n/200) ) # finding number of bins based on the size of dataset
    }
    # calculating quantiles
    if(!is.null(quantile)){
        if(quantile >= 20){
            stop("Quantile must be less than 20. Optimal size is under 10.") }
        q = quantile
    }
    else{
        q = round(25/(b)) # finding the quantile on each side based on number of bins
        q = max(1,q) # minimum on each side is 1%
    }
    quant = unname(quantile(widths, probs = c((q/100), (1-(q/100)))))
    seq_10 = seq(quant[1], quant[2], length = b+1)
    div = c(-Inf, round(seq_10), Inf)
    listOutput <- list("bins"= b,"quantile"= q, "divisions" = div)
    return(listOutput)
}
