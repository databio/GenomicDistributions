# Old, slow version based on GRanges methods
#
# Find the distance to the nearest genomic feature.
# 
# For a given query set of genomic regions, and a given feature set of 
# regions, this function will return the distance for each query region to its
# closest feature. It ignores strand and returns the distance as positive or 
# negative, depending on whether the feature is upstream or downstream.
# 
# This function is similar to the bioconductor distanceToNearest function, but
# returns negative values for downstream distances instead of absolute values.
# This allows you to assess the relative location.
# 
# @param query A GRanges or GRangesList object with query sets
# @param features A GRanges object with features to test distance to
# 
# @return A vector of genomic distances for each query region relative to its 
#         closest feature.
calcFeatureDistBioc = function(query, features) {
    .validateInputs(list(query=x("GRangesList","GRanges")))
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcFeatureDist, features)
        return(x)
    }

    precedeInd = precede(query, features)
    preIndNA = is.na(precedeInd)
    followInd = follow(query, features)
    folIndNA = is.na(followInd)
    preDist = rep(NA, length(query))

    preDist[!preIndNA] = -distance(query[!preIndNA], 
                                   features[precedeInd[!preIndNA]])

    postDist = rep(NA, length(query))
    postDist[!folIndNA] = distance(query[!folIndNA], 
                                   features[followInd[!folIndNA]])

    postHits = -preDist > postDist
    postHitsNA = is.na(postHits)
    dists = preDist
    dists[postHits[!postHitsNA]] = postDist[postHits[!postHitsNA]]
    return(dists)
}

#' Find the distance to the nearest genomic feature
#' 
#' For a given query set of genomic regions, and a given feature set of 
#' regions, this function will return the distance for each query region to its
#' closest feature. It ignores strand and returns the distance as positive or 
#' negative, depending on whether the feature is upstream or downstream
#' 
#' This function is similar to the bioconductor distanceToNearest function, but
#' returns negative values for downstream distances instead of absolute values.
#' This allows you to assess the relative location.
#' 
#' @param query A GRanges or GRangesList object with query sets
#' @param features A GRanges object with features to test distance to
#' 
#' @return A vector of genomic distances for each query region relative to its 
#'     closest feature.
#' @export
#' @examples 
#' vistaSftd = GenomicRanges::shift(vistaEnhancers, 100000)
#' calcFeatureDist(vistaEnhancers, vistaSftd) 
calcFeatureDist = function(query, features) {
    .validateInputs(list(query=c("GRangesList","GRanges")))
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcFeatureDist, features)
        return(x)
    }
    queryDT = grToDt(query)
    featureDT = grToDt(features)
    queryDTs = splitDataTable(queryDT, "chr")
    featureDTs = splitDataTable(featureDT, "chr")
   as.vector(unlist(mapply(queryDTs, featureDTs[names(queryDTs)], 
                           FUN=DTNearest)))
}

# Function uses data.table rolling join to identify the nearest features
# really quickly.
#
# @param DT1 A data.table object to be joined to a second data.table object.
# @param DT2 A second data.table object to join with DT1.
#
# @return A rolling joined data.table object.
DTNearest = function(DT1, DT2) {
    #data.table::set(DT1, j=mid, value=start + round((end-start)/2))
    #data.table::set(DT2, j=mid, value=start + round((end-start)/2))
    if (is.null(DT1)) {
        return(NULL)
    }
    if (is.null(DT2)) {
        return(rep(NA, nrow(DT1)))
    }
    DT1[, mid:=start + round((end-start)/2)]
    DT2[, mid:=start + round((end-start)/2)]
    data.table::setattr(DT1, "sorted", "mid")
    data.table::setattr(DT2, "sorted", "mid")
    DT2[J(DT1), roll="nearest"]
    DT2[J(DT1), start+round((end-start)/2)-mid, roll="nearest"]
}


#' Calculates the distribution of distances from a query set to closest TSS
#' 
#' Given a query GRanges object and an assembly string, this function will grab
#' the TSS list for the given reference assembly and then calculate the distance
#' from each query feature to the closest TSS. It is a wrapper of
#' \code{calcFeatureDist} that uses built-in TSS features for a reference
#' assembly
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes with
#'     \code{getTSSs}.
#' @return A vector of distances for each query region relative to TSSs.
#' @export
#' @examples 
#' calcFeatureDistRefTSS(vistaEnhancers, "hg19")
calcFeatureDistRefTSS = function(query, refAssembly) {
    features = getTSSs(refAssembly)
    return(calcFeatureDist(query, features))
}


# Converts a nucleotide count into a label with abbreviation
# @param x base count
# @return A label with 'kb' or 'mb' appended if appropriate
genomeLabel = function(x) {
    .validateInputs(list(x="numeric"))
    lab = x
    if (abs(x) > 1e6){
        lab = paste0(round(x/1e6), " mb")
    }
    else if (abs(x) > 1e3){
        lab = paste0(round(x/1e3), " kb")
    }
    return(lab)
}


#' Plots a histogram of distances to genomic features
#' 
#' Given the results from \code{featureDistribution}, plots a histogram of
#' distances surrounding the features of interest
#' 
#' @param dists Results from \code{featureDistribution}
#' @param bgdists Background distances. If provided, will plot a background
#'     distribution of expected distances 
#' @param featureName Character vector for plot labels (optional).
#' @param numbers a logical indicating whether the raw numbers should be 
#'     displayed, rather than percentages (optional).
#' @param nbins Number of bins on each side of the center point.
#' @param size Number of bases to include in plot on each side of the 
#'     center point.
#' @param infBins Include catch-all bins on the sides?
#' @param tile Turn on a tile mode, which plots a tiled figure 
#'     instead of a histogram.
#' @param labelOrder -- Enter "default" to order by order of user input (default); 
#' Enter TSS to order by closest proximity to the TSS in descending order (TSS)
#' @return A ggplot2 plot object
#' @export
#' @examples
#' TSSdist = calcFeatureDistRefTSS(vistaEnhancers, "hg19")
#' f = plotFeatureDist(TSSdist, featureName="TSS")
plotFeatureDist = function(dists, bgdists=NULL, featureName="features", 
                           numbers=FALSE, nbins=50, size=100000, 
                           infBins=FALSE, tile=FALSE, labelOrder = "default") {
    df = cutDists(dists, divisions=NULL, nbins, size, infBins)
    
    if(is.list(dists)){
        nplots = length(dists)
    } else {
        nplots = 1
    }

    if (!is.null(bgdists)) {
        bgDistsDF = cutDists(bgDists, divisions=NULL, nbins, size, infBins)
        # bgDistsDF$Freq= scale(bgDistsDF$Freq, center=FALSE)
        bgDistsDF$Freq = (bgDistsDF$Freq / sum(bgDistsDF$Freq)) * 100
        df$bgFreq = rep(bgDistsDF$Freq, nplots)
        df$bgX = rep(seq_len(nrow(bgDistsDF)-1), nplots)
    }

    
    if ("name" %in% names(df)){
        df$name = sortingFunction(df, labelOrder, nbins)
        if (!numbers)
            df$Freq = df[, .(Freq.Per = (Freq / sum(Freq)) * 100), 
                         by = name]$"Freq.Per"
            # It has multiple regions
            df$name = sortingFunction(df, labelOrder, nbins)
            g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
           facet_grid(. ~name)
    } else {
        if (!numbers) 
            df$Freq = (df$Freq / sum(df$Freq)) * 100
            g = ggplot(df, aes(x=cuts, y=Freq))
    }

    if (!is.null(bgdists)) {

    # bgtrack = scale(smooth(bgDistsDF$Freq), center=FALSE)
    g = g + 
        geom_line(stat="identity", aes(x=bgX,y=bgFreq), 
                  color="gray", alpha=1, size=1.5) + 
        geom_bar(stat="identity", aes(x=cuts,y=bgFreq), fill="gray", alpha=0.8)
    }

    # find midpoint
    midx = nrow(df)/2/nplots
    barcount = nrow(df)/nplots
    minlabel = genomeLabel(-size)
    maxlabel = genomeLabel(size)
    edgeLabels = c(minlabel, rep("", barcount-2), maxlabel)

    if (tile) {
        if (!"name"  %in% names(df)) {
            df$name = "Region set"
    }

    ncuts = length(unique(df$cuts))
    xs = rep(seq_len(ncuts), nplots)
    g = ggplot(df) + 
        geom_raster(aes(x=xs, y=name, fill=Freq)) +
        scale_fill_gradient(low="navy", high="orange") +
        geom_point(aes(x=midx, y=0.5), color="black", 
                   size=2, shape=17, alpha=0.8) + 
        theme_classic() + 
        labs(fill=ifelse(numbers,"Counts","Frequency (%)")) +
        theme(legend.position="bottom") + 
        xlab(paste("Distance to", featureName)) +
        theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5)) +
        scale_x_continuous(breaks=c(1, ncuts), labels=c(minlabel, maxlabel))
        return(g)
    }

    g = g +
        geom_bar(data=df, stat="identity", fill="darkblue", alpha=0.7) + 
        geom_point(aes(x=midx, y=0), color="tan2", size=2, 
                   shape=17, alpha=0.8) +
        guides(fill=FALSE) + # remove legend for geom_point
        theme_classic() + 
        theme(aspect.ratio=1) + 
        theme_blank_facet_label() + 
        xlab(paste("Distance to", featureName)) +
        ylab(ifelse(numbers,"Counts","Frequency (%)")) +
        # theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5))
        theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5)) + 
        theme(plot.title = element_text(hjust = 0.5)) + # Center title
        ggtitle(paste("Distribution relative to", featureName)) +
        theme(legend.position="bottom") + 
        theme(panel.spacing.x=unit(1, "lines")) + 
        scale_x_discrete(labels=edgeLabels) +
        scale_x_discrete(labels=edgeLabels, expand=expansion(mult=0.035))

    return(g)
}

# Internal helper function for \code{plotFeatureDist}
#
# @param dists A vector of genomic distances.
# @param divisions A vector of bin sizes to divide the dists into.
# @param nbins Number of bins on each side of the center point.
# @param size Number of bases to include in plot on 
# each side of the center point.
# @param infBins Include catch-all bins on the sides?
# @return A data.frame of the table of the frequency of dists in divisions.

sortingFunction = function(df, labelOrder = "default", nbins = 50){
    if(labelOrder == "default"){
        orderedLabels = unique(df$name)
        orderedNames = factor(df$name, levels = orderedLabels)
        return(orderedNames)
    }
    if (labelOrder == "TSS"){
        temp = df[df$cuts == "0 to 2,000",]
        orderedLabels = temp[order(temp$Freq, decreasing = TRUE), ]$name
        orderedNames = factor(df$name, levels = orderedLabels)
        return(orderedNames)
        
        # To make it based on the bins, I wasn't sure how to alternate by nbins + 1 without a for loop...
        # middleVal = nbins + 1
        # for (i in df$Freq){
        # orderedLabels = type(df[rev(order(df$Freq[i])),]$name #row index )
        # i <= middleVal + middleVal
        }

    }

}


# Internal helper function for \code{plotFeatureDist}
#
# @param dists A vector of genomic distances.
# @param divisions A vector of bin sizes to divide the dists into.
# @param nbins Number of bins on each side of the center point.
# @param size Number of bases to include in plot on 
# each side of the center point.
# @param infBins Include catch-all bins on the sides?
# @return A data.frame of the table of the frequency of dists in divisions.
cutDists = function(dists, divisions=NULL, nbins=50, 
                    size=100000, infBins=TRUE) {
    if (is.null(divisions)) {
        poscuts = seq(0, size, by=size/nbins)
        divisions = sort(unique(c(-poscuts, poscuts)))
        if (infBins) {
            divisions = c(-Inf, divisions, Inf)
        }
    }
    if (is.list(dists)) {
        x = lapply(dists, cutDists, divisions)

        # To accommodate multiple lists, we'll need to introduce a new 'name'
        # column to distinguish them.
        nameList = names(dists)
        if(is.null(nameList)) {
            nameList = seq_along(dists) # Fallback to sequential numbers
        }

    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, vapply(x, nrow, integer(1)))

    return(xb)
    }

    labels = labelCuts(sort(divisions), collapse=" to ", infBins=infBins)
    cuts = cut(dists, divisions, labels)
    df = as.data.frame(table(cuts))
    return(df)
}

