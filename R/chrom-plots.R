
#' Divide regions into roughly equal bins
#'
#' Given a start coordinate, end coordinate, and number of bins to divide, 
#' this function will split the regions into that many bins.
#' Bins will be only approximately the same size, due to rounding.
#' (they should not be more than 1 different).
#'
#' Use case: take a set of regions, like CG islands, and bin them; now you can
#' aggregate signal scores across the bins, giving you an aggregate signal
#' in bins across many regions of the same type.
#'
#' In theory, this just runs on 3 values, but you can run it inside a 
#' data.table j expression to divide a bunch of regions in the same way.
#' @param start The starting coordinate
#' @param end The ending coordinate
#' @param binSize The size of bin to divide the genome into. You must supply
#'     either binSize (priority) or binCount.
#' @param binCount The number of bins to divide. If you do not supply binSize,
#'     you must supply binCount, which will be used to calculate the binSize.
#' @param indicator A vector with identifiers to keep with your bins, in case
#'     you are doing this on a long table with multiple segments concatenated
#'
#' @return
#' A data.table, expanded to nrow = number of bins, with these id columns:
#'     id: region ID
#'     binID: repeating ID (this is the value to aggregate across)
#'     ubinID: unique bin IDs
#' @export
#' @examples
#' Rbins = binRegion(1, 3000, 100, 1000)
#' 
binRegion = function(start, end, binSize=NULL, binCount=NULL, indicator=NULL) {
    .validateInputs(list(start="numeric", end="numeric"))
    if (is.null(binSize) & is.null(binCount)) {
        stop("You must provide either binSize or binCount")
    }
    if (is.null(binSize)) {
        binSize = round(sum(end-start)/binCount)
    }
    binCountByChrom = round((end-start)/binSize)
    binCountByChrom[binCountByChrom==0]=1
    binSizeByChrom = (end-start)/(binCountByChrom)
    breaks = round(unlist(lapply(binCountByChrom, 
                            function(x) seq(from=0, to=x))) * 
                            rep(binSizeByChrom, (binCountByChrom+1)))
    endpoints = cumsum(binCountByChrom + 1) 
    startpoints = c(1, endpoints[-length(endpoints)]+1)

    dataTable = data.table(start=breaks[-endpoints]+1, 
            end=breaks[-startpoints],
            id=rep((seq_along(start)), binCountByChrom),
            binID=unlist(lapply(binCountByChrom, 
                            function(x) seq(from=1, to=x))),
            ubinID=seq_along(breaks[-startpoints]),
            key="id")

    if (!is.null(indicator)){
        idCol = rep(indicator, binCountByChrom)
        dataTable = data.table(idCol, dataTable)
    }
    return(dataTable)
}

#' Bins a BSgenome object.
#'
#' Given a BSgenome object (to be loaded via \code{loadBSgenome}), and a number
#' of bins, this will bin that genome. It is a simple wrapper of the
#' \code{binChroms} function
#' 
#' @param genome A UCSC-style string denoting reference assembly (e.g. 'hg38')
#' @param binCount number of bins per chromosome
#' @return A data.table object showing the region and bin IDs 
#'         of the reference genome.
#' @export
#' @examples
#' \dontrun{
#' binCount = 1000
#' refGenomeBins = binBSGenome("hg19", binCount)
#' }
binBSGenome = function(genome, binCount) {
    .validateInputs(list(genome="character", binCount="numeric"))
    BSG = loadBSgenome(genome)
    chromSizes = seqlengths(BSG)
    return(binChroms(binCount, chromSizes))
}

#' Naively splits a chromosome into bins
#' 
#' Given a list of chromosomes with corresponding sizes, this script will
#' produce (roughly) evenly-sized bins across the chromosomes. It does not
#' account for assembly gaps or the like.
#' 
#' @param binCount number of bins (total; *not* per chromosome)
#' @param chromSizes a named list of size (length) for each chromosome.
#' @return A data.table object assigning a bin ID to each chromosome region.
#' @export
#' @examples 
#' chromSizes = c(chr1=249250621, chr2=243199373, chr3=198022430)
#' cBins = binChroms(1000, chromSizes)
#' 
binChroms = function(binCount, chromSizes) {
    .validateInputs(list(chromSizes="numeric", binCount="numeric"))
    seqnamesColName="chr"
    rangeDT = data.table(chr=names(chromSizes), start=1, end=chromSizes)
    binnedDT = rangeDT[, binRegion(start, end, binCount=binCount,
            indicator=get(seqnamesColName))]
    return(binnedDT)
}


#' Calculates the distribution of a query set over the genome
#' 
#' Returns a data.table showing counts of regions from the query that overlap
#' with each bin.
#' In other words, where on which chromosomes are the ranges distributed?
#' You must provide binned regions. Only the midpoint of each query region is
#' used to test for overlap with the bin regions.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param bins Pre-computed bins (as a GRangesList object) to aggregate
#'    over; for example, these could be genome bins
#' @return A data.table showing where on which chromosomes 
#'    ranges are distributed.
#' @export
#' @examples
#' 
#' chromSizes = getChromSizes("hg19")
#' genomeBins  = getGenomeBins(chromSizes)
#' chromDistribution = calcChromBins(vistaEnhancers, genomeBins)
#' 
#' vistaSftd = GenomicRanges::shift(vistaEnhancers, 100000)
#' vistaSftd2 = GenomicRanges::shift(vistaEnhancers, 200000)
#' calcChromBins(vistaEnhancers, GRangesList(vistaSftd, vistaSftd2))
calcChromBins = function(query, bins) {
    .validateInputs(list(bins=c("GRanges","GRangesList"),
                           query=c("GRanges","GRangesList")))
    if (is(query, "GRangesList"))  {
        # Recurse over each GRanges object
        x = lapply(query, calcChromBins, bins)
        # To accommodate multiple regions, we'll need to introduce a new 'name'
        # column to distinguish them.
        nameList = names(query)
    if(is.null(nameList)) {
        nameList = seq_along(query) # Fallback to sequential numbers
    }
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, vapply(x, nrow, integer(1)))
    return(xb)
    }

    queryDT = grToDt(query)
    
    # This function will just count the number of regions.
    res = calcOLCount(queryDT, bins)

    # order chromosomes by current order.
    res[, chr:=factor(chr, levels=unique(res$chr))]
    return(res)
}

#' Returns the distribution of query over a reference assembly

#' Given a query set of elements (a GRanges object) and a reference assembly
#' (*e.g. 'hg38'), this will aggregate and count the distribution of the query
#' elements across bins of the reference genome. This is a helper function to
#' create features for common genomes. It is a wrapper of
#' \code{calcChromBins}, which is more general.

#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector that will be used to grab chromosome
#'     sizes with \code{getChromSizes}
#' @param binCount Number of bins to divide the chromosomes into
#' @return A data.table showing the distribution of regions across bins of the
#' reference genome.
#' @examples 
#' ChromBins = calcChromBinsRef(vistaEnhancers, "hg19")
calcChromBinsRefSlow = function(query, refAssembly, binCount=3000) {
    .validateInputs(list(refAssembly="character", 
                           query=c("GRanges","GRangesList")))
    # Bin the genome
    chromSizes = getChromSizes(refAssembly)
    binnedDT = binChroms(binCount, chromSizes)
    splitBinnedDT = splitDataTable(binnedDT, "id")
    listGR = lapply(splitBinnedDT, dtToGr, chr="idCol")
    genomeBins =  GRangesList(listGR)
    return(calcChromBins(query, genomeBins))
}


#' Returns the distribution of query over a reference assembly

#' Given a query set of elements (a GRanges object) and a reference assembly
#' (*e.g. 'hg38'), this will aggregate and count the distribution of the query
#' elements across bins of the reference genome. This is a helper function to
#' create features for common genomes. It is a wrapper of
#' \code{calcChromBins}, which is more general.

#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector that will be used to grab chromosome
#'     sizes with \code{getChromSizes}
#' @param binCount Number of bins to divide the chromosomes into
#' @return A data.table showing the distribution of regions across bins of the
#' reference genome.
#' @export
#' @examples 
#' ChromBins = calcChromBinsRef(vistaEnhancers, "hg19")
calcChromBinsRef = function(query, refAssembly, binCount=3000) {
   .validateInputs(list(refAssembly="character",
                           query=c("GRanges","GRangesList")))
    if (is(query, "GRangesList"))  {
        # Recurse over each GRanges object
        x = lapply(query, calcChromBinsRef, refAssembly, binCount)
        # To accommodate multiple regions, we'll need to introduce a new 'name'
        # column to distinguish them.
        nameList = names(query)
       if(is.null(nameList)) {
            nameList = seq_along(query) # Fallback to sequential numbers
        }
        # Append names
        xb = rbindlist(x)
        xb$name = rep(nameList, vapply(x, nrow, integer(1)))
        return(xb)
    }        
   # Bin the genome
    chromSizes = getChromSizes(refAssembly)
    binnedDT = binChroms(binCount, chromSizes)
    queryDT = grToDt(query)
    setnames(binnedDT, "idCol", "chr")
    queryDT[, midpoint:=start + (end-start)]
    # Here I use a non-equi join to get the overlaps
    res = binnedDT[queryDT, .(chr, regionID=ubinID, withinGroupID=x.binID, start=x.start, end=x.end), 
                    on=.(chr, start<=midpoint, end>=midpoint), nomatch=0L][, list(.N), by=list(chr, start, end, regionID, withinGroupID)][order(regionID),]
    res[, chr:=factor(chr, levels=unique(res$chr))]
    return(res)
}



#' Plot distribution over chromosomes
#' 
#' Plots result from \code{genomicDistribution} calculation
#' @param genomeAggregate The output from the genomicDistribution function
#' @param binCount Number of bins (should match the call to
#'     \code{genomicDistribution})
#' @param plotTitle Title for plot.
#' @param ylim Limit of y-axes. Default "max" sets limit to N of biggest bin.
#' @return A ggplot object showing the distribution of the query 
#'     regions over bins of
#' the reference genome.
#' @export
#' @examples
#' agg = data.frame("regionID"=1:5, "chr"=rep(c("chr1"), 5), 
#'                 "withinGroupID"=1:5, "N"=c(1,3,5,7,9))  
#' ChromBins = plotChromBins(agg)
#' 
plotChromBins = function(genomeAggregate, binCount=10000, 
                           plotTitle="Distribution over chromosomes", ylim="max") {
    .validateInputs(list(genomeAggregate=c("data.table","data.frame")))
    
    if ("name" %in% names(genomeAggregate)){
        # It has multiple regions
        # sort the regions labels again
        setkey(genomeAggregate, regionID)
        genomeAggregate[, chr:=factor(chr, levels=unique(genomeAggregate$chr))]
        # and plot
        g = ggplot(genomeAggregate, aes(x=withinGroupID, y=N, 
                                        fill=name, color=name))
    } else {
        # It's a single region
        g = ggplot(genomeAggregate, aes(x=withinGroupID, y=N))
    }
    g = g +
        xlab("Genome") + 
        ylab("Number of regions") +
        geom_bar(stat="identity") + # Spread out to max width
        facet_grid(chr ~ .) + # Place chromosomes one on top of another
        theme_classic() + # Clean up cruft
        theme_blank_facet_label() + # No boxes around labels
        theme(panel.spacing=unit(0, "lines")) + # Reduce whitespace
        theme(strip.text.y=element_text(size=12, angle=0)) + # Rotate labels
        geom_hline(yintercept=0, color="#EEEEEE") + # Light chrom lines
        {if (ylim == "max") {
            scale_y_continuous(breaks = c(max(genomeAggregate$N)),
                               limits = c(0, max(genomeAggregate$N)))
        } else {
            scale_y_continuous(breaks = ylim,
                               limits = c(0, ylim))
        }} +
    scale_x_continuous(breaks=c(0, binCount), labels=c("Start", "End")) +
    theme(plot.title=element_text(hjust=0.5)) + # Center title
    ggtitle(plotTitle) +
    theme(legend.position="bottom")
    return(g)
}

#' Returns bins used in `calcChromBins` function

#' Given a named vector of chromosome sizes, the function returns
#' GRangesList object with bins for each chromosome.

#' @param chromSizes a named list of size (length) for each chromosome.
#' @param binCount number of bins (total; *not* per chromosome), 
#'        defaults to 10,000
#' @return A GRangesList object with bins that separate chromosomes
#'         into equal parts.
#' @export
#' @examples 
#' chromSizes = getChromSizes("hg19")
#' chromBins  = getGenomeBins(chromSizes)
#' 
getGenomeBins = function(chromSizes, binCount=10000) {
  .validateInputs(list(chromSizes="integer"))
  
  binnedDT = binChroms(binCount, chromSizes)
  splitBinnedDT = splitDataTable(binnedDT, "id")
  listGR = lapply(splitBinnedDT, dtToGr, chr="idCol")
  genomeBins =  GRangesList(listGR)
  return(genomeBins)
}
