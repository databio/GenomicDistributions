
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
#' 		id: region ID
#' 		binID: repeating ID (this is the value to aggregate across)
#' 		ubinID: unique bin IDs
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
	breaks = round(unlist(sapply(binCountByChrom, function(x) seq(from=0, to=x))) * rep(binSizeByChrom, (binCountByChrom+1)))
	endpoints = cumsum(binCountByChrom + 1) 
	startpoints = c(1, endpoints[-length(endpoints)]+1)

	dt = data.table(start=breaks[-endpoints]+1, 
					end=breaks[-startpoints],
					id=rep((1:length(start)), binCountByChrom),
					binID=unlist(sapply(binCountByChrom, function(x) seq(from=1, to=x))),
					ubinID=1:length(breaks[-startpoints]),
					key="id")

	if (!is.null(indicator)){
		idCol = rep(indicator, binCountByChrom)
		dt = data.table(idCol, dt)
	}
	return(dt)
}

#' Bins a BSgenome object.
#'
#' Given a BSgenome object (to be loaded via \code{loadBSgenome}), and a number
#' of bins, this will bin that genome. It is a simple wrapper of the
#' \code{binChroms} function
#' 
#' @param genome A UCSC-style string denoting reference assembly (e.g. 'hg38')
#' @param binCount number of bins per chromosome
#' @return A data.table object showing the region and bin IDs of the reference genome.
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
#' hg19chromSizes = system.file("data/chromSizes_hg19.RData")
#' cBins = binChroms(1000, hg19chromSizes)
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
#' Returns a data.table showing counts of regions in GR, in the bins
#' In other words, where on which chromosomes are the ranges distributed?
#' You must provide binned regions.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param bins Pre-computed bins (as a GRangesList object) to aggregate
#'     over; for example, these could be genome bins
#' @return A data.table showing where on which chromosomes ranges are distributed.
#' @export
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
			nameList = 1:length(query) # Fallback to sequential numbers
		}
		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, sapply(x, nrow))
		return(xb)
	}

	queryDT = grToDt(query)
	
	# This jExpression will just count the number of regions.
	jExpr = ".N"
	res = BSAggregate(queryDT, bins, jExpr=jExpr)

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
#' @export
#' @examples 
#' query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' GRquery = rtracklayer::import(query)
#' ChromBins = calcChromBinsRef(GRquery, "hg19")
#' 
calcChromBinsRef = function(query, refAssembly, binCount=10000) {
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

#' Plot distribution over chromosomes
#' 
#' Plots result from \code{genomicDistribution} calculation
#' @param genomeAggregate The output from the genomicDistribution function
#' @param binCount Number of bins (should match the call to
#'     \code{genomicDistribution})
#' @param plotTitle Title for plot.
#' @return A ggplot object showing the distribution of the query regions over bins of
#' the reference genome.
#' @export
#' @examples
#' agg = data.frame("regionID"=1:5, "chr"=rep(c("chr1"), 5), "withinGroupID"=1:5, "N"=c(1,3,5,7,9))  
#' ChromBins = plotChromBins(agg)
#' 
plotChromBins = function(genomeAggregate, binCount=10000, 
                         plotTitle="Distribution over chromosomes") {
    .validateInputs(list(genomeAggregate=c("data.table","data.frame")))
    if ("name" %in% names(genomeAggregate)){
		# It has multiple regions
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
		scale_y_continuous(breaks=c(max(genomeAggregate$N)), 
		                   limits=c(0, max(genomeAggregate$N))) +
		scale_x_continuous(breaks=c(0, binCount), labels=c("Start", "End")) +
		theme(plot.title=element_text(hjust=0.5)) + # Center title
		ggtitle(plotTitle) +
		theme(legend.position="bottom")
	return(g)
}