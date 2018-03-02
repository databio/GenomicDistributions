
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
binRegion = function(start, end, binSize=NULL, binCount=NULL, indicator=NULL) {
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
#' Given a BSgenome object (possibly loaded via \code{loadBSgenome}), and a
#' number of bins for each chromosome, this will 
#' @param genome A UCSC-style string denoting reference assembly (e.g. 'hg38')
#' @param binCount number of bins per chromosome
#' @export
binBSGenome = function(genome, binCount) {
	BSG = loadBSgenome(genome)
	chromSizes = seqlengths(BSG)
	binChroms(binCount, chromSizes)
}

#' Splits a chromosome into binSize

#' Given a list of chromosomes with corresponding sizes, this script will produce
#' (roughly) evenly-sized bins across the chromosomes
#' @param binCount number of bins per chromosome
#' @param chromSizes a named list of size (length) for each chromosome
#' @export
binChroms = function(binCount, chromSizes) {
	seqnamesColName="chr"
	rangeDT = data.table(chr=names(chromSizes), start=1, end=chromSizes)
	binnedDT = rangeDT[, binRegion(start, end, binCount=binCount, indicator=get(seqnamesColName))]
	return(binnedDT)
}

getChromSizes = function(refAssembly) {
	# query available datasets
	ad = data(package="GenomicDistributions")
	adm = ad$results[,"Item"]
	chromSizesGenomeVar = paste0("chromSizes_", refAssembly)
	if (chromSizesGenomeVar %in% adm){
		# load it!
		data(list=chromSizesGenomeVar,
				package="GenomicDistributions",
				envir=environment())
		return(get(chromSizesGenomeVar))
	} else {
		message("I don't have archived chromSizes for reference assembly ",
			refAssembly)
	}

}

binGenome = function(genome, binCount) {
	chromSizes = getChromSizes(genome)
	return(binChroms(binCount, chromSizes))
}




#' Calculates the distribution of a query set over the genome
#' 
#' Returns a data.table showing counts of regions in GR, in the bins
#' In other words, where on which chromosomes are the ranges distributed?
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector that will be used to grab a BSGenome object
#'     by \code{binBSGenome}
#' @param binCount Number of bins to divide the chromosomes into
#' @param genomeBins You may supply pre-computed genome bins here (output from
#'     \code{binBSGenome}); if this is left empty, they will be computed
#' @export
genomicDistribution = function(query, refAssembly, binCount=10000, genomeBins=NULL) {
	if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, genomicDistribution, refAssembly, binCount, genomeBins)

		# To accomodate multiple regions, we'll need to introduce a new 'name'
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

	if (is.null(genomeBins)) {
		binnedDT = binGenome(refAssembly, binCount)
	}

	sdt = GenomicDistributions:::splitDataTable(binnedDT, "id")
	sdtl = lapply(sdt, GenomicDistributions:::dtToGr, chr="idCol")

	RDT = GenomicDistributions:::grToDt(query)
	# This jExpression will just count the number of regions.
	jExpr = ".N"
	res = GenomicDistributions:::BSAggregate(RDT, GRangesList(sdtl), jExpr=jExpr)

	# order chromosomes by current order.
	res[, chr:=factor(chr, levels=unique(res$chr))]
	return(res)
}

#' Plot distribution over chromosomes
#' 
#' Plots result from \code{genomicDistribution} calculation
#' @param GD The output from the genomicDistribution function
#' @param binCount Number of bins (should match the call to
#'     \code{genomicDistribution}
#' @param plotTitle Title for plot.
#' @export
plotGenomicDist = function(GD, binCount=10000, plotTitle="Distribution over chromosomes") {
	if ("name" %in% names(GD)){
		# It has multiple regions
		g = ggplot(GD, aes(x=withinGroupID, y=N, fill=name, color=name))
	} else {
		# It's a single region
		g = ggplot(GD, aes(x=withinGroupID, y=N))
	}
	g = g +
		xlab("Genome") + ylab("Number of regions") +
		geom_bar(stat="identity") + # Spread out to max width
		facet_grid(chr ~ .) + # Place chroms one on top of another
		theme_classic() + # Clean up cruft
		theme_blank_facet_label() + # No boxes around labels
		theme(panel.spacing=unit(0, "lines")) + # Reduce whitespace
		theme(strip.text.y=element_text(size=8, angle=0)) + # Rotate labels
		geom_hline(yintercept = 0, color="#EEEEEE") + # Light chrom lines
		scale_y_discrete(breaks=c(0, max(x$N)), limits=c(max(x$N))) + 
		scale_x_continuous(breaks=c(0,binCount), labels=c("Start", "End")) +
		theme(plot.title = element_text(hjust = 0.5)) + # Center title
		ggtitle(plotTitle) +
		theme(legend.position="bottom")
	return(g)
}
# genomicDistPlot(x)
