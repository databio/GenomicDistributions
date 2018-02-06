
#' Aggregating signals in bins across a set of regions
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
#' @param bins The number of bins to divide this segment
#' @param idDF An identifier vector to keep with your bins, in case you are
#'     doing this on a long table with multiple segments concatenated
#'
#' @return
#' A data.table, expanded to nrow = number of bins, with these id columns:
#' 		id: region ID
#' 		binID: repeating ID (this is the value to aggregate across)
#' 		ubinID: unique bin IDs
#' @export
#' @examples
#' loadCGData("hg19")
#' cgIslandsDT = data.table(...)
#' binnedCGI = cgIslandsDT[, binRegion(start, end, 50)]
binRegion = function(start, end, bins, idDF=NULL) {
	#if (!is.null(idDF) & ( ! "data.frame"  %in% class(idDF))) {
	#	stop("idDF should be a data.frame")
	#}
	binSize = (end-start)/(bins)
	breaks = round(rep(start, each=(bins+1)) + (0:(bins)) * rep(binSize, each=(bins+1)))

	endpoints = (bins+1) * (1:(length(start)))
	startpoints = 1 + (bins+1)  * (0:(length(start)-1))
	#TODO: remove this code split
	if (is.null(idDF)) {
		dt = data.table(start=breaks[-endpoints], 
						end=breaks[-startpoints],
						id=rep((1:length(start)), each=bins),
						binID= 1:bins,
						ubinID=1:length(breaks[-startpoints]),
						key="id")
	} else {
		chr = rep(idDF, each=bins)
		dt = data.table(chr, 
						start=breaks[-endpoints],
						end=breaks[-startpoints],
						id=rep((1:length(start)), each=bins),
						binID= 1:bins,
						ubinID=1:length(breaks[-startpoints]),
						key="id")

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
	rangeDT = data.table(chr=seqnames(BSG), start=1, end=seqlengths(BSG))
	seqnamesColName="chr"
	binnedDT = rangeDT[, binRegion(start, end, binCount, get(seqnamesColName))]
	return(binnedDT)
}


#' Calculates the distribution of a query set over the genome
#' 
#' Returns a data.table showing counts of regions in GR, in the bins
#' In other words, where on which chromosomes are the ranges distributed?
#' @param genome A character vector that will be used to grab a BSGenome object
#'     by binBSGenome
#' @export
genomicDistribution = function(GR, genome, binCount=1000, genomeBins=NULL) {
	if (is.null(genomeBins)) {
		binnedDT = binBSGenome(genome, binCount)
	}

	sdt = splitDataTable(binnedDT, "id")
	sdtl = lapply(sdt, dtToGr)
#	GRangesList(sdtl)

	RDT = grToDt(GR)
	jExpr = ".N"
	res = BSAggregate(RDT, GRangesList(sdtl), jExpr=jExpr)

	# order chromosomes by curent order.
	res[, chr:=factor(chr, levels=unique(res$chr))]
	return(res)
}

#' Plot distribution over chromosomes
#' 
#' Plots result from \code{genomicDistribution} calculation
#' @param GD The output from the genomicDistribution function
#' @param binCount Number of bins (should match the call to
#'     \code{genomicDistribution}
#' @export
plotGenomicDist = function(GD, binCount=1000,
	plotTitle="Distribution over chromosomes") {
	g = ggplot(GD, aes(x=withinGroupID, y=N)) + 
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
		ggtitle(plotTitle)
	return(g)
}
# genomicDistPlot(x)
