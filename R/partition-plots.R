
#' Calculates the distribution of overlaps for a query set to genomic partitions
#' 
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes with
#'     \code{getTSSs}.
#' @export
genomicPartitions = function(query, refAssembly) {
	geneModels = GenomicDistributions:::getGeneModels(refAssembly)
	partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
	return(assignPartitions(query, partitionList))
}


#' Create a basic genome partition list of genes, exons, introns, and intergenic
#' 
#' Given GRanges for genes, and a GRanges for exons, returns a list of GRanges
#' corresponding to various breakdown of the genome, based on the given
#' annotations; it gives you proximal and core promoters, exons, and introns. To
#' be used as a partionList for assignPartitions()
#' @param genesGR a GRanges object of gene coordinates
#' @param exonsGR a GRanges object of exons coordinates
#' @export
genomePartitionList = function(genesGR, exonsGR) {
	# Discard warnings (prompted from notifications to trim, which I do)
	withCallingHandlers({
		promCore = trim(promoters(genesGR, upstream=100, downstream=0))
		promProx = trim(promoters(genesGR, upstream=2000, downstream=0))
	}, warning=function(w) {
	    if (startsWith(conditionMessage(w), "GRanges object contains"))
	        invokeRestart("muffleWarning")
	})
	partitionList = list(promoterCore=promCore,
							promoterProx=promProx,
							exon=exonsGR,
							intron=genesGR)
	return(partitionList)
}


#' Takes a GRanges object, then assigns each element to a partition from the
#' provided partitionList, and then tallies the number of regions assigned to
#' each partition. A typical example of partitions is promoter, exon, intron,
#' etc; this function will yield the number of each for a query GRanges object
#' There will be a priority order to these, to account for regions that may
#' overlap multiple genomic partitions.
#' @param query 			 GRanges or GRangesList with regions to classify
#' @param partitionList     an ORDERED and NAMED list of genomic partitions
#'     GRanges. This list must be in priority order; the input will be assigned
#'     to the first partition it overlaps
#' @param remainder    A character vector to assign any query regions that do
#'     not overlap with anything in the partitionList. Defaults to "intergenic"
#' @export
assignPartitions = function(query, partitionList, remainder="intergenic") {
	if (methods::is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, assignPartitions, partitionList, remainder)
		nameList = names(query)
		if(is.null(nameList)) {
			nameList = 1:length(query) # Fallback to sequential numbers
		}
		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, sapply(x, nrow))
		return(xb)
	}

	#Overlap each of the partition list.
	partitionNames = names(partitionList)
	partition = rep(0, length(query));
	for (pi in 1:length(partitionList)) {
		cat(partitionNames[pi],":")
		ol = countOverlaps(query[partition==0], partitionList[[pi]])
		message("\tfound ", sum(ol>0))
		partition[partition==0][ol > 0] = partitionNames[pi]
	}
	partition[partition=="0"] = remainder
	tpartition = table(partition)
	
	return(data.frame(tpartition))
}

#' Produces a barplot showing how query regions of interest are distributed
#' across a given partition list
#' 
#' This function can be used to test a GRanges object against any arbitrary list
#' of genome partitions. The partition list is a priority-ordered list of
#' GRanges objects. Each region in the query will be assigned to a given
#' partition that it overlaps with the highest priority.
#' 
#' @param listGR a list of GRanges objects with independent region sets of interest
#' @param partitionList A priority-ordered and named list of GRanges objects to test
#'     query regions for overlaps
#'
#' @export
plotPartitions = function(assignedPartitions, labels=NULL) {
	# resAll = t(sapply(assignedPartitions, table))
	# resAllAve = sweep(resAll, 1, apply(resAll, 1, sum), FUN="/")*100
	# df = data.frame(partition=colnames(resAll), nOverlaps=t(resAll))

	if ("name" %in% names(assignedPartitions)){
		# It has multiple regions
		g = ggplot(assignedPartitions, aes(x=partition, y=Freq, fill=factor(name)))
	} else {
		g = ggplot(assignedPartitions, aes(x=partition, y=Freq))
	}

	g = g +
		geom_bar(stat="identity", position="dodge") + 
		theme_classic() + 
		theme_blank_facet_label() + # No boxes around labels
		theme(aspect.ratio=1) + 
		xlab("Genomic partition") +
		ylab("Number of overlaps") +
		theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
		theme(plot.title=element_text(hjust = 0.5)) + # Center title
		ggtitle(paste("Distribution across genomic partitions")) +
		scale_fill_discrete(name="User set") + 
		theme(legend.position="bottom")

	return(g)
}

partitionPercents = function(listGR, partitionList, backgroundGR = NULL) {
	if (! is(listGR, "list")) {
		# Try to correct for someone providing a single GR instead of a list.
		listGR = list(listGR)
	}
	res = lapply(listGR, assignPartitions, partitionList)
	classes = c(names(partitionList), 0)
	resAll = t(sapply(res, tableCount, classList=classes))

	resAllAve = sweep(resAll, 1, apply(resAll, 1, sum), FUN="/")*100
	rownames(resAllAve) = names(listGR)
	#incDecCol = c("goldenrod3", "navy", "purple", "orange")
	if (!is.null(backgroundGR)) {
		back = assignPartitions(backgroundGR, partitionList)
		backAll = table(back)
		backAllAve = sweep(backAll, 1, sum(backAll), FUN="/")*100
		resDiffAveNorm = log10(sweep(resAllAve, 2, backAllAve, FUN="/"))

		return(nlist(resAllAve, resDiffAveNorm))
	}
	return(nlist(resAllAve))
}

# A version for percentages, not yet activated
plotPartitionPercents = function(percList, labels = NULL) {
	if(is.null(labels)) {
		labels = rownames(percList$resAllAve)
	}
	colors = 1:NROW(percList$resAllAve)

	barplot(percList$resAllAve, beside=TRUE, col=colors, ylab="Percent")
	legend('topright', labels, pch=15, col=colors)
	if (! is.null(percList$resDiffAveNorm)) {
	barplot(percList$resDiffAveNorm, beside=TRUE, col=colors,
		ylab=expression('Log'[10]*'(fold change)'))
	legend('bottomright', labels, pch=15, col=colors)
	}
}



