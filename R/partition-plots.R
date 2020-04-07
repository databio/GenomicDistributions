#' Calculates the distribution of overlaps for a query set to a reference 
#' assembly
#'
#' This function is a wrapper for \code{calcPartitions} that uses built-in
#' partitions for a given reference genome assembly.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes 
#'     with \code{getTSSs}.
#' @return A data.frame indicating the number of query region overlaps in   
#'     several genomic partitions.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' calcPartitionsRef(query, "hg19")
calcPartitionsRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         refAssembly="character"))
	geneModels = getGeneModels(refAssembly)
	partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
	return(calcPartitions(query, partitionList))
}


#' Calculates the distribution of observed versus expected overlaps for a 
#' query set to a reference assembly
#'
#' This function is a wrapper for \code{calcExpectedPartitions} that uses 
#' built-in partitions for a given reference genome assembly.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes 
#'     with \code{getTSSs}.
#' @return A data.frame indicating the number of query region overlaps in   
#'     several genomic partitions.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' calcExpectedPartitionsRef(query, "hg19")
calcExpectedPartitionsRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         refAssembly="character"))
	geneModels = getGeneModels(refAssembly)
	partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
	return(calcExpectedPartitions(query, partitionList))
}


#' Calculates the cumulative distribution of overlaps for a query set to a
#' reference assembly
#' 
#' This function is a wrapper for \code{calcCumulativePartitions} that uses 
#' built-in partitions for a given reference genome assembly.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes 
#'     with \code{getTSSs}.
#' @return A data.frame indicating the number of query region overlaps in   
#'     several genomic partitions.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' calcCumulativePartitionsRef(query, "hg19")
calcCumulativePartitionsRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         refAssembly="character"))
	geneModels = getGeneModels(refAssembly)
	partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
	return(calcCumulativePartitions(query, partitionList))
}


#' Create a basic genome partition list of genes, exons, introns, and 
#' intergenic
#' 
#' Given GRanges for genes, and a GRanges for exons, returns a list of GRanges
#' corresponding to various breakdown of the genome, based on the given
#' annotations; it gives you proximal and core promoters, exons, and introns. 
#' To be used as a partionList for calcPartitions()
#' @param genesGR a GRanges object of gene coordinates
#' @param exonsGR a GRanges object of exons coordinates
#' @return A list of GRanges objects, each corresponding to a partition of the 
#' genome. Partitions include proximal and core promoters, exons and introns.
#' @export
#' @examples 
#' geneModels = getGeneModels("hg38")
#' partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
genomePartitionList = function(genesGR, exonsGR) {
    .validateInputs(list(exonsGR=c("GRanges", "GRangesList"), 
                         genesGR="GRanges"))
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


#' Calculates the distribution of overlaps between query and arbitrary genomic
#' partitions
#' 
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
#' @return A data.frame assigning each element of a GRanges object to a
#'  partition from a previously provided partitionList.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' geneModels = getGeneModels("hg38")
#' partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
#' calcPartitionsRef(query, partitionList)
calcPartitions = function(query, partitionList, remainder="intergenic") {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         partitionList="list"))
    if (methods::is(query, c("GRangesList"))) {
    	# Recurse over each GRanges object
    	x = lapply(query, calcPartitions, partitionList, remainder)
    	nameList = names(query)
    	if(is.null(nameList)) {
    		nameList = 1:length(query) # Fallback to sequential numbers
    	}
		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, vapply(x, nrow, integer(1)))
		return(xb)
	}
    #Overlap each of the partition list.
	partitionNames = names(partitionList)
	partition = rep(0, length(query))  
	for (pi in 1:length(partitionList)) {
		cat(partitionNames[pi],":")
		ol = suppressWarnings(
            countOverlaps(query[partition==0], partitionList[[pi]]))
		message("\tfound ", sum(ol>0))
		partition[partition==0][ol > 0] = partitionNames[pi]
	}
	partition[partition=="0"] = remainder
    tpartition = table(partition)
	return(data.frame(tpartition))
}


#' Calculates the distribution of overlaps between query and arbitrary genomic
#' partitions
#' 
#' Takes a GRanges object, then assigns each element to a partition from the
#' provided partitionList, and then tallies the number of regions assigned to
#' each partition. A typical example of partitions is promoter, exon, intron,
#' etc; this function will yield the number of each for a query GRanges object
#' There will be a priority order to these, to account for regions that may
#' overlap multiple genomic partitions.
#' @param query 		 GRanges or GRangesList with regions to classify.
#' @param partitionList  An ORDERED and NAMED list of genomic partitions
#'     GRanges. This list must be in priority order; the input will be assigned
#'     to the first partition it overlaps.
#' @return A data.frame assigning each element of a GRanges object to a
#'     partition from a previously provided partitionList.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' geneModels = getGeneModels("hg38")
#' partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
#' calcExpectedPartitions(query, partitionList)
calcExpectedPartitions = function(query, partitionList) {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         partitionList="list"))
    if (methods::is(query, c("GRangesList"))) {
    	# Recurse over each GRanges object
    	x = lapply(query, calcPartitions, partitionList, remainder)
    	nameList = names(query)
    	if(is.null(nameList)) {
    		nameList = 1:length(query) # Fallback to sequential numbers
    	}
		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, vapply(x, nrow, integer(1)))
		return(xb)
	}
    #Overlap each of the partition list.
	partitionNames = names(partitionList)
    partitionCounts = data.table::data.table(
        data.frame(N=elementNROWS(partitionList)), keep.rownames="name")
    partitionCounts = partitionCounts[order(partitionCounts$name)]
	partition = rep(0, length(query))
	for (pi in 1:length(partitionList)) {
		cat(partitionNames[pi],":")
		ol = suppressWarnings(
            countOverlaps(query[partition==0], partitionList[[pi]]))
		message("\tfound ", sum(ol>0))
		partition[partition==0][ol > 0] = partitionNames[pi]
	}
    # Remove remainder
	partition = partition[!partition=="0"]
    tpartition = table(partition)
	expectedPartitions = data.table::data.table(name=factor(names(tpartition)),
                                                observed=as.vector(tpartition))
    expectedPartitions[,expected:=partitionCounts$N]
    expectedPartitions[,log10OE:=log10(expectedPartitions$observed/
                                       expectedPartitions$expected)]
	return(expectedPartitions[match(partitionNames, expectedPartitions$name),])
}


#' Calculates the cumulative distribution of overlaps between query and 
#' arbitrary genomic partitions
#' 
#' Takes a GRanges object, then assigns each element to a partition from the
#' provided partitionList, and then tallies the number of regions assigned to
#' each partition. A typical example of partitions is promoter, exon, intron,
#' etc; this function will yield the number of each for a query GRanges object
#' There will be a priority order to these, to account for regions that may
#' overlap multiple genomic partitions.
#' @param query 		 GRanges or GRangesList with regions to classify.
#' @param partitionList  An ORDERED and NAMED list of genomic partitions
#'     GRanges. This list must be in priority order; the input will be assigned
#'     to the first partition it overlaps.
#' @return A data.frame assigning each element of a GRanges object to a
#'     partition from a previously provided partitionList.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' geneModels = getGeneModels("hg38")
#' partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
#' calcCumulativePartitions(query, partitionList)
calcCumulativePartitions = function(query, partitionList, remainder="intergenic") {
    .validateInputs(list(query=c("GRanges", "GRangesList"), 
                         partitionList="list"))
    if (methods::is(query, c("GRangesList"))) {
    	# Recurse over each GRanges object
    	x = lapply(query, calcPartitions, partitionList, remainder)
    	nameList = names(query)
    	if(is.null(nameList)) {
    		nameList = 1:length(query) # Fallback to sequential numbers
    	}
		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, vapply(x, nrow, integer(1)))
		return(xb)
	}

    #Overlap each of the partition list.
	partitionNames = names(partitionList)
    # Need total number of bases (not regions)
    query_total = sum(width(query))
    frif = data.table::data.table(name=as.character(),
                                  size=as.numeric(),
                                  count=as.numeric(),
                                  cumsum=as.numeric(),
                                  cumsize=as.numeric(),
                                  frif=as.numeric())
	for (pi in 1:length(partitionList)) {
		cat(partitionNames[pi],":")
        # Find overlaps
        hits  = suppressWarnings(findOverlaps(query, partitionList[[pi]]))
        olap  = suppressWarnings(pintersect(query[queryHits(hits)],
                                 partitionList[[pi]][subjectHits(hits)]))
        polap = width(olap) / width(partitionList[[pi]][subjectHits(hits)])
        hits  = data.table::data.table(xid=queryHits(hits),
                                       yid=subjectHits(hits),
                                       polap=polap)
        # Grab hits
        pHits = partitionList[[pi]][hits$yid]
        hits[, size:=width(pHits)]
        # Sum the weighted count column (polap*region size)
        hits[, count:= sum(polap*size), by=yid]
        message("\tfound ", nrow(hits))

        # Make mutually exclusive; remove hits from query
        query = query[-hits$xid]
        # Isolate hits
        hits = unique(hits, by="yid")
        # Link to positional data for feature of interest
        x = data.table::data.table(name=partitionNames[pi],
                                   size=hits$size,
                                   count=hits$count)
        x = x[order(x$count, x$size),]
        x$cumsum  = cumsum(x$count)
        x$cumsize = cumsum(x$size)
        x$frif    = x$cumsum/query_total
        frif = rbind(frif, x)
	}
    # Create remainder...
    cat(remainder,":")
    x = data.table::data.table(name=remainder, size=as.numeric(width(query)))
    message("\tfound ", length(query))
    x = x[order(x$size),]
    x$count   = x$size
    x$cumsum  = cumsum(x$count)
    x$cumsize = cumsum(x$size)
    x$frif    = x$cumsum/query_total
    return(rbind(frif, x))
}


#' Internal helper function for \code{plotCumulativePartitions}
#'
#' @param assignedPartitions Results from \code{calcCumulativePartitions}
setLabels = function(assignedPartitions) {
    if (is(assignedPartitions, "list")) {
		# It has multiple regions
		x = lapply(assignedPartitions, setLabels)
        x = data.table::rbindlist(x, idcol=TRUE)
		return(x)
	}
    return(data.table::data.table(
            xPos=0.95*max(log10(assignedPartitions$cumsize)),
            yPos=max(assignedPartitions$frif)+0.001,
            val=sprintf(max(assignedPartitions$frif),fmt="%#.2f"),
            color=NA,
            stringsAsFactors=FALSE))
}


#' Plot the cumulative distribution of regions in features
#'
#' This function plots the cumulative distribution of regions across a 
#' feature set.
#' @param assignedPartitions Results from \code{calcCumulativePartitions}
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the GenomicRanges or GenomicRangesList 
#'                      object.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' p = calcCumulativePartitionsRef(query, "hg19")
#' cumuPlot = plotCumulativePartitions(p)
plotCumulativePartitions = function(assignedPartitions, feature_names=NULL) {
    .validateInputs(list(assignedPartitions="data.frame"))
    if ("name" %in% names(assignedPartitions)){
        # It has multiple regions
        assignedPartitions = splitDataTable(assignedPartitions, "name")
    }
    plot_labels = setLabels(assignedPartitions)
    palette = colorRampPalette(c("#A6CEE3", "#025EBA", "#B2DF8A", "#05A602",
                                 "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                 "#CAB2D6", "#57069E", "#F0FC03", "#B15928"))

    plot_colors = palette(length(assignedPartitions))
   
    # Set colors
    plot_labels[,color:=plot_colors]
    
    if (is(assignedPartitions, "list")) {
        feature_lengths = data.table::data.table(num_feats=elementNROWS(assignedPartitions))
    } else {
        feature_lengths = data.table::data.table(num_feats=length(assignedPartitions))
    }
    
    if (is(assignedPartitions, "list")) {
        assignedPartitions = data.table::rbindlist(assignedPartitions,
                                                   idcol=TRUE)
    }

    # If name vector provided, update names
    if (all(!is.null(feature_names))) {
        if (length(feature_names) == nrow(feature_lengths)) {
            plot_labels[,.id:=feature_names]
            assignedPartitions[,.id:=rep(feature_names,
                               each=feature_lengths$num_feats)]
        } else {
            if (!".id" %in% colnames(plot_labels)) {
                plot_labels[,.id:=seq(1:nrow(feature_lengths))]
                assignedPartitions[,.id:=rep(seq(1:nrow(feature_lengths)),
                                   feature_lengths$num_feats)]
            }
        }
    } else {
        if (!".id" %in% colnames(plot_labels)) {
            plot_labels[,.id:=seq(1:nrow(feature_lengths))]
            assignedPartitions[,.id:=rep(seq(1:nrow(feature_lengths)),
                               feature_lengths$num_feats)]
        }
    }

    p = ggplot(assignedPartitions,
               aes(x=log10(cumsize), y=frif, group=.id, color=.id))
    p = p + 
        geom_line(size=2, alpha=0.5) +
        guides(linetype = FALSE) +
        labs(x=expression(log[10]("number of bases")),
             y="Cumulative distribution across genomic partitions") +
        theme_classic() +
        theme(axis.line = element_line(size = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank(),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent",
                                             color = NA),
              legend.background = element_rect(fill = "transparent",
                                               color = NA),
              legend.box.background = element_rect(fill = "transparent",
                                                   color = NA),
              aspect.ratio = 1,
              legend.position = "none",
              plot.title = element_text(hjust = 0.5),              
              panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )

    # Recolor and reposition legend
    p = p + scale_color_manual(
            labels=paste0(plot_labels$.id, ": ", plot_labels$val),
            values=plot_labels$color) +
        labs(color="Cumulative distribution across genomic partitions") +
        theme(legend.position=c(0.075,0.975),
              legend.justification=c(0.1,0.9),
              legend.title = element_blank(),
              legend.text = element_text(size = rel(0.65)), 
              legend.background=element_blank(),
              legend.key = element_blank(),
              axis.text.x = element_text(angle = 0, hjust = 1, vjust=0.5)
        )

    if (!exists("p")) {
        p = ggplot()
    }

    return(p)
}


#' Produces a barplot showing how query regions of interest are distributed
#' relative to the expected distribution across a given partition list
#'
#' @param expectedPartitions  A data.frame holding the frequency of assignment 
#'     to each of the partitions, the expected number of each partition, and
#'     the log10 of the observed over expected. Produced by 
#'     \code{calcExpectedPartitions}.
#' @param feature_names  Character vector with labels for the partitions 
#'     (optional). By default it will use the names from the first argument.
#' @return A ggplot object using a barplot to show the distribution of the 
#'     query regions across a given partition list.  
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' p = calcExpectedPartitionsRef(query, "hg19")
#' expectedPlot = plotExpectedPartitions(p)
plotExpectedPartitions = function(expectedPartitions, feature_names=NULL) {
    .validateInputs(list(expectedPartitions="data.frame"))
    palette = colorRampPalette(c("#A6CEE3", "#025EBA", "#B2DF8A", "#05A602",
                                 "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                 "#CAB2D6", "#57069E", "#F0FC03", "#B15928"))
    # Add 1 for the remainder
    plot_colors = palette(nrow(expectedPartitions)+1)
    
    # Set colors (drop the remainder)
    expectedPartitions[,color:=plot_colors[1:nrow(expectedPartitions)]]

    # If name vector provided, update names
    if (all(!is.null(feature_names))) {
        if (length(feature_names) == nrow(expectedPartitions)) {
            expectedPartitions[,name:=feature_names]
        } else {
            if (!"name" %in% colnames(plot_labels)) {
                expectedPartitions[,name:=seq(1:nrow(expectedPartitions))]
            }
        }
    }
    
    expectedPartitions = expectedPartitions[order(expectedPartitions$log10OE),]
    expectedPartitions$name  = factor(expectedPartitions$name,
                                      levels=expectedPartitions$name)
    expectedPartitions$color = factor(expectedPartitions$color,
                                      levels=expectedPartitions$color)

    p = ggplot(expectedPartitions, aes(x = name, y = log10OE))
    p = p + 
        geom_bar(stat="identity",
                 fill = expectedPartitions$color,
                 alpha = 0.5) + 
        geom_hline(aes(yintercept=0), linetype="dotted") +
        xlab('') +
        ylab(expression(log[10](over(Obs, Exp)))) +
        coord_flip() +
        theme(axis.line = element_line(size = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent",
                                             color = NA),
              legend.background = element_rect(fill = "transparent",
                                               color = NA),
              legend.box.background = element_rect(fill = "transparent",
                                                   color = NA),
              aspect.ratio = 1,
              legend.position = "none",
              plot.title = element_text(hjust = 0.5),
              panel.border = element_rect(colour = "black", fill=NA,
                                          size=0.5)
        )

    if (!exists("p")) {
        p = ggplot()
    }

    return(p)
}


#' Produces a barplot showing how query regions of interest are distributed
#' across a given partition list
#' 
#' This function can be used to test a GRanges object against any arbitrary 
#' list of genome partitions. The partition list is a priority-ordered list of
#' GRanges objects. Each region in the query will be assigned to a given
#' partition that it overlaps with the highest priority.
#' 
#' @param assignedPartitions  A table holding the frequency of assignment to
#'     each of the partitions. Produced by \code{calcPartitions}
#' @param labels Character vector with labels for the partitions (optional). By
#'     default it will use the names from the first argument.
#' @return A ggplot object using a barplot to show the distribution of the query 
#'  regions across a given partition list.  
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' p = calcPartitionsRef(query, "hg19")
#' partPlot = plotPartitions(p)
plotPartitions = function(assignedPartitions, labels=NULL) {
	# resAll = t(sapply(assignedPartitions, table))
	# resAllAve = sweep(resAll, 1, apply(resAll, 1, sum), FUN="/")*100
	# df = data.frame(partition=colnames(resAll), nOverlaps=t(resAll))
    .validateInputs(list(assignedPartitions="data.frame"))
	if ("name" %in% names(assignedPartitions)){
		# It has multiple regions
		g = ggplot(assignedPartitions, 
		           aes(x=partition, y=Freq, fill=factor(name)))
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
		theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + 
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
	res = lapply(listGR, calcPartitions, partitionList)
	classes = c(names(partitionList), 0)
	resAll = t(sapply(res, tableCount, classList=classes))

	resAllAve = sweep(resAll, 1, apply(resAll, 1, sum), FUN="/")*100
	rownames(resAllAve) = names(listGR)
	#incDecCol = c("goldenrod3", "navy", "purple", "orange")
	if (!is.null(backgroundGR)) {
		back = calcPartitions(backgroundGR, partitionList)
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
