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
#' 
#' To be used as a partitionList for \code{calcPartitions}.
#' 
#' @param genesGR a GRanges object of gene coordinates
#' @param exonsGR a GRanges object of exons coordinates
#' @return A list of GRanges objects, each corresponding to a partition of the 
#'     genome. Partitions include proximal and core promoters, exons and 
#'     introns. 
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
#'
#' @param query GRanges or GRangesList with regions to classify
#' @param partitionList an ORDERED and NAMED list of genomic partitions
#'     GRanges. This list must be in priority order; the input will be assigned
#'     to the first partition it overlaps
#' @param remainder A character vector to assign any query regions that do
#'     not overlap with anything in the partitionList. Defaults to "intergenic"
#' @return A data.frame assigning each element of a GRanges object to a
#'     partition from a previously provided partitionList.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' geneModels = getGeneModels("hg38")
#' partitionList = genomePartitionList(geneModels$genesGR, geneModels$exonsGR)
#' calcPartitions(query, partitionList)
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
#'
#' @param query GRanges or GRangesList with regions to classify.
#' @param partitionList An ORDERED and NAMED list of genomic partitions
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
        x = lapply(query, calcExpectedPartitions, partitionList)
        nameList = names(query)
        if(is.null(nameList)) {
            nameList = 1:length(query) # Fallback to sequential numbers
        }
        # Append names
        xb = data.table::rbindlist(x)
        xb$name = rep(nameList, vapply(x, nrow, integer(1)))
        return(xb)
    }
    #Overlap each of the partition list.
    partitionNames = names(partitionList)
    #expected scaled by number of regions in query
    query_total = length(query)
    partitionCounts = data.table::data.table(
        data.frame(N=(elementNROWS(partitionList)/
                      sum(elementNROWS(partitionList))*query_total)),
        keep.rownames="partition")
    partitionCounts = partitionCounts[order(partitionCounts$partition)]
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
    expectedPartitions = data.table::data.table(
        partition=factor(names(tpartition)), observed=as.vector(tpartition))
    expectedPartitions = merge(expectedPartitions, partitionCounts,
                               by = "partition")
    setnames(expectedPartitions,"N","expected")
    expectedPartitions[,log10OE:=log10(expectedPartitions$observed/
                                       expectedPartitions$expected)]
    return(expectedPartitions[match(partitionNames,
           expectedPartitions$partition),])
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
#' 
#' @param query GRanges or GRangesList with regions to classify.
#' @param partitionList An ORDERED and NAMED list of genomic partitions
#'     GRanges. This list must be in priority order; the input will be assigned
#'     to the first partition it overlaps.
#' @param remainder  Which partition do you want to account for 'everything 
#'     else'?
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
        x = lapply(query, calcCumulativePartitions, partitionList, remainder)
        nameList = names(query)
        if(is.null(nameList)) {
            nameList = 1:length(query) # Fallback to sequential numbers
        }
        # Append names
        xb = data.table::rbindlist(x)
        xb$name = rep(nameList, vapply(x, nrow, integer(1)))
        return(xb)
    }

    #Overlap each of the partition list.
    partitionNames = names(partitionList)
    # Need total number of bases (not regions)
    query_total = sum(width(query))
    frif = data.table::data.table(partition=as.character(),
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
        x = data.table::data.table(partition=partitionNames[pi],
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
    x = data.table::data.table(partition=remainder,
                               size=as.numeric(width(query)))
    message("\tfound ", length(query))
    x = x[order(x$size),]
    x$count   = x$size
    x$cumsum  = cumsum(x$count)
    x$cumsize = cumsum(x$size)
    x$frif    = x$cumsum/query_total
    return(rbind(frif, x))
}


# Internal helper function for \code{plotCumulativePartitions}.
#
# @param assignedPartitions Results from \code{calcCumulativePartitions}.
# @return A data.table object of partition names and values.
setLabels = function(assignedPartitions) {
    if (methods::is(assignedPartitions, c("list"))){
        # It has multiple regions
        x = lapply(assignedPartitions, setLabels)
        nameList = names(assignedPartitions)
        if(is.null(nameList)) {
            nameList = 1:length(assignedPartitions) # Fallback to sequential numbers
        }
        # Append names
        xb = data.table::rbindlist(x)
        xb$name = rep(nameList, vapply(x, nrow, integer(1)))
        return(xb)
    }
    partition = assignedPartitions[, partition, by=partition]
    xPos = assignedPartitions[, 0.95*max(log10(cumsize)), by=partition]
    yPos = assignedPartitions[, max(frif)+0.001, by=partition]
    val  = assignedPartitions[, sprintf(max(frif),fmt="%#.2f"), by=partition]
    return(data.table::data.table(partition=partition$partition,
                                  xPos=xPos$V1,
                                  yPos=yPos$V1,
                                  val=val$V1,
                                  stringsAsFactors=FALSE))
}


#' Plot the cumulative distribution of regions in features
#'
#' This function plots the cumulative distribution of regions across a 
#' feature set.
#'
#' @param assignedPartitions Results from \code{calcCumulativePartitions}
#' @param feature_names An optional character vector of feature names, in the 
#'     same order as the GenomicRanges or GenomicRangesList object.
#' @return A ggplot object of the cumulative distribution of regions in 
#'     features.
#' @export
#' @examples 
#' f = system.file("extdata", "vistaEnhancers.bed.gz",
#'     package="GenomicDistributions")
#' query = rtracklayer::import(f)
#' p = calcCumulativePartitionsRef(query, "hg19")
#' cumuPlot = plotCumulativePartitions(p)
plotCumulativePartitions = function(assignedPartitions, feature_names=NULL) {
    .validateInputs(list(assignedPartitions="data.frame"))
    palette = colorRampPalette(c("#A6CEE3", "#025EBA", "#B2DF8A", "#05A602",
                                 "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
                                 "#CAB2D6", "#57069E", "#F0FC03", "#B15928"))
    if ("name" %in% names(assignedPartitions)){
        # It has multiple regions
        p = ggplot(assignedPartitions, aes(x=log10(cumsize), y=frif,
                   group=partition, color=partition)) +
            facet_grid(. ~name)
        plot_labels = setLabels(splitDataTable(assignedPartitions, "name"))
        partition_sizes = assignedPartitions[, .N, by=.(partition, name)]
        plot_labels[, label:=sprintf(" %s:%s", plot_labels$partition,
                                     plot_labels$val)]
        label = plot_labels[, list(label = paste(label, collapse="\n"))
                            , by = name]
    } else {
        p = ggplot(assignedPartitions,
               aes(x=log10(cumsize), y=frif, group=partition, color=partition))
        plot_labels = setLabels(assignedPartitions)
        partition_sizes = assignedPartitions[, .N, by=partition]
        plot_labels[, label:=sprintf(" %s:%s", plot_labels$partition,
                                     plot_labels$val)]
        label = plot_labels[, list(label = paste(label, collapse="\n"))]
    }

    # If name vector provided, update names
    if (all(!is.null(feature_names))) {
        if (length(feature_names) == nrow(partition_sizes)) {
            plot_labels[,partition:=feature_names]
            assignedPartitions[,partition:=rep(feature_names,
                               each=partition_sizes$num_feats)]
        } else {
            if (!"partition" %in% colnames(plot_labels)) {
                plot_labels[,partition:=seq(1:nrow(partition_sizes))]
                assignedPartitions[,
                    partition:=rep(seq(1:nrow(partition_sizes)),
                                   partition_sizes$num_feats)]
            }
        }
    } else {
        if (!"partition" %in% colnames(plot_labels)) {
            plot_labels[,partition:=seq(1:nrow(partition_sizes))]
            assignedPartitions[,partition:=rep(seq(1:nrow(partition_sizes)),
                               partition_sizes$num_feats)]
        }
    }

    p = p + 
        geom_line(size=2, alpha=0.5) +
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
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5),              
              panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )

    # Add label text
    p = p +
        geom_text(data=label, mapping=aes(x=-Inf, y=Inf, label=label),
        hjust="inward", vjust=1.05, inherit.aes=FALSE)

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
    expectedPartitions = na.omit(expectedPartitions)
    if ("name" %in% names(expectedPartitions)){
        # It has multiple regions
        p = ggplot(expectedPartitions, 
                   aes(x=partition, y=log10OE, fill=factor(name)))
    } else {
        p = ggplot(expectedPartitions, aes(x=partition, y=log10OE))
    }

    # If feature name vector provided, update partition names
    if (all(!is.null(feature_names))) {
        if (length(feature_names) == nrow(expectedPartitions)) {
            expectedPartitions[,partition:=feature_names]
        } else {
            if (!"partition" %in% colnames(plot_labels)) {
                expectedPartitions[,partition:=seq(1:nrow(expectedPartitions))]
            }
        }
    }

    if ("name" %in% names(expectedPartitions)){
        expectedPartitions = expectedPartitions[
            order(expectedPartitions$name, expectedPartitions$log10OE),]
    } else {
        expectedPartitions = expectedPartitions[
            order(expectedPartitions$log10OE),]
    }
    
    p = p + 
        geom_bar(stat="identity", position="dodge") + 
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


# Calculate the percentage overlap of a list of GRanges object against a 
# list of known partitions.
#
# @param listGR  A list of GRanges objects.
# @param partitionList A character vector of partition names.
# @param backgroundGR A GRanges object to remove from listGR overlap 
#     as background.
# @return A named list of the frequency of input regions against provided
#     partitions.
calcPartitionPercents = function(listGR, partitionList, backgroundGR=NULL) {
    .validateInputs(list(listGR=c("GRanges","GRangesList"),
                         partitionList="list"))
    if(!(is(listGR, "GRangesList"))){
        listGR = GRangesList(listGR)
    }
    res = lapply(listGR, calcPartitions, partitionList)
    resAll = Reduce(function(x,y) merge(x, y, by="partition"), res)
    names(resAll) = c("Partition", names(res))
    rownames(resAll) = resAll$Partition
    resAll[, 1] = NULL
    # should transpose df to calculate percentages
    tResAll = t(resAll)
    resAllAve = sweep(tResAll, 1, apply(tResAll, 1, sum), FUN="/")*100
    resAllAve = as.data.frame(resAllAve)
    # if background GRanges is provided
    if (!is.null(backgroundGR)) {
        back = calcPartitions(backgroundGR, partitionList)
        rownames(back) = back$Partition
        back[, 1] = NULL
        backAll = t(back)
        backAllAve = sweep(backAll, 1, sum(backAll), FUN="/")*100
        resDiffAveNorm = log10(sweep(resAllAve, 2, backAllAve, FUN="/"))
        return(nlist(resAllAve, resDiffAveNorm))
    }
    return(nlist(resAllAve))
}


# A version for percentages, not yet activated
# Plot the percentage overlap of a list of GRanges objects.
# 
# @param percList A named list of percentage overlap of regions to genomic
#     partitions.
plotPartitionPercents = function(percentList, labels = NULL) {
    .validateInputs(list(percentList="list"))
    if(is.null(labels)) {
        labels = rownames(percentList$resAllAve)
    }
    # need to reshape the data to account for diff regionsets and partitions
    percData = as.matrix(t(percentList$resAllAve))
    percReshaped = reshape2::melt(percData, value.name="Percent")
    colnames(percReshaped)[colnames(percReshaped) == "Var2"] = "regionSet"
    g = ggplot2::ggplot(percReshaped, aes(x=Var1, y=Percent, fill=regionSet)) +
      geom_bar(stat="identity", position = position_dodge()) +
      xlab("Genomic Partition") +
      theme_classic() +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) +
      theme(plot.title=element_text(hjust = 0.5)) +
      theme(aspect.ratio=1) +
      ggtitle("Percentage distribution across genomic partitions") 
  
    # If multiple regionsets are provided
    if (length(rownames(percentList$resAllAve)) > 1) { 
        g = g + theme(legend.position = "bottom") 
    } else {
        # If a single regionset provided, no need to include legend
        g = g + theme(legend.position = "none")
    }
    if (!is.null(percentList$resDiffAveNorm)) {
        b = barplot(percentList$resDiffAveNorm, beside=TRUE, col=colors,
              ylab=expression('Log'[10]*'(fold change)'))
        return(b)
    } else {
        return(g)
    }
}

# # A version for percentages, not yet activated
# plotPartitionPercents = function(percList, labels = NULL) {
#     if(is.null(labels)) {
#         labels = rownames(percList$resAllAve)
#     }
#     colors = 1:NROW(percList$resAllAve)
# 
#     barplot(percList$resAllAve, beside=TRUE, col=colors, ylab="Percent")
#     legend('topright', labels, pch=15, col=colors)
#     if (! is.null(percList$resDiffAveNorm)) {
#     barplot(percList$resDiffAveNorm, beside=TRUE, col=colors,
#         ylab=expression('Log'[10]*'(fold change)'))
#     legend('bottomright', labels, pch=15, col=colors)
#     }
# }
