#' Calculates the distribution of overlaps for a query set to a reference
#' assembly
#'
#' This function is a wrapper for \code{calcPartitions}
#' and \code{calcPartitionPercents} that uses built-in
#' partitions for a given reference genome assembly.
#'
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab annotation
#'     models with \code{getGeneModels}
#' @param bpProportion logical indicating if overlaps should be calculated
#'     based on number of base pairs overlapping with each partition.
#'     bpProportion=FALSE does overlaps in priority order,
#'     bpProportion=TRUE counts number of overlapping
#'     base pairs between query and each partition.
#' @return A data.frame indicating the number of query region overlaps in
#'     several genomic partitions.
#' @export
#' @examples
#' calcPartitionsRef(vistaEnhancers, "hg19")
#' calcPartitionsRef(vistaEnhancers, "hg19", bpProportion=TRUE)
calcPartitionsRef = function(query, refAssembly, bpProportion=FALSE){
    .validateInputs(list(query=c("GRanges", "GRangesList"),
                         refAssembly="character"))
    geneModels = getGeneModels(refAssembly)
    partitionList = genomePartitionList(geneModels$genesGR,
                                        geneModels$exonsGR,
                                        geneModels$threeUTRGR,
                                        geneModels$fiveUTRGR)
    message("Calculating overlaps...")
    return(calcPartitions(query, partitionList, bpProportion=bpProportion))
}


#' Calculates the distribution of observed versus expected overlaps for a
#' query set to a reference assembly
#'
#' This function is a wrapper for \code{calcExpectedPartitions} that uses
#' built-in partitions for a given reference genome assembly.
#'
#' @param query A GenomicRanges or GenomicRangesList object with query regions
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab annotation
#'     models with \code{getGeneModels}, and chromosome sizes with\code{getChromSizes}
#' @param bpProportion logical indicating if overlaps should be calculated based
#'     on number of base pairs overlapping with each partition.
#'     bpProportion=FALSE does overlaps in priority order,
#'     bpProportion=TRUE counts number of overlapping
#'     base pairs between query and each partition.
#' @return A data.frame indicating the number of query region overlaps in
#'     several genomic partitions.
#' @export
#' @examples
#' \dontrun{
#' calcExpectedPartitionsRef(vistaEnhancers, "hg19")
#' calcExpectedPartitionsRef(vistaEnhancers, "hg19", bpProportion=TRUE)
#' }
calcExpectedPartitionsRef = function(query, refAssembly,
                                     bpProportion=FALSE) {
    .validateInputs(list(query=c("GRanges", "GRangesList"),
                         refAssembly="character"))
    geneModels = getGeneModels(refAssembly)
    chromSizes = getChromSizes(refAssembly)
    genomeSize = sum(chromSizes)
    partitionList = genomePartitionList(geneModels$genesGR,
                                        geneModels$exonsGR,
                                        geneModels$threeUTRGR,
                                        geneModels$fiveUTRGR)
    expectedPartitions = calcExpectedPartitions(query, partitionList,
                           genomeSize, bpProportion=bpProportion)
    return(expectedPartitions)
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
#' calcCumulativePartitionsRef(vistaEnhancers, "hg19")
calcCumulativePartitionsRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges", "GRangesList"),
                         refAssembly="character"))
    geneModels = getGeneModels(refAssembly)
    partitionList = genomePartitionList(geneModels$genesGR,
                                        geneModels$exonsGR,
                                        geneModels$threeUTRGR,
                                        geneModels$fiveUTRGR)
    return(calcCumulativePartitions(query, partitionList))
}


#' Create a basic genome partition list of genes, exons, introns, UTRs, and
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
#' @param threeUTRGR a GRanges object of 3' UTRs
#' @param fiveUTRGR a GRanges object of 5' UTRs
#' @return A list of GRanges objects, each corresponding to a partition of the
#'     genome. Partitions include proximal and core promoters, exons and
#'     introns.
#' @export
#' @examples
#' partitionList = genomePartitionList(geneModels_hg19$genesGR,
#'                                     geneModels_hg19$exonsGR,
#'                                     geneModels_hg19$threeUTRGR,
#'                                     geneModels_hg19$fiveUTRGR)
genomePartitionList = function(genesGR, exonsGR, threeUTRGR=NULL,
                               fiveUTRGR=NULL) {
    .validateInputs(list(exonsGR=c("GRanges", "GRangesList"),
                         genesGR="GRanges",
                         threeUTRGR=c("GRanges", "GRangesList", "NULL"),
                         fiveUTRGR=c("GRanges", "GRangesList", "NULL")))
    # Discard warnings (prompted from notifications to trim, which I do)
    withCallingHandlers({
        promCore = trim(promoters(genesGR, upstream=100, downstream=0))
        promProx = trim(promoters(genesGR, upstream=2000, downstream=0))
    }, warning=function(w) {
        if (startsWith(conditionMessage(w), "GRanges object contains"))
            invokeRestart("muffleWarning")
    })

    # subtract overlaps (promoterCore lies within PromoterProx)
    promoterProx = GenomicRanges::setdiff(promProx, promCore)

    # remove any possible overlaps between classes
    fiveUTRGR = GenomicRanges::setdiff(fiveUTRGR, threeUTRGR)
    exonsGR = GenomicRanges::setdiff(exonsGR, threeUTRGR)
    exonsGR = GenomicRanges::setdiff(exonsGR, fiveUTRGR)

    #   introns = gene - (5'UTR, 3'UTR, exons)
    nonThree = GenomicRanges::setdiff(genesGR, threeUTRGR)
    nonThreeFive = GenomicRanges::setdiff(nonThree, fiveUTRGR)
    intronGR = GenomicRanges::setdiff(nonThreeFive, exonsGR)

    partitionList = list(promoterCore=promCore,
                         promoterProx=promoterProx,
                         threeUTR=threeUTRGR,
                         fiveUTR=fiveUTRGR,
                         exon=exonsGR,
                         intron=intronGR)


    return(Filter(Negate(is.null), partitionList))
}


#' Calculates the distribution of overlaps between
#' query and arbitrary genomic partitions
#'
#' Takes a GRanges object, then assigns each element to a partition from the
#' provided partitionList, and then tallies the number of regions assigned to
#' each partition. A typical example of partitions is promoter, exon, intron,
#' etc; this function will yield the number of each for a query GRanges object
#' There will be a priority order to these, to account for regions that may
#' overlap multiple genomic partitions.
#'
#' @param query GRanges or GRangesList with regions to classify
#' @param partitionList an ORDERED (if bpProportion=FALSE) and NAMED list of
#'     genomic partitions GRanges. This list must be in priority order; the
#'     input will be assigned to the first partition it overlaps.
#'     bpProportion=TRUE, the list does not need ordering.
#' @param remainder A character vector to assign any query regions that do
#'     not overlap with anything in the partitionList. Defaults to "intergenic"
#' @param bpProportion logical indicating if overlaps should be calculated based
#'     on number of base pairs overlapping with each partition.
#'     bpProportion=FALSE does overlaps in priority order,
#'     bpProportion=TRUE counts number of overlapping
#'     base pairs between query and each partition.
#' @return A data.frame assigning each element of a GRanges object to a
#'     partition from a previously provided partitionList.
#' @export
#' @examples
#' partitionList = genomePartitionList(geneModels_hg19$genesGR,
#'                                     geneModels_hg19$exonsGR,
#'                                     geneModels_hg19$threeUTRGR,
#'                                     geneModels_hg19$fiveUTRGR)
#' calcPartitions(vistaEnhancers, partitionList)
#' calcPartitions(vistaEnhancers, partitionList, bpProportion=TRUE)
calcPartitions = function(query, partitionList,
                          remainder="intergenic", bpProportion=FALSE) {
  .validateInputs(list(query=c("GRanges", "GRangesList"),
                       partitionList="list"))
  if (methods::is(query, c("GRangesList"))) {
    # Recurse over each GRanges object
    x = lapply(query, calcPartitions, partitionList, remainder, bpProportion)
    nameList = names(query)
    if(is.null(nameList)) {
      newnames = seq_along(query) # Fallback to sequential numbers
      nameList = names
    }
    # Append names
    xb = rbindlist(x)
    xb$name = rep(nameList, vapply(x, nrow, integer(1)))
    return(xb)
  }

  # proportional partitions
  if (bpProportion){
    # do overlap with each partition and record the overlap widths
    totalOverlap = lapply(partitionList, overlapWidths, query)

    # calculate the number of bases that did not fall anywhere - remainder
    remainderBases = sum(width(query)) - sum(unlist(totalOverlap))

    # there are some remainder overlaps between classes (e.g. promoter and 5'UTR)
    # their elimination might not make functional sense, but the double class can
    # lead to negative numbers in intergenic regions (always a very small number)
    # to eliminate this - set negative numbers to 0
    if (remainderBases < 0){
      remainderBases = 0
    }

    # gather all overlaps into data.frame
    propPartitions = data.frame(partition = c(names(totalOverlap),
                                              remainder),
                                bpOverlap = c(unlist(totalOverlap),
                                              remainderBases))
    propPartitions$frequency = propPartitions$bpOverlap / sum(propPartitions$bpOverlap)
    return(propPartitions)
  } else {
    # priority overlap partitions
    #Overlap each of the partition list.
    partitionNames = names(partitionList)
    partition = rep(0, length(query))
    for (pi in seq_along(partitionList)) {
      #message(partitionNames[pi],":")
      ol = suppressWarnings(
        countOverlaps(query[partition==0], partitionList[[pi]]))
      #message("\tfound ", sum(ol>0))
      partition[partition==0][ol > 0] = partitionNames[pi]
    }
    partition[partition=="0"] = remainder
    tpartition = table(partition)
    tpartition = data.frame(tpartition)

    partitionNamesFull = c(partitionNames, remainder)
    if (!all(partitionNamesFull %in% tpartition$partition)){
      notIncluded = partitionNamesFull[!(partitionNamesFull %in%
                                           tpartition$partition)]
      addRows = data.frame(partition = notIncluded,
                           Freq = rep(0, length(notIncluded)))
      tpartition = rbind(tpartition, addRows)
    }

    return(data.frame(tpartition))
  }
}

#' Calculates expected partiton overlap based on contribution of each
#' feature (partition) to genome size. Expected and observed overlaps
#' are then compared.
#'
#' @param query GRanges or GRangesList with regions to classify.
#' @param partitionList An ORDERED (if bpProportion=FALSE) and NAMED
#'     list of genomic partitions GRanges. This list must be in
#'     priority order; the input will be assigned
#'     to the first partition it overlaps. However, if bpProportion=TRUE,
#'     the list does not need ordering.
#' @param genomeSize The number of bases in the query genome. In other words,
#'     the sum of all chromosome sizes.
#' @return A data.frame assigning each element of a GRanges object to a
#'     partition from a previously provided partitionList.
#' @param remainder  Which partition do you want to account for 'everything
#'     else'?
#' @param bpProportion logical indicating if overlaps should be calculated based
#'     on number of base pairs overlapping with each partition.
#'     bpProportion=FALSE does overlaps in priority order,
#'     bpProportion=TRUE counts number of overlapping
#'     base pairs between query and each partition.
#' @export
#' @examples
#' partitionList = genomePartitionList(geneModels_hg19$genesGR,
#'                                     geneModels_hg19$exonsGR,
#'                                     geneModels_hg19$threeUTRGR,
#'                                     geneModels_hg19$fiveUTRGR)
#' chromSizes = getChromSizes('hg19')
#' genomeSize = sum(chromSizes)
#' calcExpectedPartitions(vistaEnhancers, partitionList, genomeSize)
#' calcExpectedPartitions(vistaEnhancers, partitionList, genomeSize, bpProportion=TRUE)
calcExpectedPartitions = function(query, partitionList,
                                  genomeSize=NULL, remainder="intergenic",
                                  bpProportion=FALSE) {
  .validateInputs(list(query=c("GRanges", "GRangesList"),
                       partitionList="list"))
  if (methods::is(query, c("GRangesList"))) {
    # Recurse over each GRanges object
    x = lapply(query, calcExpectedPartitions, partitionList,
               genomeSize, remainder,bpProportion)
    nameList = names(query)
    if(is.null(nameList)) {
      nameList = seq_along(query) # Fallback to sequential numbers
    }
    # Append names
    xb = data.table::rbindlist(x)
    xb$name = rep(nameList, vapply(x, nrow, integer(1)))
    return(xb)
  }

  # Get expected partitions - total number of bp each element
  # contributes to the genome
  #(intron = gene - (exon + 3'UTR + 5'UTR))
  #intergenic = genomeSize - other
  partitionNames = names(partitionList)

  if (bpProportion){
    # get the total number of base pairs in query
    query_total = sum(width(query))
  } else {
    # get the number of elements in query
    query_total = length(query)
  }

  # calculate the number of bp in each partition
  widths = lapply(partitionList, width)
  elements_total = lapply(widths, sum)

  partitionCounts = data.table::data.table(
    plyr::ldply(elements_total, data.frame))
  colnames(partitionCounts) = c("partition", "N")

  if (!is.null(genomeSize)) {
    # Calculate remainder
    partitionCounts = rbind(partitionCounts,
                            data.table::data.table(partition=remainder,
                                                   N=(genomeSize-sum(partitionCounts$N))))
  }

  partitionCounts$N = partitionCounts$N / genomeSize * query_total
  # these are the final expected partition counts based on element sizes
  partitionCounts = partitionCounts[order(partitionCounts$partition)]

  observedPartition = calcPartitions(query, partitionList, remainder, bpProportion)
  if (bpProportion) {
    # if flagged os bpProportions=TRUE - the output has different column names
    # now create data frame with observed and expected overlaps
    expectedPartitions = data.table::data.table(
      partition=observedPartition$partition,
      observed=observedPartition$bpOverlap)
  } else {
    expectedPartitions = data.table::data.table(
      partition=observedPartition$partition,
      observed=observedPartition$Freq)
  }
  expectedPartitions = expectedPartitions[partitionCounts, on = "partition", nomatch=0]
  data.table::setnames(expectedPartitions,"N","expected")
  expectedPartitions[,log10OE:=log10(expectedPartitions$observed/
                                       expectedPartitions$expected)]
  partitionNames = c(partitionNames, remainder)
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
#' partitionList = genomePartitionList(geneModels_hg19$genesGR,
#'                                     geneModels_hg19$exonsGR,
#'                                     geneModels_hg19$threeUTRGR,
#'                                     geneModels_hg19$fiveUTRGR)
#' calcCumulativePartitions(vistaEnhancers, partitionList)
calcCumulativePartitions = function(query, partitionList,
                                    remainder="intergenic") {
    .validateInputs(list(query=c("GRanges", "GRangesList"),
                         partitionList="list"))
    if (methods::is(query, c("GRangesList"))) {
        # Recurse over each GRanges object
        x = lapply(query, calcCumulativePartitions, partitionList, remainder)
        nameList = names(query)
        if(is.null(nameList)) {
            nameList = seq_along(query) # Fallback to sequential numbers
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
    for (pi in seq_along(partitionList)) {
        #message(partitionNames[pi],":")
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
        h = nrow(hits)
        #message("\tfound ", h)


        # Make mutually exclusive; remove hits from query
        query = query[-hits$xid]
        # Isolate hits
        hits = unique(hits, by="yid")
        # Link to positional data for feature of interest
        x = data.table::data.table(partition=partitionNames[pi],
                                size=if (h!=0) size=hits$size else size=0,
                                count=if (h!=0) count=hits$count else count=0)
        x = x[order(x$count, x$size),]
        x$cumsum  = cumsum(x$count)
        x$cumsize = cumsum(x$size)
        x$frif    = x$cumsum/query_total
        frif = rbind(frif, x)
    }
    # Create remainder...
    #message(remainder,":")
    x = data.table::data.table(partition=remainder,
                               size=as.numeric(width(query)))
    #message("\tfound ", length(query))
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
            # Fallback to sequential numbers
            nameList = seq_along(assignedPartitions)
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
#' p = calcCumulativePartitionsRef(vistaEnhancers, "hg19")
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
            facet_wrap(. ~name)
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
                plot_labels[,partition:=seq_len(nrow(partition_sizes))]
                assignedPartitions[,
                    partition:=rep(seq_len(nrow(partition_sizes)),
                                   partition_sizes$num_feats)]
            }
        }
    } else {
        if (!"partition" %in% colnames(plot_labels)) {
            plot_labels[,partition:=seq_len(nrow(partition_sizes))]
            assignedPartitions[,partition:=rep(seq_len(nrow(partition_sizes)),
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
#' p = calcExpectedPartitionsRef(vistaEnhancers, "hg19")
#' expectedPlot = plotExpectedPartitions(p)
#' pp = calcExpectedPartitionsRef(vistaEnhancers, "hg19", bpProportion=TRUE)
#' expectedPropPlot = plotExpectedPartitions(pp)
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
                expectedPartitions[,
                    partition:=seq_len(nrow(expectedPartitions))]
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
        theme_blank_facet_label() +
        theme(axis.line = element_line(size = 0.5),
              axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.5),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "transparent", color = NA),
              aspect.ratio = 1,
              legend.position = "bottom",
              plot.title = element_text(hjust = 0.5),
              panel.border = element_rect(colour = "black", fill=NA,
                                          size=0.5)
        ) +
        scale_fill_discrete(name="User set")

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
#' @param numbers logical indicating whether raw overlaps should be
#'     plotted instead of the default percentages
#' @param stacked logical indicating that data should be plotted as stacked
#'     bar plot
#' @return A ggplot object using a barplot to show the distribution
#'     of the query
#'  regions across a given partition list.
#' @export
#' @examples
#' p = calcPartitionsRef(vistaEnhancers, "hg19")
#' partPlot = plotPartitions(p)
#' partCounts = plotPartitions(p, numbers=TRUE)
#' partPlot = plotPartitions(p, stacked=TRUE)
#' pProp = calcPartitionsRef(vistaEnhancers, "hg19", bpProportion=TRUE)
#' partPlot = plotPartitions(pProp)
#' partCounts = plotPartitions(pProp, numbers=TRUE)
#' partPlot = plotPartitions(pProp, stacked=TRUE)
plotPartitions = function(assignedPartitions, numbers=FALSE, stacked=FALSE) {
    .validateInputs(list(assignedPartitions="data.frame"))

    if ("bpOverlap" %in% colnames(assignedPartitions)){
      df = data.frame(partition=assignedPartitions$partition,
                      Freq=assignedPartitions$bpOverlap)
      if("name" %in% names(assignedPartitions)){
        df$name = assignedPartitions$name
      }
    }else{
      df = assignedPartitions
    }
    # stacked bar option
    if (stacked){
      if ("name" %in% names(assignedPartitions)){
        # multiple datasets
        g = ggplot(df, aes(x=name, y=Freq, fill=partition))

        if (numbers) {
          g = g +
            geom_bar(stat="identity", position="stack")
        } else {
          g = g +
            geom_bar(stat="identity", position="fill")
        }

      } else {
        # single dataset stacked
        df$regionSet = "regionSet"
        g = ggplot(df, aes(x=regionSet, y=Freq, fill=partition))
        if (numbers) {
          g = g +
            geom_bar(stat="identity", position="stack")
        } else {
          g = g +
            geom_bar(stat="identity", position="fill")
        }
      }

      g = g +
        xlab("Region set") +
        ylab(ifelse(numbers,"Counts","Frequency"))
    } else {
      # not stacked
      # For multiple regions
      if ("name" %in% names(assignedPartitions)) {
        # percentages are to be set as the default instead of raw overlaps
        if (numbers == FALSE) {
          # assigned partitions is a data table
          df$Freq = (df$Freq / sum(df$Freq)) * 100
        }
        g = ggplot(df, aes(x=partition, y=Freq, fill=factor(name)))
      } else {
        # not a data table, a single regionset df
        if (numbers == FALSE) {
          df$Freq = (df$Freq / sum(df$Freq)) * 100
        }
        g = ggplot(df, aes(x=partition, y=Freq))
      }

      g = g +
        geom_bar(stat="identity", position=position_dodge()) +
        theme_blank_facet_label() + # No boxes around labels
        xlab("Genomic partition") +
        ylab(ifelse(numbers & "bpOverlap" %in% colnames(assignedPartitions),
                    "Counts [bp]",
                    ifelse(numbers, "Counts [regions]", "Frequency (%)")))  +
        scale_fill_discrete(name="regionSet") +
        theme(legend.position="bottom")
    }

    g = g +
      theme_classic()+
      theme(aspect.ratio=1) +
      theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) +
      theme(plot.title=element_text(hjust = 0.5)) + # Center title
      ggtitle(paste("Distribution across genomic partitions"))

    return(g)
}


# Internal helper function to overlap two data.table objects and
# get the total sum of their overlaps in base pairs
#
# @param partiton GRanges object1
# @param query GRanges object2
# @return A numeric value - sum of overlaps between
# partition and query
overlapWidths = function(partition, query){
  .validateInputs(list(query=c("GRanges"),
                       partition="GRanges"))
  partitionDT = grToDt(partition)
  queryDT = grToDt(query)
  setkey(queryDT, chr, start, end)

  # find overlaps
  hits = foverlaps(partitionDT, queryDT, nomatch=NULL)

  # get the widths of the overlaps - get maximum start
  # value and minumum end value, get their difference
  hits[, maxStart:=max(start, i.start), by=seq_len(nrow(hits))]
  hits[, minEnd:=min(end, i.end), by=seq_len(nrow(hits))]
  hits[, overlap:=minEnd-maxStart+1]

  # get total number of overlapping bases
  totalOverlap = sum(hits[, overlap])
  return(totalOverlap)
}





