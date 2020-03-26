#' Calculate the Fraction of Regions in Features (FRiF)
#'
#' This function calculates the fraction of regions in a feature set and
#' returns the cumulative sum of reads, cumulative size of covered features, 
#' the fraction of reads in those features, and the number of total features.
#'
#' @param query A GenomicRanges object with query regions and a metadata
#'              column named 'counts' which represents the number of reads
#'              or bases that aligned to the corresponding regions.
#' @param total A numeric value presenting the total number of aligned 
#'              reads/bases.
#' @keywords FRiF
#' @export
calcFRiF = function(query, total) {
    if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, calcFRiF, total)
		return(x)
	}

    # Total must be some number greater than 0
    if (total <= 0) {
        warning("The total number (total:", total, ") of aligned reads or ",
                "bases must be greater than 0.")
        return(NULL)
    }

    # This is past logic recreated, just starting with a GRanges object now
    bedFile  = data.table::data.table(chr=as.character(seqnames(query)),
                                      start=start(query),
                                      end=end(query),
                                      counts=mcols(query)$count)
    data.table::setkey(bedFile, chr, start, end)
    # Create a reduced region set to remove overlap
    reduceGR = reduce(query)
    # Build the reduced region set with combined counts
    reduceDT = data.table::data.table(chr=as.character(seqnames(reduceGR)),
                                      start=start(reduceGR),
                                      end=end(reduceGR))
    data.table::setkey(reduceDT, chr, start, end)
    # Identify which regions overlap the reduced region set
    hitsGR   = findOverlaps(query=reduceGR, subject=query)
    hits     = data.table::data.table(xid=queryHits(hitsGR),
                                      yid=subjectHits(hitsGR))
    # hits   = foverlaps(reduceDT, bedFile, by.x=c("chr", "start", "end"),
    #                    type="any", which=TRUE, nomatch=0)
    counts = data.table::data.table(index=rep(1:nrow(bedFile)),
                                    counts=bedFile$counts)
    data.table::setkey(hits, yid)
    data.table::setkey(counts, index)
    out = hits[counts, nomatch=0]
    # For reduced regions with multiple hits, sum the counts
    out[, countsSum:= sum(counts), by=xid]
    
    # Need an index column to combine
    reduceDT$index = rep(1:nrow(reduceDT))
    data.table::setkey(reduceDT, index)
    counts = data.table::data.table(xid=out$xid, counts=out$countsSum)
    data.table::setkey(counts, xid)
    # Add the combined counts column
    reduceDT = reduceDT[counts, nomatch=0]
    # Remove copied rows
    reduceDT = reduceDT[reduceDT[,.I[which.max(counts)],by=index]$V1]
    # Remove index column
    reduceDT = reduceDT[,index:=NULL]
    # Add a region size column
    reduceDT$size = (reduceDT$end - reduceDT$start)
    # Sort by the counts column
    reduceDT = reduceDT[order(-reduceDT$counts),]
    # Remove regions that contain no counts
    reduceDT = reduceDT[apply(reduceDT != 0, 1, all),]
    # Establish the cumulative sum of counts for the regions in this feature
    reduceDT = cbind(reduceDT, cumsum=cumsum(reduceDT$counts))
    # Calculate the cumulative size of the regions in a feature
    reduceDT = cbind(reduceDT, cumsize=cumsum(reduceDT$size))
    # Calculate the cumulative fraction of counts in this feature
    reduceDT = cbind(reduceDT, frip=reduceDT$cumsum/as.numeric(total))
    # Populate the cumulative number of observed regions of this feature type
    reduceDT = cbind(reduceDT, numfeats=as.numeric(1:nrow(reduceDT)))

    return(reduceDT)
}

#' Internal helper function for \code{plotcFRiF} and 
#' \code{calcExpectedDistribution}
#' @param frif Results from \code{calcFRiF}
setLabels = function(frif) {
    if (is(frif, "list")) {
		# Recurse over each GRanges object
		x = lapply(frif, setLabels)
        x = data.table::rbindlist(x, idcol=TRUE)
		return(x)
	}
    return(data.table::data.table(xPos=0.95*max(log10(frif$cumsize)),
                                  yPos=max(frif$frip)+0.001,
                                  val=round(max(frif$frip),2),
                                  color=NA,
                                  stringsAsFactors=FALSE))
}

#' Calculate the expected fraction of reads in features.
#'
#' @param frif Results from \code{calcFRiF}
#' @param genome_size Numeric value representing the size of the genome in
#'                    bases.
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the GenomicRanges or GenomicRangesList 
#'                      object.
#' @export
calcExpectedDistribution = function(frif, genome_size, feature_names = NA) {
    if (is(frif, "list")) {
		# Recurse over each GRanges object
		x = lapply(frif, calcExpectedDistribution, genome_size)
        x = data.table::rbindlist(x, idcol=TRUE)
		return(x)
	}
    feature_dist = data.table::data.table(numfeats=max(frif$numfeats),
                                          numbases=max(frif$cumsize),
                                          expected=(max(frif$cumsize)/
                                                    genome_size),
                                          stringsAsFactors=FALSE)
    plot_labels = setLabels(frif)
    # Add observed fractions
    feature_dist$observed = as.numeric(plot_labels$val)
    # Calculate the log observed/expected
    feature_dist$logOE    = log10(feature_dist$observed/feature_dist$expected)
    feature_dist$logOE    = ifelse(feature_dist$logOE < 0,0,feature_dist$logOE)
    feature_dist = feature_dist[order(feature_dist$logOE),]
    
    return(feature_dist)
}

#' Plot Fraction of Reads in Features (FRiF)
#'
#' This function plots the log10 of the observed/expected distribution of
#' reads/bases in a set of features.
#'
#' @param dists Results from \code{calcExpectedDistribution}
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the dists object.
#' @keywords FRiF
#' @export
plotFRiF = function(dists, feature_names = NA) {
    palette = colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    plot_colors = palette(nrow(dists))

    if (is(dists, "list")) {
        dists = data.table::rbindlist(dists, idcol=TRUE)
    }
    
    # Set colors
    dists[,color:=plot_colors]

    # If name vector provided, update names
    if (all(!is.na(feature_names))) {
        if (length(feature_names) == nrow(dists)) {
            dists[,.id:=feature_names]
        } else {
            if (!".id" %in% colnames(plot_labels)) {
                dists[,.id:=seq(1:nrow(dists))]
            }
        }
    }
    
    dists = dists[order(dists$logOE),]
    dists$.id   = factor(dists$.id, levels=dists$.id)
    dists$color = factor(dists$color, levels=dists$color)

    p = ggplot(dists, aes(x = .id, y = logOE)) +
            geom_bar(stat="identity",
                     fill = dists$color,
                     alpha = 0.5) + 
            geom_hline(aes(yintercept=0), linetype="dotted") +
            xlab('') +
            ylab(expression(log[10](over(Obs, Exp)))) +
            coord_flip() +
            theme(axis.line = element_line(size = 0.5),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
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

#' Plot cumulative Fraction of Reads in Features (cFRiF)
#'
#' This function plots the cumulative fraction of reads in a set of features
#'
#' @param frif  Results from \code{calcFRiF}
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the GenomicRanges or GenomicRangesList 
#'                      object.
#' @keywords cFRiF
#' @export
plotcFRiF = function(frif, feature_names = NA) {
    plot_labels = setLabels(frif)
    palette = colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    plot_colors = palette(length(frif))
    # Set colors
    plot_labels[,color:=plot_colors]
    
    if (is(frif, "list")) {
        feature_lengths = data.table::data.table(num_feats=elementNROWS(frif))
    } else {
        feature_lengths = data.table::data.table(num_feats=length(frif))
    }
    
    if (is(frif, "list")) {
        frif = data.table::rbindlist(frif, idcol=TRUE)
    }

    # If name vector provided, update names
    if (all(!is.na(feature_names))) {
        if (length(feature_names) == nrow(feature_lengths)) {
            plot_labels[,.id:=feature_names]
            frif[,.id:=rep(feature_names, each=feature_lengths$num_feats)]
        } else {
            if (!".id" %in% colnames(plot_labels)) {
                plot_labels[,.id:=seq(1:nrow(feature_lengths))]
                frif[,.id:=rep(seq(1:nrow(feature_lengths)),
                     feature_lengths$num_feats)]
            }
        }
    } else {
        if (!".id" %in% colnames(plot_labels)) {
            plot_labels[,.id:=seq(1:nrow(feature_lengths))]
            frif[,.id:=rep(seq(1:nrow(feature_lengths)),
                     feature_lengths$num_feats)]
        }
    }

    p = ggplot(frif, aes(x=log10(cumsize), y=frip, group=.id, color=.id)) +
        geom_line(size=2, alpha=0.5) +
        guides(linetype = FALSE) +
        labs(x=expression(log[10]("number of bases")), y="cFRiF") +
        theme(axis.line = element_line(size = 0.5),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
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
              panel.border = element_rect(colour = "black", fill=NA, size=0.5)
        )

    # Recolor and reposition legend
    p = p + scale_color_manual(
            labels=paste0(plot_labels$.id, ": ", plot_labels$val),
            values=plot_labels$color) +
        labs(color="cFRiF") +
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