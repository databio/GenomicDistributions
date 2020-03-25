# Old, slow version based on GRanges methods
calcFeatureDistBioc = function(query, features) {
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

	preDist[!preIndNA] = -distance(query[!preIndNA], features[precedeInd[!preIndNA]])

	postDist = rep(NA, length(query))
	postDist[!folIndNA] = distance(query[!folIndNA], features[followInd[!folIndNA]])

	postHits = -preDist > postDist
	postHitsNA = is.na(postHits)
	dists = preDist
	dists[postHits[!postHitsNA]] = postDist[postHits[!postHitsNA]]
	return(dists)
}

#' Find the distance to the nearest genomic feature
#' 
#' For a given query set of genomic regions, and a given feature set of regions,
#' this function will return the distance for each query region to its closest
#' feature. It ignores strand and returns the distance as positive or negative,
#' depending on whether the feature is upstream or downstream
#' 
#' This function is similar to the bioconductor distanceToNearest function, but
#' returns negative values for downstream distances instead of absolute values.
#' This allows you to assess the relative location.
#' 
#' @param query A GRanges or GRangesList object with query sets
#' @param features A GRanges object with features to test distance to
#' 
#' @export
calcFeatureDist = function(query, features) {
	if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, calcFeatureDist, features)
		return(x)
	}
	queryDT = grToDt(query)
	featureDT = grToDt(features)
	queryDTs = splitDataTable(queryDT, "chr")
	featureDTs = splitDataTable(featureDT, "chr")
	as.vector(unlist(mapply(queryDTs, featureDTs[names(queryDTs)], FUN=DTNearest)))
}

# Function uses dat.table rolling join to identify the nearest features
# really quickly.
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
#' @export
calcFeatureDistRefTSS = function(query, refAssembly) {
	features = getTSSs(refAssembly)
	return(calcFeatureDist(query, features))
}

#' Plots a histogram of distances to genomic features
#' 
#' Given the results from \code{featureDistribution}, plots a histogram of
#' distances surrounding the features of interest
#' 
#' @param dists Results from \code{featureDistribution}
#' @param featureName Character vector for plot labels (optional).
#' @export
plotFeatureDist = function(dists, featureName="features") {

	df = cutDists(dists)
	if ("name" %in% names(df)){
		# It has multiple regions
		g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
			facet_grid(. ~name)
	} else {
		g = ggplot(df, aes(x=cuts, y=Freq))
	}

	g = g +
		geom_bar(stat="identity") + 
		geom_vline(xintercept = (length(unique(df$cuts))+1)/2, color="darkgreen") +
		theme_classic() + 
		theme(aspect.ratio=1) + 
		theme_blank_facet_label() + 
		xlab(paste("Distance to", featureName)) +
		ylab("Number of regions") +
		theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
		theme(plot.title = element_text(hjust = 0.5)) + # Center title
		ggtitle(paste("Distribution relative to", featureName)) +
		theme(legend.position="bottom")

	return(g)
}

# Internal helper function for \code{plotFeatureDist}
cutDists = function(dists, divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0, 100, 1000, 10000, 1e6, Inf)) {
	if (is.list(dists)) {
		x = lapply(dists, cutDists)

		# To accommodate multiple lists, we'll need to introduce a new 'name'
		# column to distinguish them.
		nameList = names(dists)
		if(is.null(nameList)) {
			nameList = 1:length(query) # Fallback to sequential numbers
		}

		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, sapply(x, nrow))

		return(xb)
	}

	labels = labelCuts(divisions, collapse=" to ", infBins=TRUE)
	cuts = cut(dists, divisions, labels)
	df = as.data.frame(table(cuts))
	return(df)
}


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
#' @examples
#' calcFRiF()
calcFRiF = function(query, total) {
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
    # Create a reduced region set to remove overlap
    reduceGR = reduce(query)
    # Identify which regions overlap the reduced region set
    hitsGR   = findOverlaps(query=reduceGR, subject=query)
    hits     = data.table::data.table(xid=queryHits(hitsGR),
                                      yid=subjectHits(hitsGR))
    counts   = data.table::data.table(index=rep(1:nrow(bedFile)),
                                      counts=bedFile$counts)
    setkey(hits, yid)
    setkey(counts, index)
    out = hits[counts, nomatch=0]
    # For reduced regions with multiple hits, sum the counts
    out[, countsSum:= sum(counts), by=xid]
    # Build the reduced region set with combined counts
    reduceDT = data.table::data.table(chr=as.character(seqnames(reduceGR)),
                                      start=start(reduceGR),
                                      end=end(reduceGR))
    # Need an index column to combine
    reduceDT$index = rep(1:nrow(reduceDT))
    setkey(reduceDT, index)
    counts = data.table::data.table(xid=out$xid, counts=out$countsSum)
    setkey(counts, xid)
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

#' Internal helper function for \code{plotFRiF}
#' @param query A GRanges object
setLabels = function(query) {
    return(data.table::data.table(
            xPos=0.95*max(log10(query$cumsize)),
            yPos=max(query$frip)+0.001,
            val=round(max(query$frip),2),
            color=NA,
            stringsAsFactors=FALSE)
    )
}

#' Internal helper function for \code{plotFRiF}
#' @param query A GRanges object.
#' @param genome_size Numeric value representing the size of the genome in
#'                    bases.
getExpectedFeatures = function(query, genome_size) {
    return(data.table::data.table(
            numfeats=max(query$numfeats),
            numbases=max(query$cumsize),
            expected=(max(query$cumsize)/genome_size),
            stringsAsFactors=FALSE)
    )
}

#' Plot Fraction of Reads in Features (FRiF)
#'
#' This function plots the fraction of reads in a set of features
#'
#' @param sample_name A character vector representing the name of a sample
#' @param num_reads Numeric value representing the number of aligned
#'                  reads/bases
#' @param genome_size Numeric value representing the size of a genome in bp
#' @param type A character vector representing the plot type to produce
#' @param output_name A character vector of the desired output file name
#' @param query A GenomicRangesList object with query regions.
#'              The name of each GRanges object is assumed to the feature name.
#'              Each GRanges object must include a counts metadata column.
#' @keywords cFRiF FRiF
#' @export
#' @examples
#' data("promoter")
#' data("promoter_flanking")
#' data("exon")
#' data("intron")
#' data("utr3")
#' data("utr5")
#' plotFRiF(sample_name="example", num_reads=87520,
#'          output_name="example_frif.pdf",
#'          query = c("promoter", "promoter_flanking", "exon",
#'                    "intron", "utr3", "utr5"))
#' @export
plotFRiF = function(sample_name, num_reads, genome_size,
                     type = c("cfrif", "frif", "both"),
                     reads=TRUE, output_name, query) {

    feature_dist  = data.table::data.table(feature=character(),
                                            numfeats=numeric(),
                                            numbases=numeric(),
                                            expected=numeric(),
                                            stringsAsFactors=FALSE)
    palette = colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    plot_colors = palette(length(query))

    # Calculate the FRiF for each feature type
    if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		frif = lapply(query, calcFRiF, num_reads)
        feature_dist = lapply(query, calcFRiF, num_reads)
	} else {
        frif = calcFRiF(query, num_reads)
    }
    
    # Generate plot labels and options
    labels = lapply(frif, setLabels)
    # Collapse to list
    labels = rbindlist(labels, idcol=TRUE)
    # Set colors
    labels[,color:=plot_colors]
    
    feature_dist = lapply(frif, getExpectedFeatures, genome_size)
    # Collapse to list
    feature_dist = rbindlist(feature_dist, idcol=TRUE)
    
    # Collapse frif list to data.table with column for feature names
    frif = rbindlist(frif, idcol=TRUE)

    # Finalize feature_dist table
    feature_dist$observed = as.numeric(labels$val)
    feature_dist$logOE = log10(feature_dist$observed/feature_dist$expected)
    feature_dist$logOE = ifelse(feature_dist$logOE < 0, 0, feature_dist$logOE)
    feature_dist = merge(feature_dist, labels, by=".id")
    feature_dist = feature_dist[order(feature_dist$logOE),]
    feature_dist$.id = factor(feature_dist$.id, levels=feature_dist$.id)
    feature_dist$color = factor(feature_dist$color, levels=feature_dist$color)

    if (tolower(type) == "both") {
        # Produce plot with bed files
        p = ggplot(frif, aes(x=log10(cumsize), y=frip,
                    group=.id, color=.id)) +
            geom_line(size=2, alpha=0.5) +
            guides(linetype = FALSE) +
            labs(x=expression(log[10]("number of bases")),
                 y="FRiF") +
            theme_PEPATAC()

        # Recolor and reposition legend
        p = p + scale_color_manual(labels=paste0(labels$.id, ": ",
                                                  labels$val),
                                    values=labels$color) +
            labs(color="FRiF") +
            theme(legend.position="right",
                  legend.justification=c(0.1,0.9),
                  legend.background=element_blank(),
                  legend.text = element_text(size = rel(0.65)),
                  legend.key = element_blank(),
                  axis.text.x = element_text(angle = 0, hjust = 1,
                                             vjust=0.5))

        p2 = ggplot(feature_dist, aes(x = .id, y = logOE)) +
            geom_bar(stat="identity", fill=labels$color, alpha=0.5) + 
            geom_hline(aes(yintercept=0), linetype="dotted") +
            xlab('') +
            ylab(expression(log[10](over(Obs, Exp)))) +
            coord_flip() +
            scale_x_discrete(position="top") +
            theme_PEPATAC() +
            theme(plot.background = element_rect(fill = "transparent",
                                                 color = NA,),
                  panel.background = element_rect(fill = "transparent"),
                  rect = element_rect(fill = "transparent"),
                  plot.margin = unit(c(0,0,-6.5,-6.5),"mm"))

        g   = ggplotGrob(p2)
        min_x = min(layer_scales(p)$x$range$range)
        max_x = max(layer_scales(p)$x$range$range)
        min_y = min(layer_scales(p)$y$range$range)
        max_y = max(layer_scales(p)$y$range$range)

        p = p + annotation_custom(grob = g, xmin = 1.05*min_x,
                                   xmax=min_x*2.05, ymin=max_y/2,
                                   ymax=max_y)
    } else if (tolower(type) == "cfrif") {
        p = ggplot(frif, aes(x=log10(cumsize), y=frip,
                    group=.id, color=.id)) +
            geom_line(size=2, alpha=0.5) +
            guides(linetype = FALSE) +
            labs(x=expression(log[10]("number of bases")), y="FRiF") +
            theme_PEPATAC()

        # Recolor and reposition legend
        p = p + scale_color_manual(labels=paste0(labels$.id, ": ",
                                                  labels$val),
                                    values=labels$color) +
            labs(color="FRiF") +
            theme(legend.position=c(0.075,0.975),
                  legend.justification=c(0.1,0.9),
                  legend.title = element_blank(),
                  legend.text = element_text(size = rel(0.65)), 
                  legend.background=element_blank(),
                  legend.key = element_blank(),
                  axis.text.x = element_text(angle = 0, hjust = 1,
                                             vjust=0.5))
    } else if (tolower(type) == "frif") {
        p = ggplot(feature_dist, aes(x = .id, y = logOE)) +
            geom_bar(stat="identity",
                     fill = feature_dist$color,
                     alpha = 0.5) + 
            geom_hline(aes(yintercept=0), linetype="dotted") +
            xlab('') +
            ylab(expression(log[10](over(Obs, Exp)))) +
            coord_flip() +
            theme_PEPATAC()
    } else {
        # default to both
        # Produce plot with bed files
        p = ggplot(frif,
                    aes(x=log10(cumsize), y=frip,
                        group=.id, color=.id)) +
            geom_line(aes(linetype=.id), size=2, alpha=0.5) +
            guides(linetype = FALSE) +
            labs(x=expression(log[10]("number of bases")),
                 y="FRiF") +
            theme_PEPATAC()

        # Recolor and reposition legend
        p = p + scale_color_manual(labels=paste0(labels$.id, ": ",
                                                  labels$val),
                                    values=labels$color) +
            labs(color="FRiF") +
            theme(legend.position="right",
                  legend.justification=c(0.1,0.9),
                  legend.background=element_blank(),
                  legend.key = element_blank(),
                  axis.text.x = element_text(angle = 0, hjust = 1,
                                             vjust=0.5))

        p2 = ggplot(feature_dist, aes(x = .id, y = logOE)) +
            geom_bar(stat="identity", fill=labels$color, alpha=0.5) + 
            geom_hline(aes(yintercept=0), linetype="dotted") +
            xlab('') +
            ylab(expression(log[10](over(Obs, Exp)))) +
            coord_flip() +
            scale_x_discrete(position="top") +
            theme_PEPATAC() +
            theme(plot.background = element_rect(fill = "transparent",
                                                 color = NA,),
                  panel.background = element_rect(fill = "transparent"),
                  rect = element_rect(fill = "transparent"),
                  plot.margin = unit(c(0,0,-6.5,-6.5),"mm"))

        g   = ggplotGrob(p2)
        min_x = min(layer_scales(p)$x$range$range)
        max_x = max(layer_scales(p)$x$range$range)
        min_y = min(layer_scales(p)$y$range$range)
        max_y = max(layer_scales(p)$y$range$range)

        p = p + annotation_custom(grob = g, xmin = 1.05*min_x,
                                   xmax=min_x*2.05, ymin=max_y/2,
                                   ymax=max_y)
    }

    if (!exists("p")) {
        p = ggplot()
    }

    return(p)
}