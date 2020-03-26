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

# Function uses data.table rolling join to identify the nearest features
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


#' Converts a nucleotide count into a label with abbreviation
#' @param x base count
#' @return A label with 'kb' or 'mb' appended if appropriate
genomeLabel = function(x) {
	lab = x
	if (abs(x) > 1e6){

		lab = paste0(round(x/1e6), " mb")
	}
	else if (abs(x) > 1e3){
		lab = paste0(round(x/1e3), " kb")
	}
	return(lab)
}


#' Plots a histogram of distances to genomic features
#' 
#' Given the results from \code{featureDistribution}, plots a histogram of
#' distances surrounding the features of interest
#' 
#' @param dists Results from \code{featureDistribution}
#' @param bgdists Background distances. If provided, will plot a background
#'     distribution of expected distances 
#' @param featureName Character vector for plot labels (optional).
#' @param divisions A vector that defines the bin boundaries for the histogram
#' @param numbers a logical indicating whether the raw numbers should be 
#'     displayed, rather than percentages (optional).
#' @return A ggplot2 plot object
#' @export
#' @examples
#' queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' query = rtracklayer::import(queryFile)
#' TSSdist = calcFeatureDistRefTSS(query, "hg19")
#' plotFeatureDist(TSSdist, featureName="TSS")
plotFeatureDist = function(dists, bgdists=NULL, featureName="features", divisions=NULL, 
                           numbers=FALSE) {
	if (is.null(divisions)) {
		poscuts = seq(0,100000, by=2000)
		divisions = sort(unique(c(-poscuts, poscuts)))
	}
	df = cutDists(dists, divisions)
	# We could scale
	# df$Freq = scale(df$Freq, center=FALSE)
	if(is.list(dists)){
		nplots = length(dists)
	} else {
		nplots = 1
	}

	if (!is.null(bgdists)) {
		bgDistsDF = cutDists(bgDists, divisions)
		# bgDistsDF$Freq= scale(bgDistsDF$Freq, center=FALSE)
		bgDistsDF$Freq = (bgDistsDF$Freq / sum(bgDistsDF$Freq)) * 100
		df$bgFreq = rep(bgDistsDF$Freq, nplots)
		df$bgX = rep(seq_len(length(divisions)-1), nplots)
	}

	if ("name" %in% names(df)){
	    if (!numbers)
	        df$Freq = df[, .(Freq.Per = (Freq / sum(Freq)) * 100), 
	                     by = name][, "Freq.Per"]
		# It has multiple regions
		g = ggplot(df, aes(x=cuts, y=Freq, fill=name)) + 
			facet_grid(. ~name)
	} else {
	    if (!numbers)
	        df$Freq = (df$Freq / sum(df$Freq)) * 100
		g = ggplot(df, aes(x=cuts, y=Freq))
	}

	if (!is.null(bgdists)) {

		# bgtrack = scale(smooth(bgDistsDF$Freq), center=FALSE)
		g = g + 
			geom_line(stat="identity", aes(x=bgX,y=bgFreq), color="gray", alpha=1, size=1.5) + 
			geom_bar(stat="identity", aes(x=cuts,y=bgFreq), fill="gray", alpha=0.8)
	}

	# find midpoint
	midx = length(divisions)/2
	barcount = length(divisions)
	minlabel = genomeLabel(min(divisions))
	maxlabel = genomeLabel(max(divisions))
	edgeLabels = c(minlabel, rep("", barcount-3), maxlabel)
	g = g +
		geom_bar(data=df, stat="identity", fill="darkblue", alpha=0.7) + 
		geom_point(aes(x=midx, y=0), color="tan2", size=2, shape=17, alpha=0.8) +
		guides(fill=FALSE) + # remove legend for geom_point
		theme_classic() + 
		theme(aspect.ratio=1) + 
		theme_blank_facet_label() + 
		xlab(paste("Distance to", featureName)) +
		ylab(ifelse(numbers,"Counts","Frequency (%)")) +
		# theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
		theme(axis.text.x=element_text(angle = 0, hjust = c(0,1), vjust=0.5)) + # vlab()
		theme(plot.title = element_text(hjust = 0.5)) + # Center title
		ggtitle(paste("Distribution relative to", featureName)) +
		theme(legend.position="bottom") + 
		scale_x_discrete(labels=edgeLabels)

	return(g)
}


# Internal helper function for \code{plotFeatureDist}
cutDists = function(dists, divisions=NULL) {
	if (is.null(divisions)) {
		message("hi")
		poscuts = seq(0,100000, by=2000)
		divisions = sort(unique(c(-poscuts, poscuts)))
	}
	if (is.list(dists)) {
		x = lapply(dists, cutDists, divisions)

		# To accommodate multiple lists, we'll need to introduce a new 'name'
		# column to distinguish them.
		nameList = names(dists)
		if(is.null(nameList)) {
			nameList = 1:length(query) # Fallback to sequential numbers
		}

		# Append names
		xb = rbindlist(x)
		xb$name = rep(nameList, vapply(x, nrow, integer(1)))

		return(xb)
	}

	labels = labelCuts(sort(divisions), collapse=" to ", infBins=TRUE)
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

#' Internal helper function for \code{plotcFRiF} and \code{plotFRiF}
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
#' @param query A GenomicRanges or GenomicRangesList object with query regions.
#'              The name of each GRanges object is the default feature name.
#'              Each GRanges object must include a counts metadata column.
#' @param num_reads Numeric value representing the total number of aligned
#'                  reads/bases
#' @param genome_size Numeric value representing the size of a genome in bp
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the GenomicRanges or GenomicRangesList 
#'                      object.
#' @keywords FRiF
#' @export
#' @examples
#' data("promoter")
#' data("promoter_flanking")
#' data("exon")
#' data("intron")
#' data("utr3")
#' data("utr5")
#' plotFRiF(query = c("promoter", "promoter_flanking", "exon",
#'                     "intron", "utr3", "utr5"),
#'          num_reads=87520, genome_size = 3099922541, 
#'          feature_names = c("promoter", "promoter_flanking", "exon",
#'                            "intron", "utr3", "utr5")
#'          )
#' @export
plotFRiF = function(query, num_reads, genome_size, feature_names = NA) {
    palette = colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    # Calculate the FRiF for each feature type
    if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		frif = lapply(query, calcFRiF, num_reads)

        # Generate plot labels and options
        labels = lapply(frif, setLabels)
        # Collapse to list
        labels = data.table::rbindlist(labels, idcol=TRUE)
        feature_lengths = data.table::data.table(num_feats=elementNROWS(frif))

        # Identify expected values
        feature_dist = lapply(frif, getExpectedFeatures, genome_size)
        # Collapse to data.table with column for feature names
        feature_dist = data.table::rbindlist(feature_dist, idcol=TRUE)     
	} else {
        frif   = calcFRiF(query, num_reads)
        labels = setLabels(frif)
        feature_dist    = getExpectedFeatures(frif, genome_size)
        feature_lengths = data.table::data.table(num_feats=length(query))
    }

    plot_colors = palette(nrow(feature_lengths))

    # If name vector provided, update names
    if (all(!is.na(feature_names))) {
        if (length(feature_names) == nrow(feature_lengths)) {
            labels[,.id:=feature_names]
            feature_dist[,.id:=feature_names]
        } else {
            if (!".id" %in% colnames(labels)) {
                labels[,.id:=seq(1:nrow(feature_lengths))]
                feature_dist[,.id:=seq(1:nrow(feature_lengths))]
            }
        }
    } else {
        if (!".id" %in% colnames(labels)) {
            labels[,.id:=seq(1:nrow(feature_lengths))]
            feature_dist[,.id:=seq(1:nrow(feature_lengths))]
        }
    }
    
    # Set colors
    labels[,color:=plot_colors]

    # Finalize feature_dist table
    feature_dist$observed = as.numeric(labels$val)
    feature_dist$logOE    = log10(feature_dist$observed/feature_dist$expected)
    feature_dist$logOE    = ifelse(feature_dist$logOE < 0,0,feature_dist$logOE)
    feature_dist = merge(feature_dist, labels, by=".id")
    feature_dist = feature_dist[order(feature_dist$logOE),]
    feature_dist$.id   = factor(feature_dist$.id, levels=feature_dist$.id)
    feature_dist$color = factor(feature_dist$color, levels=feature_dist$color)

    p = ggplot(feature_dist, aes(x = .id, y = logOE)) +
            geom_bar(stat="identity",
                     fill = feature_dist$color,
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
#' @param query A GenomicRanges or GenomicRangesList object with query regions.
#'              The name of each GRanges object is the default feature name.
#'              Each GRanges object must include a counts metadata column.
#' @param num_reads Numeric value representing the total number of aligned
#'                  reads/bases
#' @param feature_names An optional character vector of feature names, in the 
#'                      same order as the GenomicRanges or GenomicRangesList 
#'                      object.
#' @keywords cFRiF
#' @export
#' @examples
#' data("promoter")
#' data("promoter_flanking")
#' data("exon")
#' data("intron")
#' data("utr3")
#' data("utr5")
#' plotcFRiF(query = c("promoter", "promoter_flanking", "exon",
#'                    "intron", "utr3", "utr5"),
#'           num_reads=87520,
#'           names = c("promoter", "promoter_flanking", "exon",
#'                    "intron", "utr3", "utr5")
#'          )
#' @export
plotcFRiF = function(query, num_reads, feature_names = NA) {
    palette = colorRampPalette(c("#999999", "#FFC107", "#27C6AB", "#004D40",
                                  "#B97BC8", "#009E73", "#C92404", "#E3E550",
                                  "#372B4C", "#E3DAC7", "#27CAE6", "#B361BC",
                                  "#897779", "#6114F8", "#19C42B", "#56B4E9"))
    # Calculate the FRiF for each feature type
    if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		frif = lapply(query, calcFRiF, num_reads)
        # Generate plot labels and options
        labels = lapply(frif, setLabels)
        # Collapse to list
        labels = data.table::rbindlist(labels, idcol=TRUE)
        feature_lengths = data.table::data.table(num_feats=elementNROWS(frif))
        # Collapse frif list to data.table with column for feature names
        frif = data.table::rbindlist(frif, idcol=TRUE)
	} else {
        frif = calcFRiF(query, num_reads)
        labels = setLabels(frif)
        feature_lengths = data.table::data.table(num_feats=length(query))
    }

    plot_colors = palette(nrow(feature_lengths))

    # If name vector provided, update names
    if (all(!is.na(feature_names))) {
        if (length(feature_names) == nrow(feature_lengths)) {
            labels[,.id:=feature_names]
            frif[,.id:=rep(feature_names, each=feature_lengths$num_feats)]
        } else {
            if (!".id" %in% colnames(frif)) {
                labels[,.id:=seq(1:nrow(feature_lengths))]
                frif[,.id:=rep(seq(1:nrow(feature_lengths)),
                     feature_lengths$num_feats)]
            }
        }
    } else {
        if (!".id" %in% colnames(frif)) {
            labels[,.id:=seq(1:nrow(feature_lengths))]
            frif[,.id:=rep(seq(1:nrow(feature_lengths)),
                 feature_lengths$num_feats)]
        }
    }

    # Set colors
    labels[,color:=plot_colors]

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
    p = p + scale_color_manual(labels=paste0(labels$.id, ": ", labels$val),
                               values=labels$color) +
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