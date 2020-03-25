

# Checks to make sure a BSgenome object is installed,
# and if so, returns it. If the genome is not installed, it issues a warning
# and returns NULL.
.requireAndReturn = function(BSgenomeString) {
	if (requireNamespace(BSgenomeString))
		return(utils::getAnywhere(BSgenomeString)$objs[[1]])
	else
		warning(BSgenomeString, " is not installed")
		return(NULL)
}


# Efficiently split a data.table by a column in the table
# 
# @param DT Data.table to split
# @param splitFactor Column to split, which can be a character vector
#	or an integer.
# @return	List of data.table objects, split by column
# @examples
# DT = data.table::data.table(letters, grp = rep(c("group1", "group2"), 13))
# splitDataTable(DT, "grp")
# splitDataTable(DT, 2)
splitDataTable = function(DT, splitFactor) {
	if (is.numeric(splitFactor)) {
		splitFactor = colnames(DT)[splitFactor]
		message("Integer splitFactor, changed to: ", splitFactor)
	}
	lapply( split(1:nrow(DT), DT[, get(splitFactor)]), function(x) DT[x])
}



#Two utility functions for converting data.tables into GRanges objects
#genes = dtToGR(gModels, "chr", "txStart", "txEnd", "strand", "geneId");
dtToGrInternal = function(DT, chr, start, end=NA, strand=NA, name=NA, metaCols=NA) {
	if (is.na(end)) {
		if ("end" %in% colnames(DT)) {
			end = "end"
		} else {
			end = start
		}
	}
	if (is.na(strand)) {
		gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), strand="*")
	} else {
		# GRanges can only handle '*' for no strand, so replace any non-accepted
		# characters with '*'
		DT[,strand:=as.character(strand)]
		DT[strand=="1", strand:="+"]
		DT[strand=="-1", strand:="-"]
		DT[[`strand`]] =  gsub("[^+-]", "*", DT[[`strand`]])
		gr=GRanges(seqnames=DT[[`chr`]], ranges=IRanges(start=DT[[`start`]], end=DT[[`end`]]), strand=DT[[`strand`]])
	}
	if (! is.na(name) ) {
		names(gr) = DT[[`name`]]
	} else {
		names(gr) = 1:length(gr)
	}
	if(! is.na(metaCols)) {
		for(x in metaCols) {
			elementMetadata(gr)[[`x`]]=DT[[`x`]]
		}
	}
	gr
}


# Converts a data.table (DT) object to a GenomicRanges (GR) object. Tries to be
# intelligent, guessing chr and start, but you have to supply end or other
# columns if you want them to be carried into the GR.
dtToGr = function(DT, chr="chr", start="start", end=NA, strand=NA, name=NA, splitFactor=NA, metaCols=NA) {
	if(is.na(splitFactor)) {
		return(dtToGrInternal(DT, chr, start, end, strand, name,metaCols))
	}

	if ( length(splitFactor) == 1 ) { 
		if( splitFactor %in% colnames(DT) ) {
			splitFactor = DT[, get(splitFactor)]
		}
	}

	lapply(split(1:nrow(DT), splitFactor), function(x) { 
			dtToGrInternal(DT[x,], chr, start, end, strand, name,metaCols)
		}
	)


}

# Convert a GenomicRanges into a data.table
grToDt = function(GR) {
	DF=as.data.frame(elementMetadata(GR))
	if( ncol(DF) > 0) {
		DT = data.table(chr=as.vector(seqnames(GR)), start=start(GR), end=end(GR), DF)
	} else {
		DT = data.table(chr=as.vector(seqnames(GR)), start=start(GR), end=end(GR))
	}
	return(DT)
}


# Converts a list of data.tables (From BSreadbeds) into GRanges.
BSdtToGRanges = function(dtList) {
	gList = list()
	for (i in 1:length(dtList)) {
		#dt = dtList[[i]]
		setkey(dtList[[i]], chr, start)
		#convert the data into granges object
		gList[[i]] = GRanges(seqnames=dtList[[i]]$chr, ranges=IRanges(start=dtList[[i]]$start, end=dtList[[i]]$start), strand=rep("*", nrow(dtList[[i]])), hitCount=dtList[[i]]$hitCount, readCount=dtList[[i]]$readCount)
#I used to use end=start+1, but this targets CG instead of just a C, and it's causing edge-effects problems when I assign Cs to tiled windows using (within). Aug 2014 I'm changing to start/end at the same coordinate.
	}
	return(gList)
}


# Clear ggplot face label
# Usually ggplot2 facets are labeled with boxes surrounding the label. This
# function removes the box, so it's a simple label for each facet.
theme_blank_facet_label = function() {
	return(theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		strip.background = element_blank()
		)
	)
}


# Creates labels based on a discretization definition.
# 
# If you are building a histogram of binned values, you want to have labels for
# your bins that correspond to the ranges you used to bin. This function takes
# the breakpoints that define your bins and produces nice-looking labels for
# your histogram plot.
# 
# \code{labelCuts} will take a cut group, (e.g., a quantile division of 
# some signal), and give you clean labels (similar to the cut method).
# @param breakPoints The exact values you want as boundaries for your bins
# @param round_digits Number of digits to cut round labels to. 
# @param signif_digits Number of significant digits to specify. 
# @param collapse Character to separate the labels
# @param infBins use >/< as labels on the edge bins
# @examples 
# labelCuts(seq(0,100,by=20))
labelCuts = function(breakPoints, round_digits=1, signif_digits=3, collapse="-", infBins=FALSE) {
      roundedLabels = signif(round(
      	cbind( breakPoints[-length(breakPoints)],breakPoints[-1]),digits), signif_digits)
      # set the Inf values to NA so formatC can add commas
      is.na(roundedLabels) = sapply(roundedLabels, is.infinite) 
      labelsWithCommas = formatC(roundedLabels, format="d", big.mark=",")
      labels = apply(labelsWithCommas, 1, paste0, collapse=collapse) 
      if (infBins) {
        labels[1] = paste0("<", formatC(breakPoints[2], format="d", big.mark=","))
        labels[length(labels)] = paste0(">", formatC(breakPoints[length(breakPoints)-1], format="d", big.mark=","))
      }
      return(labels)
}

# labelCuts = function(breakPoints, digits=1, collapse="-", infBins=FALSE) {
# 	labels = 
# 	apply(round(cbind( breakPoints[-length(breakPoints)],	
# 		breakPoints[-1]),digits), 1, paste0, collapse=collapse) 

# 	if (infBins) {
# 		labels[1] = paste0("<", breakPoints[2])
# 		labels[length(labels)] = paste0(">", breakPoints[length(breakPoints)-1])
# 	}
# 	return(labels)
# }
