
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
#' @param feats A GRanges object with features to test distance to

#' @export
featureDistribution = function(query, feats) {
	if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, featureDistribution, feats)
		return(x)
	}

	precedeInd = precede(query, feats)
	preIndNA = is.na(precedeInd)
	followInd = follow(query, feats)
	folIndNA = is.na(followInd)
	preDist = rep(NA, length(query))

	preDist[!preIndNA] = -distance(query[!preIndNA], feats[precedeInd[!preIndNA]])

	postDist = rep(NA, length(query))
	postDist[!folIndNA] = distance(query[!folIndNA], feats[followInd[!folIndNA]])

	postHits = -preDist > postDist
	postHitsNA = is.na(postHits)
	dists = preDist
	dists[postHits[!postHitsNA]] = postDist[postHits[!postHitsNA]]
	return(dists)
}


#' @export
featureDistOverGenome = function(query, refAssembly) {
	if (is(query, "GRangesList")) {
		# Recurse over each GRanges object
		x = lapply(query, featureDistOverGenome, refAssembly)
		return(x)
	}

	feats = getTSSs(refAssembly)
	return(featureDistribution(query, feats))
}

#' Plots a histogram of distances to genomic features
#' 
#' Given the results from \code{featureDistribution}, plots a histogram of
#' distances surrounding the features of interest
#' 
#' @param dists Results from \code{featureDistribution}
#' @param plotTitle Title for plot.
#' @export
plotFeatureDist = function(dists, plotTitle="Distribution relative to features") {

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
		xlab("Distance to feature") +
		ylab("Number of regions") +
		theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
		theme(plot.title = element_text(hjust = 0.5)) + # Center title
		ggtitle(plotTitle) +
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

	labels = GenomicDistributions:::labelCuts(divisions, collapse=" to ", infBins=TRUE)
	cuts = cut(dists, divisions, labels)
	df = as.data.frame(table(cuts))
	return(df)
}