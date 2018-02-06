
#' Find the distance to the nearest genomic feature

#' For a given query set of genomic regions, and a given feature set of regions,
#' this function will return the distance for each query region to its closest
#' feature. It ignores strand and returns the distance as positive or negative,
#' depending on whether the feature is upstream or downstream


#' This function is similar to the bioconductor distanceToNearest function, but
#' returns negative values for downstream distances instead of absolute values.
#' This allows you to assess the relative location.
#' @param query A GenomicRanges object with query sets
#' @param feats A GenomicRanges object with features to test distance to

#' @export
featureDistribution = function(query, feats) {
	precedeInd = precede(query, feats)
	followInd = follow(query, feats)
	preDist = -distance(query, feats[precedeInd])
	postDist = distance(query, feats[followInd])
	postHits = -preDist > postDist

	dists = preDist
	dists[postHits] = postDist[postHits]
	return(dists)
}

#' Plots a histogram of distances to genomic features

#' Given the results from \code{distanceToNearestSymmetrical}, plots a histogram of
#' distances surrounding the features of interest

#' @param dists Results from \code{distanceToNearestSymmetrical}
#' @export
plotFeatureDist = function(dists,
	plotTitle="Distribution relative to features") {
	divisions = c(-Inf, -1e6, -1e4, -1000, -100, 0, 100, 1000, 10000, 1e6, Inf)
	labelCuts(divisions)
	labels = labelCuts(divisions, collapse=" to ", infBins=TRUE)

	cuts = cut(dists, divisions, labels)
	df = as.data.frame(table(cuts))
	g = ggplot(df, aes(x=cuts, y=Freq)) + 
		geom_bar(stat="identity") + 
		geom_vline(xintercept = length(divisions)/2, color="darkgreen") +
		theme_classic() + 
		theme(aspect.ratio=1) + 
		theme_blank_facet_label() + 
		xlab("Distance to feature") +
		ylab("Number of regions") +
		theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust=0.5)) + # vlab()
		theme(plot.title = element_text(hjust = 0.5)) + # Center title
		ggtitle(plotTitle)
	return(g)
}

