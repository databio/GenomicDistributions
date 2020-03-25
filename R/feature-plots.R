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

#' Plots a histogram of distances to genomic features
#' 
#' Given the results from \code{featureDistribution}, plots a histogram of
#' distances surrounding the features of interest
#' 
#' @param dists Results from \code{featureDistribution}
#' @param featureName Character vector for plot labels (optional).
#' @export
plotFeatureDist = function(dists, featureName="features", divisions=NULL) {

	df = cutDists(dists, divisions)
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
cutDists = function(dists, divisions=NULL) {
	if (is.null(divisions)) {
		poscuts = c(10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000, 1e5, 1e6, Inf)
		divisions = c(-poscuts, 0, poscuts)
	}
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

	labels = labelCuts(sort(divisions), collapse=" to ", infBins=TRUE)
	cuts = cut(dists, divisions, labels)
	df = as.data.frame(table(cuts))
	return(df)
}