
#' Loads BSgenome objects from UCSC-style character vectors.
#'
#' This function will let you use a simple character vector (e.g. 'hg19') to
#' load and then return BSgenome objects. This lets you avoid having to use the
#' more complex annotation for a complete BSgenome object (e.g.
#' BSgenome.Hsapiens.UCSC.hg38.masked)
#' 
#' @param genomeBuild	One of 'hg19', 'hg38', 'mm10', 'mm9', or 'grch38'
#' @param masked	Should we used the masked version? Default:TRUE
#' @export
#' @examples
#' \dontrun{
#' bsg = loadBSgenome('hg19')
#' }
# TODO: allow for grabbing either masked or unmasked BSgenome objects
loadBSgenome = function(genomeBuild, masked=TRUE) {
	# Convert the given string into the BSgenome notation
	if (!requireNamespace("BSgenome", quietly=TRUE)) {
		message("BSgenome package is not installed.")
	}
	databasePkgString = switch (genomeBuild,
		grch38 = "BSgenome.Hsapiens.UCSC.hg38",
		hg38 = "BSgenome.Hsapiens.UCSC.hg38",
		hg19 = "BSgenome.Hsapiens.UCSC.hg19",
		mm10 = "BSgenome.Mmusculus.UCSC.mm10",
		mm9 = "BSgenome.Mmusculus.UCSC.mm9",
		bogus = "bogus" # a bogus (uninstalled) genome for unit tests
	)
	if (masked) {
		databasePkgString = paste0(databasePkgString, ".masked")
	}

	if (is.null(databasePkgString)) {
		stop("I don't know how to map the string ", genomeBuild,
			" to a BSgenome")
	}

	return(.requireAndReturn(databasePkgString))
}


#' Loads Ensemble gene annotations (EnsDb objects) from UCSC-style character vectors
#' 
#' This function will let you use a simple character vector (e.g. 'hg19') to
#' load and then return BSgenome objects. This lets you avoid having to use the
#' more complex annotation for a complete BSgenome object (e.g.
#' BSgenome.Hsapiens.UCSC.hg38.masked)
#' 
#' @param genomeBuild	One of 'hg19', 'hg38', 'mm10', 'mm9', or 'grch38'
#' @export
#' @examples
#' \dontrun{
#' ensDb = loadEnsDb('hg19')
#' }
loadEnsDb = function(genomeBuild) {

	databasePkgString = switch (genomeBuild,
		grch38 = "EnsDb.Hsapiens.v86",
		hg38 = "EnsDb.Hsapiens.v86",
		hg19 = "EnsDb.Hsapiens.v75",
		mm10 = "EnsDb.Mmusculus.v79",
		bogus = "bogus" # a bogus (uninstalled) db for unit tests
	)

	if (is.null(databasePkgString)) {
		stop("I don't know how to map the string ", genomeBuild,
			" to a EnsDb")
	}

	return(.requireAndReturn(databasePkgString))
}



# Generate chromSizes for common genome assemblies
buildChromSizes = function(assemblyList = list("hg38", "hg19", "mm10", "mm9")) {
	for (refAssembly in assemblyList) {
		message(refAssembly)
		chromSizesGenomeVar = paste0("chromSizes_", refAssembly)
		BSG = loadBSgenome(refAssembly)
		chromSizesGenome = seqlengths(BSG)
		assign(chromSizesGenomeVar, chromSizesGenome, envir=environment())
		save(file=paste0(chromSizesGenomeVar, ".RData"), list=chromSizesGenomeVar)
	}
}


# Generates Rdata objects for TSSs, to be included in the package for example.
buildTSSs = function(assemblyList = list("hg38", "hg19", "mm10", "mm9")) {
	if (!requireNamespace("ensembldb", quietly=TRUE)) {
		message("ensembldb package is not installed.")
	}
	for (refAssembly in assemblyList) {
		message(refAssembly)
		TSSGenomeVar = paste0("TSS_", refAssembly)
		# Using EnsDb
		EnsDb = NULL
		tryCatch( {
			EnsDb = loadEnsDb(refAssembly)
			codingFilter = AnnotationFilter::AnnotationFilter(
				~ gene_biotype == "protein_coding")
			featsWide = ensembldb::genes(EnsDb, filter=codingFilter)

			# Grab just a single base pair at the TSS
			feats = promoters(featsWide, 1, 1)
			# Change from ensembl-style chrom annotation to UCSC_style
			seqlevels(feats) = paste0("chr", seqlevels(feats))
		}, error=function(err) NA)

		if (is.null(EnsDb)) {
			# Using TxDb
			txdb = loadTxDb(refAssembly)
			# Grab just a single base pair at the TSS
			feats = GenomicFeatures::promoters(txdb, 1, 1)
		}
		assign(TSSGenomeVar, feats, envir=environment())
		save(file=paste0(TSSGenomeVar, ".RData"), list=TSSGenomeVar)
	}
}


loadTxDb = function(genomeBuild) {
	databasePkgString = switch (genomeBuild,
		grch38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
		hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
		hg19 = "TxDb.Hsapiens.UCSC.hg19.knownGene",
		mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
		mm9 = "TxDb.Mmusculus.UCSC.mm9.knownGene",
		bogus = "bogus" # a bogus (uninstalled) db for unit tests
	)

	if (is.null(databasePkgString)) {
		stop("I don't know how to map the string ", genomeBuild,
			" to a TxDb")
	}

	return(.requireAndReturn(databasePkgString))
}




# Generates RData objects for gene models, which can then be included in the 
# package
buildGeneModels = function(assemblyList = list("hg38", "hg19", "mm10", "mm9")) {

	if (!requireNamespace("ensembldb", quietly=TRUE)) {
		message("ensembldb package is not installed.")
	}
	for (refAssembly in assemblyList) {
		message(refAssembly)
		tagline = "geneModels_"
		storedObjectName = paste0(tagline, refAssembly)

		tryCatch( { 
			EnsDb = loadEnsDb(refAssembly)
			codingFilter = AnnotationFilter::AnnotationFilter(
				~ gene_biotype == "protein_coding")
			geneFeats = ensembldb::genes(EnsDb, filter = codingFilter, columns=NULL)
			exonFeats = ensembldb::exons(EnsDb, filter = codingFilter, columns=NULL)

			# Smash 
			exonFeats = reduce(exonFeats)
			# Since we're storing this data, we want it to be small.
			elementMetadata(geneFeats) = NULL
			elementMetadata(exonFeats) = NULL
			# Change from ensembl-style chrom annotation to UCSC_style
			seqlevels(geneFeats) = paste0("chr", seqlevels(geneFeats))
			seqlevels(exonFeats) = paste0("chr", seqlevels(exonFeats))
		}, error=function(err) NA)
		if (is.null(EnsDb)) {
			# Try a TxDb instead
			txdb = loadTxDb(refAssembly)
			exonFeats = GenomicFeatures::exons(txdb)
			geneFeats = GenomicFeatures::transcripts(txdb)
			geneFeats = reduce(geneFeats)
			exonFeats = reduce(exonFeats)
		}

		geneModels = list(genesGR=geneFeats, exonsGR=exonFeats)
		assign(storedObjectName, geneModels, envir=environment())
		save(file=paste0(storedObjectName, ".RData"), list=storedObjectName)
	}
}

getChromSizes = function(refAssembly) {
	getReferenceData(refAssembly, tagline="chromSizes_")
}

getTSSs = function(refAssembly) { 
	getReferenceData(refAssembly, tagline="TSS_")
}

getGeneModels = function(refAssembly) { 
	getReferenceData(refAssembly, tagline="geneModels_")
}

# This is a generic getter function that will return a data object requested,
# if it is included in the built-in data with the package. Data objects can 
# be requested for different reference assemblies and data types (specified by a
# tagline, which is a unique string identifying the data type).
# @refAssembly Reference assembly string (like hg38)
# @tagline the string that was used to identify data of a given type in the 
# data building step. It's used for the filename so we know what to load, and is 
# what makes this function generic (so it can load different data types).
getReferenceData = function(refAssembly, tagline) {
	# query available datasets and convert the packageIQR object into a vector
	datasetListIQR = utils::data(package="GenomicDistributions")
	datasetList = datasetListIQR$results[,"Item"]
	dataObjectVar = paste0(tagline, refAssembly)
	if (dataObjectVar %in% datasetList){
		# load it!
		utils::data(list=dataObjectVar,
				package="GenomicDistributions",
				envir=environment())
		return(get(dataObjectVar))
	} else {
		error("I don't have built-in data for reference assembly ",
			refAssembly, "(looking for ", tagline, ")")
	}
}

