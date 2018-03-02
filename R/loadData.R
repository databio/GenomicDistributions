
#' Loads BSgenome objects from UCSC-style character vectors.
#'
#' This function will let you use a simple character vector (e.g. 'hg19') to
#' load and then return BSgenome objects. This lets you avoid having to use the
#' more complex annotation for a complete BSgenome object (e.g.
#' BSgenome.Hsapiens.UCSC.hg38.masked)
#' 
#' @param genomeBuild	One of 'hg19', 'hg38', 'mm10', 'mm9', or 'grch38'
#' @export
#' @examples
#' bsg = loadBSgenome('hg19')
# TODO: allow for grabbing either masked or unmasked BSgenome objects
loadBSgenome = function(genomeBuild) {
	# Convert the given string into the BSgenome notation
	if (!requireNamespace("BSgenome", quietly=TRUE)) {
		message("BSgenome package is not installed.")
	}
	databasePkgString = switch (genomeBuild,
		grch38 = "BSgenome.Hsapiens.UCSC.hg38.masked",
		hg38 = "BSgenome.Hsapiens.UCSC.hg38.masked",
		hg19 = "BSgenome.Hsapiens.UCSC.hg19.masked",
		mm10 = "BSgenome.Mmusculus.UCSC.mm10.masked",
		mm9 = "BSgenome.Mmusculus.UCSC.mm9.masked",
		bogus = "bogus" # a bogus (uninstalled) genome for unit tests
	)

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
#' ensDb = loadEnsDb('hg19')
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
			" to a BSgenome")
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

buildTSSs = function(assemblyList = list("hg38", "hg19", "mm10")) {
	for (refAssembly in assemblyList) {
		message(refAssembly)
		TSSGenomeVar = paste0("TSS_", refAssembly)
		EnsDb = loadEnsDb(refAssembly)
		featsWide = ensembldb::genes(EnsDb, columns=c("gene_biotype"))

		# Now, restrict to protein-coding genes and grab just a single base pair at the TSS
		feats = promoters(featsWide[featsWide$gene_biotype == "protein_coding"], 1, 1)
		# Change from ensembl-style chrom annotation to UCSC_style
		seqlevels(feats) = paste0("chr", seqlevels(feats))
		
		assign(TSSGenomeVar, feats, envir=environment())
		save(file=paste0(TSSGenomeVar, ".RData"), list=TSSGenomeVar)
	}
}

getChromSizes = function(refAssembly) {
	getReferenceData(refAssembly, tagline="chromSizes_")
}

getTSSs = function(refAssembly) { 
	getReferenceData(refAssembly, tagline="TSS_")
}

getReferenceData = function(refAssembly, tagline) {
	# query available datasets
	ad = data(package="GenomicDistributions")
	adm = ad$results[,"Item"]
	chromSizesGenomeVar = paste0(tagline, refAssembly)
	if (chromSizesGenomeVar %in% adm){
		# load it!
		data(list=chromSizesGenomeVar,
				package="GenomicDistributions",
				envir=environment())
		return(get(chromSizesGenomeVar))
	} else {
		message("I don't have data for reference assembly ",
			refAssembly)
	}
}





