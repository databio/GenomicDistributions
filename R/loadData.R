
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
