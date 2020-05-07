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
    s
    return(.requireAndReturn(databasePkgString))
}

# Generate chromSizes for selected genome assemblies
storedObjectName = "chromSizes_hg19"
message(paste("building", storedObjectName))
BSG = loadBSgenome(refAssembly)
chromSizesGenome = seqlengths(BSG)
assign(storedObjectName, chromSizesGenome)
do.call("use_data", list(as.name(storedObjectName), overwrite = TRUE))
rm(chromSizesGenome, storedObjectName)

