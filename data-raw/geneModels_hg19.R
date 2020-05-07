.requireAndReturn = function(pkg) {
    if (requireNamespace(pkg))
        return(utils::getAnywhere(pkg)$objs[[1]])
    else
        warning(pkg, " is not installed")
    return(NULL)
}

# Set reference-DB mappings
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


if (!requireNamespace("ensembldb", quietly=TRUE)) {
    message("ensembldb package is not installed.")
}

storedObjectName = "geneModels_hg19"
message(paste("building", storedObjectName, "using Ensembl"))
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
geneModels = list(genesGR=geneFeats, exonsGR=exonFeats)
assign(storedObjectName, geneModels)
do.call("use_data", list(as.name(storedObjectName), overwrite = TRUE))