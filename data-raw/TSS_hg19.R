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
storedObjectName = "TSS_hg19"
# Using EnsDb
message(paste("building", storedObjectName, "using Ensembl"))
EnsDb = loadEnsDb(refAssembly)
codingFilter = AnnotationFilter::AnnotationFilter(
    ~ gene_biotype == "protein_coding")
featsWide = ensembldb::genes(EnsDb, filter=codingFilter)

# Grab just a single base pair at the TSS
feats = promoters(featsWide, 1, 1)
# Change from ensembl-style chrom annotation to UCSC_style
seqlevels(feats) = paste0("chr", seqlevels(feats))
assign(storedObjectName, feats)
do.call("use_data", list(as.name(storedObjectName), overwrite = TRUE))
rm(feats, storedObjectName)