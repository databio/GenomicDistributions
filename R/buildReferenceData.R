#' Read local or remote file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where the downloaded GTF file should be stored
#'
#' @return data.frame retrieved file path
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' CElegansGtf = retrieveFile(CElegansGtfUrl, getwd())
retrieveFile <- function(source, destDir){
    # download file, if not local
    if (!file.exists(source)) {
        destFile = paste(destDir, basename(source), sep = "/")
        message("File will be saved in: ", destFile)
        download.file(url = source, destfile = destFile)
    } else{
        destFile = source
        message("Got local file: ", destFile)
    }
    
    return(destFile)
}


#' Get transcription start sites (TSSs) from a remote or local GTF file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where the downloaded GTF file should be stored
#' @param convertEnsemblUCSC a logical indicating whether Ensembl style chromosome annotation should be changed to UCSC style
#'
#' @return a list of GRanges objects
#'
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' CElegansTss = getTssFromGTF(CElegansGtfUrl, getwd(), TRUE)
getTssFromGTF <- function(source, destDir, convertEnsemblUCSC){
    
    message("Reading GTF file: ", destFile)
    GtfDf = as.data.frame(rtracklayer::import(retrieveFile(source, destDir)))
    
    subsetGtfDf = GtfDf %>% 
        filter(gene_biotype == "protein_coding", type == "gene")
    gr = makeGRangesFromDataFrame(subsetGtfDf, keep.extra.columns = T)
    feats = promoters(gr, 1, 1)
    if(convertEnsemblUCSC)
        seqlevels(feats) = paste0("chr", seqlevels(feats))
    feats
}


#' Get gene models from a remote or local GTF file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where the downloaded GTF file should be stored
#' @param features a vector of strings with feature identifiers that to include in the result list
#' @param convertEnsemblUCSC a logical indicating whether Ensembl style chromosome annotation should be changed to UCSC style
#'
#' @return a list of GRanges objects
#'
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' features = c("gene", "exon", "three_prime_utr", "five_prime_utr")
#' CElegansGeneModels = getGeneModelsFromGTF(CElegansGtfUrl, getwd(), features, TRUE)
getGeneModelsFromGTF <- function(source,
                                 destDir,
                                 features,
                                 convertEnsemblUCSC = FALSE) {
    message("Reading GTF file: ", destFile)
    GtfDf = as.data.frame(rtracklayer::import(retrieveFile(source, destDir)))
    subsetGtfDf = GtfDf %>%
        filter(gene_biotype == "protein_coding")
    retList = list()
    message("Extracting features: ", paste(features, collapse = ", "))
    for (feat in features) {
        featGR = unique(keepStandardChromosomes(reduce(
            makeGRangesFromDataFrame(
                subsetGtfDf %>% filter(type == feat),
                keep.extra.columns = T
            )
        ),
        pruning.mode = "coarse"))
        # change from Ensembl style chromosome annotation to UCSC style
        if (convertEnsemblUCSC)
            seqlevels(featGR) =  paste0("chr", seqlevels(featGR))
        retList[[feat]] = featGR
    }
    retList
}


#' Get gene models from a remote or local FASTA file
#'
#' @param source a string that is either a path to a local or remote FASTA
#' @param destDir a string that indicates the path to the directory where the downloaded FASTA file should be stored
#'
#' @return a named vector of sequence lengths
#' @export
#'
#' @examples
#' CElegansUrl = "http://ftp.ensembl.org/pub/release-103/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
#' CElegansChromSizes = getChromSizesFromFasta(CElegansUrl, getwd())
getChromSizesFromFasta <- function(source, destDir) {
    fastaPath = retrieveFile(source, destDir)
    fastaStringSet = readDNAStringSet(fastaPath)
    oriNames = fastaStringSet@ranges@NAMES
    names = sapply(oriNames, function(x){
        strsplit(x, " ")[[1]][1]
    })
    chromSizes = fastaStringSet@ranges@width
    names(chromSizes) = names
    chromSizes
}