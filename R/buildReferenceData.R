#' Read local or remote file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where
#'       the downloaded GTF file should be stored. If not provided, 
#'       a temporary directory will be used.
#'
#' @return data.frame retrieved file path
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' CElegansGtf = retrieveFile(CElegansGtfUrl)
retrieveFile <- function(source, destDir=NULL){
  if (is.null(destDir)) destDir = tempdir()
    # download file, if not local
  if (!file.exists(source)) {
    destFile = paste(destDir, basename(source), sep = "/")
    if (file.exists(destFile)){
      message("File exists: ", destFile)
    }else{
      message("File will be saved in: ", destFile)
      download.file(url = source, destfile = destFile)    
    }
  }else{
    destFile = source
    message("Got local file: ", destFile)
  }
  
    return(destFile)
}


#' Get transcription start sites (TSSs) from a remote or local GTF file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where 
#'        the downloaded GTF file should be stored
#' @param convertEnsemblUCSC a logical indicating whether Ensembl style 
#'        chromosome annotation should be changed to UCSC style
#'
#' @return a list of GRanges objects
#'
#' @import dplyr
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' CElegansTss = getTssFromGTF(CElegansGtfUrl, TRUE)
getTssFromGTF <- function(source, convertEnsemblUCSC=FALSE, destDir=NULL){
    GtfDf = as.data.frame(rtracklayer::import(retrieveFile(source, destDir)))
    subsetGtfDf = GtfDf %>% 
    dplyr::filter(gene_biotype == "protein_coding", type == "gene")
    gr = makeGRangesFromDataFrame(subsetGtfDf, keep.extra.columns = TRUE)
    feats = promoters(gr, 1, 1) 
    if(convertEnsemblUCSC)
      seqlevels(feats) = paste0("chr", seqlevels(feats))
    feats
}


#' Get gene models from a remote or local GTF file
#'
#' @param source a string that is either a path to a local or remote GTF
#' @param destDir a string that indicates the path to the directory where
#'        the downloaded GTF file should be stored
#' @param features a vector of strings with feature identifiers that to 
#'        include in the result list
#' @param convertEnsemblUCSC a logical indicating whether Ensembl style 
#'        chromosome annotation should be changed to UCSC style
#'
#' @return a list of GRanges objects
#'
#' @import dplyr
#' @export
#'
#' @examples
#' CElegansGtfUrl = "http://ftp.ensembl.org/pub/release-103/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.103.gtf.gz"
#' features = c("gene", "exon", "three_prime_utr", "five_prime_utr")
#' CElegansGeneModels = getGeneModelsFromGTF(CElegansGtfUrl, features, TRUE)
getGeneModelsFromGTF <- function(source,
                                 features,
                                 convertEnsemblUCSC = FALSE,
                                 destDir = NULL) {
  GtfDf = as.data.frame(rtracklayer::import(retrieveFile(source, destDir)))
  subsetGtfDf = GtfDf %>%
    filter(gene_biotype == "protein_coding")
  retList = list()
  message("Extracting features: ", paste(features, collapse = ", "))
  for (feat in features) {
    featGR = unique(GenomeInfoDb::keepStandardChromosomes(reduce(
      GenomicRanges::makeGRangesFromDataFrame(
        subsetGtfDf %>% filter(type == feat),
        keep.extra.columns = TRUE
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
#' @param source a string that is either a path to a  
#'        local or remote FASTA
#' @param destDir a string that indicates the path to the 
#'        directory where the downloaded FASTA file should be stored
#'
#' @return a named vector of sequence lengths
#' @importFrom Biostrings readDNAStringSet
#' @export
#'
#' @examples
#' CElegansUrl = "http://ftp.ensembl.org/pub/release-103/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz"
#' CElegansChromSizes = getChromSizesFromFasta(CElegansUrl)
getChromSizesFromFasta <- function(source, destDir=NULL) {
  fastaPath = retrieveFile(source, destDir)
  fastaStringSet = readDNAStringSet(fastaPath)
  oriNames = fastaStringSet@ranges@NAMES
  names = vapply(oriNames, function(x){
    strsplit(x, " ")[[1]][1]
  }, character(1))
  chromSizes = fastaStringSet@ranges@width
  names(chromSizes) = names
  chromSizes
}
