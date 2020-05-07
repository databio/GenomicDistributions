#' hg19 chromosome sizes
#'
#' A dataset containing chromosome sizes for Homo Sapiens hg38 genome assembly
#'
#' @format A named vectors of lengths with one item per chromosome
#' @source BSgenome.Hsapiens.UCSC.hg19 package
#' @name chromSizes_hg19
#' @docType data
#' @keywords datasets
#' @usage data(chromSizes_hg19)
NULL


#' hg19 TSS locations
#'
#' A dataset containing chromosome sizes for Homo Sapiens hg38 genome assembly
#'
#' @format A named vectors of lengths with one item per chromosome
#' @source EnsDb.Hsapiens.v75 package
#' @name TSS_hg19
#' @docType data
#' @keywords datasets
#' @usage data(TSS_hg19)
NULL

#' hg38 gene models
#'
#' A dataset containing gene models for Homo Sapiens hg38 genome assembly. 
#'
#' @format A list of two GRanges objects, with genes and exons locations
#' @source EnsDb.Hsapiens.v75 package
#' @name geneModels_hg19
#' @docType data
#' @keywords datasets
#' @usage data(geneModels_hg19)
NULL


#’ Example hg19 open signal matrix 
#' 
#' A dataset containing a subset of open chromatin regions across all cell types defined by ENCODE for Homo Sapiens hg19
#'
#' Preparation steps:
#' \enumerate{
#'    \item{made a universe of regions by merging regions across cell types defined as opened in ENCODE}
#'    \item{took bigwig files from ENCODE for individual cell types, merged replicates, filtered out blacklisted sites}
#'    \item{evaluated the signal above regions defined by previous step}
#'    \item{performed quantile normalization}
#'    \item{subsetted it}
#' }
#'
#' @format data.frame, rows represent whole selection of open 
#' chromatin regions across all cell types defined by ENCODE, columns are 
#' individual cell types and values are normalized open chromatin signal values.
#' @source \url{http://big.databio.org/open_chromatin_matrix/openSignalMatrix_hg19_quantileNormalized_round4.txt.gz}
#' @name exampleOpenSignalMatrix_hg19
#' @docType data
#' @keywords datasets
#' @usage data(exampleOpenSignalMatrix_hg19)
NULL


#’ Example BED file
#' 
#' Example BED file read with rtracklayer::import
#'
#' @format GenomicRanges::GRanges
#' @name vistaEnhancers
#' @docType data
#' @keywords datasets
#' @usage data(vistaEnhancers)
NULL


#’ Example BED file
#' 
#' Example BED file read with rtracklayer::import
#'
#' @format GenomicRanges::GRanges
#' @name setB_100
#' @docType data
#' @keywords datasets
#' @usage data(setB_100)
NULL


#’ Cell type metadata matrix
#' 
#' Table the maps cell types to tissues and groups
#'
#' @format data.table with 3 columns (cellType, tissue and group) and 74 rows (one per cellType)
#' @source self-curated dataset
#' @name cellTypeMetadata
#' @docType data
#' @keywords datasets
#' @usage data(cellTypeMetadata)
NULL

