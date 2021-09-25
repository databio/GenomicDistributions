#' Find overlapping annotations for a given query and annotation
#' set. Annotate the query based on the distance to the nearest
#' gene, the gene name, and the gene type.
#'
#' @param query A GRanges or GRangesList object with query sets
#' @param annotatations A GRanges or GRangesList object with annotation sets
#' @param removeUnknowns Boolean value to indicate if you'd like to remove unannotaed ranges in the query.
#' 
#' @return A data table that contains observations for each genomic region
#'         and the associated aforementioned annotations.
#' @export
#' @examples
#' queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' query = rtracklayer::import(queryFile)
#' annotations = rtracklayer::import("data/hg38.refGene.gtf")
#' 
#' queryAnnotated = calcNearestGenesRef(query, annotations)
calcNearestGenes = function(query, annotations, removeUnknowns=TRUE)
{
  .validateInputs(list(query=c("GRanges","GRangesList")))
  # find overlaps between the query and
  # the annotations and use to annotate
  # our query
  oLaps = findOverlaps(query, annotations)
  
  # init the type column
  # as unknown. It could be
  # possible that our query had
  # a GRange that wasn't present
  # in the annotation file
  query$type = "unknown"
  query$gene_name = "unknown"
  
  # map hits to gene type and the gene name
  query[queryHits(oLaps)]$gene_name = annotations[subjectHits(oLaps)]$gene_name
  
  # mapped as.character to convert from factor to string
  query[queryHits(oLaps)]$type = as.character(annotations[subjectHits(oLaps)]$type)
  
  # sort the query so that we can then run the analysis
  query = sort(query)
  query$ng = "unknown"
  query$ng_type = "unknown"
  query$ng_distance = "unknown"
  
  # go through each chromosome
  for(chr in unique(seqnames(query))) {
    # extract regions for chromosome
    chrQuery = query[seqnames(query) == chr]
    
    # stuff back into query
    query[seqnames(query) == chr]$ng = chrQuery[nearest(chrQuery)]$gene_name
    query[seqnames(query) == chr]$ng_type = chrQuery[nearest(chrQuery)]$type
    query[seqnames(query) == chr]$ng_distance = as.data.frame(distanceToNearest(chrQuery))$distance
  }
  
  query$ng_distance
  query$ng_type
  query$ng
  
  if(removeUnknowns) {
   query = query[query$ng != "unknown"] 
  }
  
  # dump a query to a data table and return
  return(grToDt(query))
}