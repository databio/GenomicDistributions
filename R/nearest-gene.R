#' Find overlapping annotations for a given query and annotation
#' set. Annotate the query based on the distance to the nearest
#' gene, the gene name, and the gene type.
#'
#' @param query A GRanges or GRangesList object with query sets
#' @param annotations A GRanges or GRangesList object with annotation sets
#' @param removeUnknowns Boolean value to indicate if you'd like to remove unannotaed ranges in the query.
#' 
#' @return A data table that contains observations for each genomic region
#'         and the associated aforementioned annotations.
#' @export
#' @examples
#' queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' query = rtracklayer::import(queryFile)
#' data(TSS_hg19)
#' 
#' queryAnnotated = calcNearestGenes(query, TSS_hg19)
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
  query$gene_id = "unknown"
  
  # map hits to gene type and the gene name
  query[queryHits(oLaps)]$gene_id = annotations[subjectHits(oLaps)]$gene_id
  
  # sort the query so that we can then run the analysis
  query = sort(query)
  query$nearest_gene = "unknown"
  query$nearest_gene_distance = 0
  
  # go through each chromosome
  for(chr in unique(seqnames(query))) {
    # extract regions for chromosome
    chrQuery = query[seqnames(query) == chr]
    
    # stuff back into query
    query[seqnames(query) == chr]$nearest_gene = chrQuery[nearest(chrQuery)]$gene_id
    query[seqnames(query) == chr]$nearest_gene_distance = as.data.frame(distanceToNearest(chrQuery))$distance
  }
  
  
  if(removeUnknowns) {
   query = query[query$nearest_gene != "unknown"] 
  }
  
  # dump a query to a data table and return
  dt = grToDt(query)
  
  # remove the gene_id column and return
  return(dt[,gene_id:=NULL])
}
