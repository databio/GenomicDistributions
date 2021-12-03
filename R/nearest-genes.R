#' Given a query and set of annotations, this function will calculate
#' the nearest annotation to each region in the region set, as well
#' as the nearest gene type and the distance to the nearest gene.
#'
#' @param query A GRanges or GRangesList object with query sets
#' @param annotations A GRanges or GRangesList object with annotation sets
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
calcNearestGenes =  function(query, annotations, gene_name_key="gene_id", gene_type_key="gene_biotype") {
  .validateInputs(list(query=c("GRanges","GRangesList")))
  
  # calculate the nearest annotations to given query
  nearestIds = nearest(query, annotations)
  
  # annotate neaest gene and type
  query$nearest_gene = annotations[nearestIds]$gene_id
  query$nearest_gene_type = annotations[nearestIds]$gene_biotype
  
  # annotate on the distance as well
  query$nearest_distance = distance(query, annotations[nearestIds])
  
  # dump a query to a data table and return
  dt = grToDt(query)
  
  return(dt)
}