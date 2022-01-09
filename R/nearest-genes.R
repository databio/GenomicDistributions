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
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    x = lapply(query, calcFeatureDist, features)
    return(x)
  }
  
  # calculate the nearest annotations to given query
  nearestIds = nearest(query, annotations)
  
  # annotate nearest gene and type
  nearestGenes = annotations[nearestIds]
  
  #
  # convert nearestGenes GRange object to data-table
  # and dynamically access the column that way...
  # this is used to circumvent the fact that we cannot
  # dynamically access metadata columns inside a GRanges
  # object like we can a datatable:
  #   col = "gene_id"
  #   dt[[col]]
  #   ^^^ This doesnt work in a GRanges object.
  #
  query$nearest_gene = grToDt(nearestGenes)[[gene_name_key]]
  query$nearest_gene_type = grToDt(nearestGenes)[[gene_type_key]]
  
  # annotate on the distance as well
  query$nearest_distance = distance(query, annotations[nearestIds])
  
  # dump a query to a data table and return
  dt = grToDt(query)
  
  return(dt)
}