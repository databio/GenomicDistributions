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
  query$ng_distance = 0
  
  # go through each chromosome
  for(chr in unique(seqnames(query))) {
    # extract regions for chromosome
    chrQuery = query[seqnames(query) == chr]
    
    # stuff back into query
    query[seqnames(query) == chr]$ng = chrQuery[nearest(chrQuery)]$gene_name
    query[seqnames(query) == chr]$ng_type = chrQuery[nearest(chrQuery)]$type
    query[seqnames(query) == chr]$ng_distance = as.data.frame(distanceToNearest(chrQuery))$distance
  }
  
  # convert the distances to log_10
  # values
  query$ng_distance = log10(query$ng_distance)
  
  if(removeUnknowns) {
   query = query[query$ng != "unknown"] 
  }
  
  # dump a query to a data table and return
  return(grToDt(query))
}

plotNearestGenes = function(df) {
  .validateInputs(list(df=c("data.frame")))
  
  g = ggplot2::ggplot(df, aes(x=ng_distance, group=ng_type, fill=ng_type)) +
    geom_density(alpha=0.4) + 
    theme_classic()
  
  g = g + 
    xlab(expression(log[10]*("bp distance by nearest gene type"))) +
    xlim(0, 10) +
    theme(aspect.ratio=1) +
    theme_blank_facet_label() +
    ggtitle("Neighboring regions distance distribution by gene type") +
    theme(plot.title = element_text(hjust=0.5))
  
  return(g)
}
