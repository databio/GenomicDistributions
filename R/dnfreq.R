#' Calculate Dinuclotide content over genomic ranges
#' 
#' Given a reference genome (BSgenome object) and ranges on the
#' reference, this function returns a data.table with 
#' counts of dinucleotides within the GRanges object.
#' 
#' @param query A GRanges object with query sets
#' @param ref Reference genome BSgenome object
#' @return A data.table with counts of dinucleotides across the GRanges object
#' @export
#' @examples
#' \dontrun{ 
#' bsg = loadBSgenome('hg19')
#' DNF = calcDinuclFreq(vistaEnhancers, bsg)
#' }

calcDinuclFreq = function(query, ref) {
  
  .validateInputs(list(query=c("GRanges","GRangesList"),
                                              ref="BSgenome"))
   if (is(query, "GRangesList")) {
     
    x = lapply(query, calcdinuclfreq, ref)
    
    namelist=names(query)
    
    if (is.null(namelist)){
      
      newnames=seq_along(query)
      
      namelist=newnames
      
      names(x)=namelist
    }
    return(x)
  }
  
  seqlevels(query, pruning.mode="coarse")=seqlevels(ref)
  
  v = IRanges::Views(ref, query)
  
  dnvec= as.data.table(Biostrings::dinucleotideFrequency(v))
                       
  return(dnvec)
 }



#' Plot dinuclotide content over genomic ranges
#' 
#' Given \code{calcDinuclFreq} results, this function 
#' generates a violin plot of dinucleotide frequency
#' 
#' @param DNF A data.table of dinucleotide counts 
#' @return A ggplot object plotting distribution of dinucleotide content per dinucleotide
#' @export
#' @examples
#' 
#' diNucDT = data.table(GC = rnorm(400, mean=0.5, sd=0.1), 
#' CG = rnorm(400, mean=0.5, sd=0.5), 
#' AT = rnorm(400, mean=0.5, sd=1), 
#' TA = rnorm(400, mean=0.5, sd=1.5))
#' 
#' DNFPlot =  plotDinuclFreq(diNucDT)
#' 
plotDinuclFreq = function(DNFDataTable) {
  
  .validateInputs(list(DNFDataTable=c("matrix", "array", 
                                      "data.frame", "data.table")))
  library(ggplot2)
  
  g = DNFDataTable
  
  ## for violinplot
  
  g=reshape2::melt(g) 
  
  names(g)[names(g)=="variable"]="dinucleotide"
  
  names(g)[names(g)=="value"]="frequency"
  
  g$frequency=as.numeric(g$frequency)
  
  g$dinucleotide=as.character(g$dinucleotide)
  
  plot=ggplot(data=g, aes(dinucleotide, frequency)) + 
    geom_violin(scale="width", trim=TRUE) + 
    geom_boxplot(width=0.1, color="grey", alpha=0.2) + 
    coord_flip() + ggtitle("Dinucleotide Frequency") 
  
  return(plot)
}

#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome (BSgenome object) and ranges on the
#' reference, this function returns a data.table with 
#' counts of dinucleotides within the GRanges object.
#' 
#' @param query A GRanges object with query sets
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes with
#'     \code{getTSSs}.
#' @return A numeric vector or list of vectors with the GC percentage of 
#'     the query regions.
#' @export
#' @examples
#' \dontrun{
#' refAssembly = 'hg19'
#' GCcontent = calcGCContentRef(vistaEnhancers, refAssembly)
#' } 

calcDinuclFreqRef= function(query, ref) {
  
  .validateInputs(list(query=c("GRanges","GRangesList"),
                       
                       refAssembly="character"))
  
  ref = loadBSgenome(refAssembly)
  
  return(calcDinuclFreqRef(query, ref))
}
