#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome as a BSgenome object and some ranges on that
#' reference, this function will return a vector of the same length as the
#' granges object, with percent of Cs and Gs.
#' 
#' @param query  A GenomicRanges or GenomicRangesList object with query regions.
#' @param ref Reference genome BSgenome object.
#' @return A numeric vector of list of vectors with the GC percentage of 
#'     the query regions.
#' @export
#' @examples
#' \dontrun{
#' bsg = loadBSgenome('hg19')
#' gcvec = calcGCContent(vistaEnhancers, bsg)
#' }
calcGCContent = function(query, ref) {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         ref="BSgenome"))
    if (is(query, "GRangesList")) {
        # Recurse over each GRanges object
        x = lapply(query, calcGCContent, ref)
        namelist = names(query)
        if (is.null(namelist)) {
            newnames = seq_along(query)
            namelist = newnames
            # Append names
            names(x) = namelist
        }
        return(x)
    }
    # Restrict the seqnames to known chromosomes
    query = GenomeInfoDb::keepStandardChromosomes(query, pruning.mode="coarse")
    v = IRanges::Views(ref, query)
    gcvec = apply(Biostrings::alphabetFrequency(v)[,c("C","G")],1, sum)/width(v)
    return(gcvec)
}


#' Calculate GC content over genomic ranges
#' 
#' Given a reference genome as a BSgenome object and some ranges on that
#' reference, this function will return a vector of the same length as the
#' granges object, with percent of Cs and Gs.
#' 
#' @param query A GenomicRanges or GenomicRangesList object with query regions
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
calcGCContentRef = function(query, refAssembly) {
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         refAssembly="character"))
    ref = loadBSgenome(refAssembly)
    return(calcGCContent(query, ref))
}

#' Plots a density distribution of GC vectors

#' Give results from the \code{calcGCContent} function, this will produce a
#' density plot
#' @param gcvectors A numeric vector or list of numeric vectors of GC contents.
#' @return A ggplot object plotting distribution of GC content in query regions.
#' @export
#' @examples
#' numVector = rnorm(400, mean=0.5, sd=0.1)
#' GCplot = plotGCContent(numVector)
#' vecs = list(example1 = rnorm(400, mean=0.5, sd=0.1), 
#'             example2 = rnorm(600, mean=0.5, sd=0.1))
#' GCplot = plotGCContent(vecs)
#' 
plotGCContent = function(gcvectors) {
    .validateInputs(list(gcvectors=c("numeric", "list")))
    
    if (is(gcvectors, "list")) {
        nameList = names(gcvectors)
        vectorLengths = unlist(lapply(gcvectors, length))
        gcdfReshaped = data.frame(value = unlist(gcvectors),
                                  regionSet = rep(nameList, vectorLengths))
        meansdf = aggregate(gcdfReshaped$value, 
                            list(gcdfReshaped$regionSet), mean)
        g = ggplot2::ggplot(gcdfReshaped, aes(x=value, colour=regionSet)) +
            geom_density() +
            geom_vline(data=meansdf, aes(xintercept=x, colour=Group.1),
                       linetype="dashed", size=0.5) +
            theme_classic() +
            theme(legend.position = "bottom")
    } else {
        # plot a single regionset
        gcdfReshaped = data.frame(value = gcvectors)
        g = ggplot2::ggplot(gcdfReshaped, aes(x=value)) + 
            geom_density() + 
            geom_vline(aes(xintercept=mean(value)),
                       color="red", linetype="dashed", size=0.5) + 
            theme_classic()
    }    
    g = g + 
        ggtitle("GC content distribution") + 
        theme(plot.title = element_text(hjust=0.5)) +
        xlab("GC content") + 
        xlim(0,1) 
    return(g)
}

#' Calculate Dinuclotide content over genomic ranges
#' 
#' Given a reference genome (BSgenome object) and ranges on the
#' reference, this function returns a data.table with 
#' counts of dinucleotides within the GRanges object.
#' 
#' @param query A GRanges object with query sets
#' @param ref Reference genome BSgenome object
#' @param rawCounts a logical indicating whether the raw numbers should be 
#'     displayed, rather than percentages (optional).
#' @return A data.table with counts of dinucleotides across the GRanges object
#' @export
#' @examples
#' \dontrun{ 
#' bsg = loadBSgenome('hg19')
#' DNF = calcDinuclFreq(vistaEnhancers, bsg)
#' }

calcDinuclFreq = function(query, ref, rawCounts=FALSE) {
    
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         ref="BSgenome"))
    if (is(query, "GRangesList")) {
        
        # Recurse over each GRanges object
        x = lapply(query, calcDinuclFreq, ref, rawCounts=rawCounts)
        
        # return a list of dinucleotide dataframes across each GRanges object
        return(x)
    }
    # Restrict the seqnames to known chromosomes
    query = GenomeInfoDb::keepStandardChromosomes(query, pruning.mode="coarse")
    v = IRanges::Views(ref, query)
    regionNames = data.frame(region = paste(seqnames(query), 
                                            start(query), 
                                            end(query), sep="_"))
    dnvec= Biostrings::dinucleotideFrequency(v)
    # claculate frequencies if raw counts not required
    if(!rawCounts){
      dnvec = prop.table(dnvec, margin = 1)*100
    }
    dnvec = cbind(regionNames, as.data.frame(dnvec))
    return(dnvec)
}


#' Calculate dinucleotide content over genomic ranges
#' 
#' Given a reference genome (BSgenome object) and ranges on the
#' reference, this function returns a data.table with 
#' counts of dinucleotides within the GRanges object.
#' 
#' @param query A GRanges object with query sets
#' @param refAssembly A character vector specifying the reference genome
#'     assembly (*e.g.* 'hg19'). This will be used to grab chromosome sizes with
#'     \code{getTSSs}.
#' @param rawCounts a logical indicating whether the raw numbers should be 
#'     displayed, rather than percentages (optional).
#' @return A numeric vector or list of vectors with the GC percentage of 
#'     the query regions.
#' @export
#' @examples
#' \dontrun{
#'query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#'GRquery = rtracklayer::import(query)
#'refAssembly = 'hg19'
#'DNF = calcDinuclFreqRef(GRquery, refAssembly)
#' } 

calcDinuclFreqRef= function(query, refAssembly, rawCounts=FALSE) {
    
    .validateInputs(list(query=c("GRanges","GRangesList"),
                         
                         refAssembly="character"))
    
    ref = loadBSgenome(refAssembly)
    
    return(calcDinuclFreq(query, ref, rawCounts=rawCounts))
}


#' Plot dinuclotide content within region set(s)
#' 
#' Given \code{calcDinuclFreq} or \code{calcDinuclFreqRef} results, this function 
#' generates a violin plot of dinucleotide frequency
#' 
#' @param DNFDataTable A data.table, data.frame, or a list of dinucleotide counts - 
#'                    results from \code{calcDinuclFreq} or \code{calcDinuclFreqRef}
#' @return A ggplot object plotting distribution of dinucleotide content in query regions
#' @export
#' @examples
#' 
#' DNFDataTable = data.table::data.table(GC = rnorm(400, mean=0.5, sd=0.1), 
#' CG = rnorm(400, mean=0.5, sd=0.5), 
#' AT = rnorm(400, mean=0.5, sd=1), 
#' TA = rnorm(400, mean=0.5, sd=1.5))
#' DNFPlot =  plotDinuclFreq(DNFDataTable)
#' 
#' \dontrun{
#' query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
#' GRquery = rtracklayer::import(query)
#' refAssembly = 'hg19'
#' DNF = calcDinuclFreqRef(GRquery, refAssembly)
#' DNFPlot2 =  plotDinuclFreq(DNF)
#' } 

plotDinuclFreq = function(DNFDataTable) {
  .validateInputs(list(DNFDataTable=c("data.table","data.frame","list")))
  
  # reshape the data for plotting
  if (is(DNFDataTable, "list") && 
      any(vapply(DNFDataTable, function(x) any(names(x) == "region"), logical(1)))){
    g = reshape2::melt(DNFDataTable,id.vars="region", 
                       variable.name="dinucleotide", value.name="frequency") 
  } else if ((is(DNFDataTable, "data.frame") | is(DNFDataTable, "data.table"))&& 
             ("region" %in% colnames(DNFDataTable))){
    g = reshape2::melt(DNFDataTable,id.vars="region", 
                       variable.name="dinucleotide", value.name="frequency") 
  } else {
    g = reshape2::melt(DNFDataTable, id.vars=NULL,
                       variable.name="dinucleotide", value.name="frequency") 
  }
  
  # plot data as violin plots
  # if multiple inuts - make a facet for each dinucleotide to make the plot easier to read
  if (is(DNFDataTable, "list")){
    plot = ggplot2::ggplot(data=g, ggplot2::aes(x=L1, y=frequency, fill=L1)) +
      facet_wrap(~dinucleotide, nrow=4) +
      theme_bw() +
      theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_text(face = "bold")) +
      theme(axis.text.x = element_text(angle=90, hjust=1)) +
      xlab(" ")
  } else{
    plot = ggplot2::ggplot(data=g, ggplot2::aes(x=dinucleotide, y=frequency))+
      xlab("Dinucleotide")+
      theme_bw()
  }
  plot = plot +
    geom_violin(trim=TRUE, scale = "width") +
    geom_boxplot(alpha=0.2, outlier.shape = NA)+
    ggtitle("Dinucleotide Frequency") +
    guides(fill="none") + 
    theme(plot.title = element_text(hjust = 0.5))
  # check if we have raw counts or frequencies
  if (is(g[,"frequency"], "integer")){
    plot = plot + 
      ylab("Dinucleotide counts per region [n]")
    
  } else {
    plot = plot + 
      ylab("Dinucleotide frequency per region [%]")
  }
  return(plot)
}
