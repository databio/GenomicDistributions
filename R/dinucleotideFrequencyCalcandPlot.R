## dependencies
library("BSgenome")
library("BSgenome.Hsapiens.UCSC.hg19")
library("GenomicDistributions")
library("Biostrings")
genome <- BSgenome.Hsapiens.UCSC.hg19
bsg = GenomicDistributions::loadBSgenome('hg19')
queryFile = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
query = rtracklayer::import(queryFile)
## end dependencies

## calc function

calcDinuclFreq <- function(query, ref) {
  
  # here I used the same function structure from calcGCcontent
  # query and ref as input, gvec as output
  
  # if class(query) != "GRanges" & != "GRangesList" {
  #   stop("Error: query must be GRanges or GrangesList object")
  # else 
  #   # Recurse over each GRanges object
  #   x = lapply(query, calcdinuclfreq, ref)
  #   
  # }
  
  GenomicDistributions:::.validateInputs(list(query=c("GRanges","GRangesList"),
                                              ref="BSgenome"))
  
  # I believe this function validates the input as GRanges
  
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    x = lapply(query, calcdinuclfreq, ref)
    namelist=names(query)
    if (is.null(namelist)){
      newnames=seq_along(query)
      namelist=newnames
      
      names(x)=namelist
    }
    return(x)
  }
  
  
  # This returns a list of the inputs with lapply
  seqlevels(query, pruning.mode="coarse")=seqlevels(ref)
  v = IRanges::Views(ref, query)
  
  # The Views virtual class is a general container for storing a set of views 
  # on an arbitrary Vector object, called the "subject"
  # The subject in this case are inputs
  dnvec= as.data.frame(Biostrings::dinucleotideFrequency(v, simplify.as="collapsed", as.matrix = FALSE))
  
  
  dnvec$ID<-rownames(dnvec)
  colnames(dnvec)<-c("frequency", "dinucleotide")
  library(tidyverse)
  dnvec <- tibble::rowid_to_column(dnvec, "ID")
  dnvec<-dnvec[,-1]
  dnvec[,2]<-as.factor(dnvec[,2])
  dnvec[,1]<-as.numeric(dnvec[,1])
  
  ## adding columns for plotting
  # dnvec$x = dnvec$dinucleotide
  # dnvec$x = factor(dnvec$x, levels=dnvec$x)
  library(ggplot2)
  # I swapped out alphabetfrequency function for dinucleotideFrequency
  
  return(dnvec)
  
}
## end calc function

## plot function

plotDinuclContent = function(gcvectors) {
  
  GenomicDistributions:::.validateInputs(list(gcvectors=c("matrix", "array", "data.frame")))
  library(ggplot2)
  g = gcvectors
  plot=ggplot(data=g, aes(x=dinucleotide, y=frequency)) + ggtitle("Dinucleotide Frequency Plot")+
    geom_bar(stat="identity", colour="black", show.legend=FALSE) +
    
    
    geom_text(aes(label=ifelse(frequency>700, paste0(frequency, " (",sprintf("%1.1f", frequency/sum(frequency)*100),"%)"),
                               ifelse(frequency>300, frequency, "")), y=0.5*frequency), colour="white", size=2.5)+coord_flip() + theme_bw() + labs(x="")  
  return(plot)
}

## end plot function


## test functions

DNF<-calcDinuclFreq(query, bsg)
plotDinuclContent(DNF)

## end test functions