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
  dnvec= as.data.frame(Biostrings::dinucleotideFrequency(v))
                       
  return(dnvec)
 }

## test function

DNF<-calcDinuclFreq(query, bsg)

## end calc function

## plot function

plotDinuclContent = function(gcvectors) {
  
  GenomicDistributions:::.validateInputs(list(gcvectors=c("matrix", "array", "data.frame")))
  library(ggplot2)
  g = gcvectors
  
  ## for violinplot
  
  g=reshape2::melt(g) 
  names(g)[names(g)=="variable"]<-"dinucleotide"
  names(g)[names(g)=="value"]<-"frequency"
  g$frequency<-as.numeric(g$frequency)
  g$dinucleotide<-as.character(g$dinucleotide)
  plot=ggplot(data=g, aes(dinucleotide, frequency)) + geom_violin(scale="width", trim=TRUE) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + coord_flip() + ggtitle("Dinucleotide Frequency") 
  
  return(plot)
}

# test

plotDinuclContent(DNF)

## end plot function 
 

### check max function




## try violin plot



###
v = IRanges::Views(bsg, query)
dnvec= as.data.frame(Biostrings::dinucleotideFrequency(v))
dnvec=reshape2::melt(dnvec)
library('ggplot2')
names(dnvec)[names(dnvec)=="variable"]<-"dinucleotide"
names(dnvec)[names(dnvec)=="value"]<-"frequency"
# dnvec$ID <- rownames(dnvec)


dnvec$frequency<-as.numeric(dnvec$frequency)
dnvec$dinucleotide<-as.character(dnvec$dinucleotide)

ggplot(dnvec, aes(dinucleotide, frequency)) +
  geom_violin(scale="width", trim=TRUE) + geom_boxplot(width=0.1, color="grey", alpha=0.2) + coord_flip() + ggtitle("Dinucleotide Frequency")



## check max  function output

# DNF$density
# ggtitle("Dinucleotide Frequency Plot")
# 
# geom_density()
# 
# 
# dotplot = function(sigList, violin=TRUE, nameList=NULL, title=NULL, ylab="", xlab="") {
#   if(is.null(nameList)) {
#     if (is.null(names(sigList))) {
#       nameList = 1:length(sigList) #fallback to sequential numbers
#     } else {
#       nameList=names(sigList) #use names if available
#     }
#   }
#   
#   sigListLengths = sapply(sigList, length)
#   if (sum(sigListLengths) > 1e4) {
#     warning("You have more than 10k points, consider using a density plot")
#   }
#   DT = data.table(name=factor(rep(nameList, sigListLengths), levels=nameList), value= unlist(sigList))
#   
#   g = ggplot(DT, aes(factor(name), value, fill=factor(name))) 
#   if (violin) { 
#     # Add a violin plot on bottom.
#     g = g + geom_violin(color="gray", width = .66, aes(fill=factor(name)),  alpha=.4) 
#   }
#   g = g + geom_boxplot(width=.33, aes(alpha=0), outlier.shape = NA) + geom_point(position="jitter", aes(fill=factor(name)), alpha=.5, color="black", pch=21) + scale_fill_brewer(palette="Set1") + scale_color_brewer(palette="Set1") +  theme_classic() + ylab(ylab) + xlab(xlab) + ggtitle(title)  + theme(aspect.ratio=1) + theme(legend.position="none") #+ stat_summary(fun.y="median", geom="point", pch=3) + geom_violin(aes(fill="gray", color="gray"))
#   
#   return(g)
# }

##### atherocluster

### build scatter plot 1 axis frip score other axis accessibility ever single peak 

## raster graphic/contour plot