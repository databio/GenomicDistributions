% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/content-plots.R
\name{plotDinuclFreq}
\alias{plotDinuclFreq}
\title{Plot dinuclotide content within region set(s)}
\usage{
plotDinuclFreq(DNFDataTable)
}
\arguments{
\item{DNFDataTable}{A data.table, data.frame, or a list of dinucleotide counts - 
results from \code{calcDinuclFreq} or \code{calcDinuclFreqRef}}
}
\value{
A ggplot object plotting distribution of dinucleotide content in query regions
}
\description{
Given \code{calcDinuclFreq} or \code{calcDinuclFreqRef} results, this function 
generates a violin plot of dinucleotide frequency
}
\examples{

DNFDataTable = data.table::data.table(GC = rnorm(400, mean=0.5, sd=0.1), 
CG = rnorm(400, mean=0.5, sd=0.5), 
AT = rnorm(400, mean=0.5, sd=1), 
TA = rnorm(400, mean=0.5, sd=1.5))
DNFPlot =  plotDinuclFreq(DNFDataTable)

\dontrun{
query = system.file("extdata", "vistaEnhancers.bed.gz", package="GenomicDistributions")
GRquery = rtracklayer::import(query)
refAssembly = 'hg19'
DNF = calcDinuclFreqRef(GRquery, refAssembly)
DNFPlot2 =  plotDinuclFreq(DNF)
} 
}
