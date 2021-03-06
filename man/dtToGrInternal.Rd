% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{dtToGrInternal}
\alias{dtToGrInternal}
\title{Two utility functions for converting data.tables into GRanges objects}
\usage{
dtToGrInternal(DT, chr, start, end = NA, strand = NA, name = NA, metaCols = NA)
}
\arguments{
\item{DT}{A data.table representing genomic regions.}

\item{chr}{A string representing the chromosome column.}

\item{start}{A string representing the name of the start column.}

\item{end}{A string representing the name of the end column.}

\item{strand}{A string representing the name of the strand column.}

\item{name}{A string representing the name of the name column.}

\item{metaCols}{A string representing the name of the metadata column(s)
to include in the returned GRanges object.}
}
\value{
A GRanges object.
}
\description{
Two utility functions for converting data.tables into GRanges objects
}
