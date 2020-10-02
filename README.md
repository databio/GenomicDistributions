[![Build Status](https://travis-ci.org/databio/GenomicDistributions.svg?branch=dev)](https://travis-ci.org/databio/GenomicDistributions)

# GenomicDistributions

Functions for calculating and plotting the distribution of query features (*e.g.* genomic ranges) across the genome. If you have a set of genomic ranges, the GenomicDistributions R package can help you with visualizations and comparison.

## Installing

### Main package

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicDistributions")
```

### Data package

[GenomicDistributuionsData](https://github.com/databio/GenomicDistributionsData): includes full data files, too large to include in GenomicDistributions


## Quick start

See the vignettes for more information: http://code.databio.org/GenomicDistributions

## Building long vignettes

In the [long_vignettes](/long_vignettes) are vignettes that require large external data and take a long time to run. Therefore, they should be pre-built. You can render them manually by running [long_vignettes/render-long-vignettes.R](long_vignettes/render-long-vignettes.R). This will use `knitr` to run the vignette and put the result into the `vignettes` folder, along with output figures.
