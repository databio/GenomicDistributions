[![Build Status](https://travis-ci.org/databio/GenomicDistributions.svg?branch=dev)](https://travis-ci.org/databio/GenomicDistributions)

# GenomicDistributions

Functions for calculating and plotting the distribution of query features (*e.g.* genomic ranges) across the genome. If you have a set of genomic ranges, the GenomicDistributions R package can help you with some simple visualizations.

## Installing

### Main package

```r
devtools::install_github("databio/GenomicDistributions")
```

### Data package

Includes full data files, too large to include in GenomicDistributions

```r
install.packages("http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.1.tar.gz", repos=NULL)
```

## Quick start

See the vignettes for more information: http://code.databio.org/GenomicDistributions

## Building long vignettes

In the [long_vignettes](/long_vignettes) are vignettes that require large external data and take a long time to run. Therefore, they should be pre-built. You can render them manually by running [long_vignettes/render-long-vignettes.R](long_vignettes/render-long-vignettes.R). This will use `knitr` to run the vignette and put the result into the `vignettes` folder, along with output figures.
