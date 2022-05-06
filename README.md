| Master | Dev |
|--------|-----|
|[![Build Status](https://travis-ci.org/databio/GenomicDistributions.svg?branch=master)](https://travis-ci.org/databio/GenomicDistributions) | [![Build Status](https://travis-ci.org/databio/GenomicDistributions.svg?branch=dev)](https://travis-ci.org/databio/GenomicDistributions) |



# GenomicDistributions

An R package that provides functions for 1) calculating and 2) visualizing a variety of statistics for a collection of genomic ranges. If you have a set of genomic ranges, such as a BED file the GenomicDistributions R package can help you to explore, annotate, visualize,and compare it.

## Installing

### Main package

With Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("GenomicDistributions")
```

Or from GitHub:

```r
devtools::install_github("databio/GenomicDistributions")
```

### Data package

[GenomicDistributionsData](https://github.com/databio/GenomicDistributionsData): includes full data files, too large to include in GenomicDistributions


## Quick start

See the vignettes for more information: http://code.databio.org/GenomicDistributions

## Building long vignettes

In the [long_vignettes](/long_vignettes) are vignettes that require large external data and take a long time to run. Therefore, they should be pre-built. You can render them manually by running [long_vignettes/render-long-vignettes.R](long_vignettes/render-long-vignettes.R). This will use `knitr` to run the vignette and put the result into the `vignettes` folder, along with output figures.

**Cite GenomicDistributions:**

Kupkova, K., Mosquera, J.V., Smith, J.P., Stolarczyk M, Danehy T., Lawson J.T., Rogers S., LeRoy N., Sheffield N.C. GenomicDistributions: fast analysis of genomic intervals with Bioconductor. *BMC Genomics* 23, 299 (2022). https://doi.org/10.1186/s12864-022-08467-y
