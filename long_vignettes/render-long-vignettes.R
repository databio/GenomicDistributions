knitr::opts_knit$set(base.dir = 'vignettes/', progress = TRUE, verbose = TRUE, fig.path="figures-full-power/")
knitr::knit("long_vignettes/full-power.Rmd", "vignettes/full-power.Rmd")
knitr::opts_knit$set(base.dir = 'vignettes/', progress = TRUE, verbose = TRUE, fig.path="figures-GDData/")
knitr::knit("long_vignettes/GenomicDistributionsData.Rmd", "vignettes/GenomicDistributionsData.Rmd")
