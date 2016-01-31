isomiRs
=======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45382.svg)](http://dx.doi.org/10.5281/zenodo.45382)

analyze isomiRs from seqbuster tool 

It needs R > 3.2

###Installation

Inside R:

```
library(devtools)
library(roxygen2)
devtools::install_github("lpantano/isomiRs", ref="master")
library(isomiRs)
openVignette(isomiRs)
```

There is a R script and set of example data at folder `$PATH2SEQBUSTER/R/isomiR_package/test`

Use last version with (only for BioC develop or R > 2.1):

```
devtools::install_github("lpantano/isomiRs", ref="develop")
library(isomiRs)
openVignette(isomiRs)
```

or look at the vignette here: http://lpantano.github.io/isomiRs/vignettes/isomiR-intro.pdf
