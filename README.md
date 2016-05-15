isomiRs
=======

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45382.svg)](http://dx.doi.org/10.5281/zenodo.45382)

analyze isomiRs from seqbuster tool  or any BAM file after using [seqcluster seqbuster](http://seqcluster.readthedocs.io/mirna_annotation.html#mirna-isomirs-annotation-with-python)

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

or look at the vignette here: http://lpantano.github.io/isomiRs/vignettes/isomiR-intro.pdf
