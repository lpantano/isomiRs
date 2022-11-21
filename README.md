isomiRs
=======

<img src="https://raw.githubusercontent.com/lpantano/isomiRs/main/inst/stickers/isomirs.png" width="150" height="150" align="right"/>

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.45382.svg)](http://dx.doi.org/10.5281/zenodo.45382)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build status](https://secure.travis-ci.org/lpantano/isomiRs.png)](https://travis-ci.org/lpantano/isomiRs)
[![coverage](https://img.shields.io/codecov/c/github/lpantano/isomiRs/master.svg)](https://codecov.io/github/lpantano/isomiRs?branch=master)


Analyze isomiRs from seqbuster tool  or any BAM file after using [seqbuster miraligner](http://seqcluster.readthedocs.io/mirna_annotation.html#mirna-isomirs-annotation-with-python)


[Bioconductor]: https://bioconductor.org
[devtools]: https://cran.r-project.org/package=devtools
[R]: https://www.r-project.org

## Installation

This is an [R][] package.

### [Bioconductor][] stable version

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
install.packages("BiocManager")
BiocManager::install("isomiRs")
```

### [Bioconductor][] latest version

```r
devtools::install_git("https://git@git.bioconductor.org/packages/isomiRs")
```

### [devtools][] development version

```r
install.packages("devtools")
devtools::install_github("lpantano/isomiRs")
```
