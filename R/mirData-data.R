#' @name mirData
#' @title Example of IsomirDataSeq with human brain miRNA counts data
#' @description This data set is the object return by \code{\link{IsomirDataSeqFromFiles}}.
#' It contains miRNA count data from 14 samples: 7 control individuals (pc) and
#' 7 patients with Parkinson's disease in early stage (Pantano et al, 2016).
#' Use \code{colData} to see the experiment design.
#' @docType data
#' @aliases mirData
#' @usage
#' data("mirData")
#' @source
#' Data is available from GEO dataset under accession number GSE97285
#'
#' Every sample was analyzed with seqbuster tool, see
#' \url{http://seqcluster.readthedocs.org/mirna_annotation.html} for
#' more details. You can get same files running the small RNA-seq
#' pipeline from \url{https://github.com/bcbio/bcbio-nextgen}.
#'
#' bcbio_nextgen was used for the full analysis.
#'
#' See \code{raw-data.R} to know how to recreate the object. 
#' This script is inside "extra" folder of the package.
#'
#' @references
#' Pantano L, Friedlander MR, Escaramis G, Lizano E et al.
#' Specific small-RNA signatures in the amygdala at premotor and motor stages
#' of Parkinson's disease revealed by deep sequencing analysis.
#' Bioinformatics 2016 Mar 1;32(5):673-81. PMID: 26530722
#' @format a \code{\link{IsomirDataSeq}} class.
#' @author Lorena Pantano, 2018-04-27
NULL