#' @name mirData
#' @title Example of IsomirDataSeq with human brain miRNA counts data
#' @description This data set is the object return by \code{\link{IsomirDataSeqFromFiles}}.
#' It contains miRNA count data from 6 samples: 3 newborns and 3 elderly
#' human individuals (Somel et al, 2010).
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
#' pipeline from \url{https://github.com/chapmanb/bcbio-nextgen}.
#'
#' bcbio_nextgen was used for the full analysis.
#'
#' \code{library(isomiRs)}
#' \code{files = list.files(file.path(root_path),pattern = "mirbase-ready",
#' recursive = T,full.names = T)}
#' \code{metadata_fn =  list.files(file.path(root_path),
#' pattern = "summary.csv$",recursive = T, full.names = T)}
#' \code{metadata = read.csv(metadata_fn, row.names="sample_id")}
#' \code{condition = names(metadata)[1]}
#' \code{mirData <- IsomirDataSeqFromFiles(files[rownames(design)], metadata)}
#'
#' @references
#' Pantano L, Friedlander MR, Escaramis G, Lizano E et al.
#' Specific small-RNA signatures in the amygdala at premotor and motor stages
#' of Parkinson's disease revealed by deep sequencing analysis.
#' Bioinformatics 2016 Mar 1;32(5):673-81. PMID: 26530722
#' @format a \code{\link{IsomirDataSeq}} class.
#' @author Lorena Pantano, 2016-04-07
NULL