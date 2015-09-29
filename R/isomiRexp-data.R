#' @name isomiRexp
#' @title IsomirDataSeq
#' @description This data set is the object return by \code{\link{IsomirDataSeqFromFiles}}. It mainly
#' contains miRNA sequencing raw data from 6 samples: 3 newborns and 3 elderly individuals.
#' human individuals (Somel et al, 2010). 
#' Use \code{colData} to see the experiment design of the data.
#' @docType data
#' @aliases isomiRexp
#' @usage 
#' data("isomiRexp")
#' @source
#' Data is available from GEO dataset under accession number GSE18069. 
#' Samples used are: GSM450597, GSM450598, GSM450609 and GSM450608.
#' 
#' Every sample was analyzed with seqbuster tool, see 
#' \url{http://seqcluster.readthedocs.org/mirna_annotation.html} for
#' more details.
#' 
#' Adapter removal with the following parameters:
#' \code{adrec fastq_file TCGTATGCCGTCTT 8 1 0.34}
#' 
#' miRNAs annotation with the following parameters:
#' \code{miraligner adrec_output mirbase_files hsa 1 3 3 out_prefix}
#' 
#' The data was created with \code{isomiRs-package} package:
#' 
#' \code{library(isomiRs)}
#' 
#' \code{fns <- c("GSM450597.mirna", "GSM450598.mirna", 
#' "GSM450609.mirna", "GSM450608.mirna")}
#' 
#' \code{design <- data.frame(condition=c('new', 'new', 'new', 
#' 'old', 'old', 'old'))}
#' 
#' \code{isd <- IsomirDataSeqFromFiles(fns, design)}
#' 
#' @references 
#' Mehmet Somel et al. MicroRNA, mRNA, and protein expression link 
#' development and aging in human and macaque brain.
#' \code{Genome Research,20(9):1207-1218,2010.doi:10.1101/gr.106849.110}
#' @format a \code{\link{IsomirDataSeq}} class.
#' @author Lorena Pantano, 2015-05-19
NULL