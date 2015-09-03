#' @name isomiRexp
#' @title IsomirDataSeq
#' @description This data set is the object return by \code{IsomirDataSeqFromFiles}. It mainly
#' contains miRNA sequencing raw data from 6 samples: 3 newnborns and 3 elderly 
#' human individuals (Somel et al, 2010). 
#' Use \code{colData} to see the experiment design of the data.
#' @docType data
#' @aliases isomiRexp
#' @usage 
#' data("isomiRexp")
#' @source
#' Data is available from GEO dataset under accesion number GSE18069. 
#' Samples used are: GSM450597, GSM450598, GSM450609 and GSM450608.
#' 
#' Sequence data was analyzed with seqbuster tool, see 
#' \url{http://seqcluster.readthedocs.org/mirna_annotation.html} for
#' more details.
#' 
#' Adapter removal with the following parameters:
#' \code{adrec fastq_file TCGTATGCCGTCTT 8 1 0.34}
#' 
#' miRNAs annotation with the following parameters:
#' \code{miraligner adrec_output mirbase_files hsa 1 3 3 out_prefix}
#' 
#' @references 
#' Mehmet Somel et al. MicroRNA, mRNA, and protein expression link 
#' development and aging in human and macaque brain. 
#' Genome Research, 20(9):1207–1218, 2010.doi:10.1101/gr.106849.110.
#' @format a \code{IsomirDataSeq} class.
#' @author Lorena Pantano, 2015-05-19
NULL