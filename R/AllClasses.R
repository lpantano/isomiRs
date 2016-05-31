#' Class that contains all isomiRs annotation for all samples
#'
#' The \code{\link{IsomirDataSeq}} is a subclass of
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' used to store the raw data, intermediate calculations and results of an
#' miRNA/isomiR analysis. This class stores all raw isomiRs
#' data for each sample, processed information,
#' summary for each isomiR type,
#' raw counts, normalized counts, and table with
#' experimental information for each sample.
#'
#' \code{\link{IsomirDataSeqFromFiles}} creates this object using seqbuster output files.
#'
#' Methods for this objects are \code{\link[isomiRs]{counts}} to get count matrix
#' and \code{\link[isomiRs]{isoSelect}}
#' for miRNA/isomiR selection. Functions
#' available for this object are \code{\link[isomiRs]{isoCounts}} for count matrix creation,
#' \code{\link[isomiRs]{isoNorm}} for normalization, \code{\link[isomiRs]{isoDE}} for
#' differential expression and \code{\link{isoPLSDA}} for clustering.
#' \code{\link[isomiRs]{isoPlot}} helps with basic expression plot.
#'
#' \code{metadata} contains two lists: \code{rawList} is a list with same
#' length than number of samples and stores the input files
#' for each sample; \code{isoList} is a list with same length than
#' number of samples and stores information for each isomiR type summarizing
#' the different changes for the different isomiRs (trimming at 3',
#' trimming a 5', addition and substitution). For instance, you can get
#' the data stored in \code{isoList} for sample 1 and 5' changes
#' with this code \code{metadata(ids)[['isoList']][[1]]$t5sum}.
#'
#' @aliases IsomirDataSeq-class
#' @examples
#' path <- system.file("extra", package="isomiRs")
#' fn_list <- list.files(path, full.names = TRUE)
#' de <- data.frame(row.names=c("f1" , "f2"), condition = c("newborn", "newborn"))
#' ids <- IsomirDataSeqFromFiles(fn_list, design=de)
#'
#' head(counts(ids))
#'
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq <- setClass("IsomirDataSeq",
                          contains = "RangedSummarizedExperiment")

setValidity( "IsomirDataSeq", function( object ) {
    if (!("counts" %in% names(assays(object))))
        return( "the assays slot must contain a matrix named 'counts'" )
    if ( !is.numeric( counts(object) ) )
        return( "the count data is not numeric" )
    if ( any( is.na( counts(object) ) ) )
        return( "NA values are not allowed in the count matrix" )
    if ( any( counts(object) < 0 ) )
        return( "the count data contains negative values" )
    TRUE
} )

# Constructor
.IsomirDataSeq <- function(se, rawList=NULL, isoList=NULL){
    if (!is(se, "RangedSummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
                  se <- as(se, "RangedSummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
                  # only to help transition from SummarizedExperiment to new
                  # RangedSummarizedExperiment objects,
                  # remove once transition is complete
                  se <- as(se, "RangedSummarizedExperiment")
        } else {
                  stop("'se' must be a SummarizedExperiment object")
        }
    }
    ids <- new("IsomirDataSeq", se)
    metadata(ids) <- list(rawList = rawList, isoList = isoList)
    ids
}


#' \code{IsomirDataSeqFromFiles} loads miRNA annotation from seqbuster tool
#'
#' This function parses
#' output of seqbuster tool to allow isomiRs/miRNAs analysis of samples
#' in different groups such as
#' characterization, differential expression and clustering. It creates an
#' \code{\link[isomiRs]{IsomirDataSeq}} object.
#'
#' @param files files with the output of seqbuster tool
#' @param design data frame containing groups for each sample
#' @param header boolean to indicate files contain headers
#' @param skip skip first line when reading files
#' @param quiet boolean indicating to print messages
#'  while reading files. Default \code{FALSE}.
#' @param ... arguments provided to \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' including rowData.
#' @details
#' This function parses the output of \url{http://seqcluster.readthedocs.org/mirna_annotation.html}
#' for each sample to create a count matrix for isomiRs, miRNAs or isomiRs grouped in
#' types (i.e all sequences with variations at 5' but ignoring any other type). It creates
#' \code{\link[isomiRs]{IsomirDataSeq}} object (see link to example usage of this class)
#' to allow visualization, queries, differential
#' expression analysis and clustering.
#' To create the \code{\link[isomiRs]{IsomirDataSeq}}, it parses the isomiRs files, and generates
#' an initial matrix having all isomiRs detected among samples. As well, it creates
#' a summary for each isomiR type (trimming, addition and substitution) to
#' visualize general isomiRs distribution.
#'
#' @rdname IsomirDataSeqFromFiles
#' @name IsomirDataSeqFromFiles
#' @return
#' \code{\link{IsomirDataSeq}} class object.
#' @examples
#' path <- system.file("extra", package="isomiRs")
#' fn_list <- list.files(path, full.names = TRUE)
#' de <- data.frame(row.names=c("f1" , "f2"), condition = c("newborn", "newborn"))
#' ids <- IsomirDataSeqFromFiles(fn_list, design=de)
#'
#' head(counts(ids))
#'
#' @export
IsomirDataSeqFromFiles <- function(files, design,
                                   header=FALSE, skip=1, quiet=TRUE, ...){
    listSamples <- vector("list")
    listIsomirs <- vector("list")
    idx <- 0
    if (header == TRUE)
      skip = 0
    for (f in files){
        idx <- idx + 1
        d <- read.table(f, header=header, skip=skip, stringsAsFactors = FALSE)
        if (quiet == FALSE)
          cat("reading file: ", f, "\n")
        if (nrow(d) < 2){
            warning(paste0("This sample hasn't any lines: ", f))
        }else{
            d <- .filter_table(d)
            out <- list(summary = 0,
                        t5sum = .isomir_position(d, 6),
                        t3sum = .isomir_position(d, 7),
                        subsum = .subs_position(d, 4),
                        addsum = .isomir_position(d, 5)
            )
            listSamples[[row.names(design)[idx]]] <- d
            listIsomirs[[row.names(design)[idx]]] <- out
        }
    }
    design = design[names(listSamples),,drop=FALSE]
    countData <- IsoCountsFromMatrix(listSamples, design)
    se <- SummarizedExperiment(assays = SimpleList(counts=countData),
                               colData = DataFrame(design), ...)
    # ids <- new("IsomirDataSeq", se)
    ids <- .IsomirDataSeq(se, listSamples, listIsomirs)
    return(ids)
}
