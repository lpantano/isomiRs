#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq <- setClass("IsomirDataSeq",
                            contains = "RangedSummarizedExperiment",
                            representation = representation(
                                isoList="list",
                                rawList="list",
                                statsList="list"
                        ))

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

#' Class that contain all isomiRs annotation for all samples
#' 
#' The \code{\link{IsomirDataSeq}} is a subclass of 
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' used to store the raw data, intermediate calculations and results of an
#' isomiR analysis.  The \code{\link{IsomirDataSeq}} class stores all raw isomiRs
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
#' @param se \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param expList list of samples with seqbuster output.
#' @param isoList list of samples with summarized isomiR 
#' information for each type.
#' @param statsList list of samples with general isomiR information.
#' Could be empty list.
#' 
#' @aliases IsomirDataSeq IsomirDataSeq-class
#' @examples 
#' \dontrun{
#' fn_list = c("url1", "url2")
#' de = data.frame(row.names=c("f1" , "f2"), condition = c("n1", "o1")) 
#' ids <- IsomirDataSeqFromFiles(fn_list, design=de)
#' 
#' select(ids, "hsa-let-7a-5p")
#' counts(ids)[1:5, ]
#' }
#' @docType class
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq <- function(se, expList, isoList, statsList){
    if (!is(se, "RangedSummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
                  se <- as(se, "RangedSummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
                  # only to help transition from SummarizedExperiment to new
                  # RangedSummarizedExperiment objects,
                  # remove once transition is complete
                  se <- as(se, "RangedSummarizedExperiment")
        } else {
                  stop("'se' must be a RangedSummarizedExperiment object")
        }
    }
    new("IsomirDataSeq", se, rawList=expList,
        isoList=isoList, statsList=statsList)
}


#' \code{IsomirDataSeqFromFiles} loads miRNA annotation from seqbuster tool
#'
#' This function parses \url{http://seqcluster.readthedocs.org/mirna_annotation.html}
#' output to allow isomiRs/miRNAs analysis of samples in different groups such as
#' characterization, differential expression and clustering. It creates
#' \code{\link[isomiRs]{IsomirDataSeq}} object.
#'
#' @param files files with the output of seqbuster tool
#' @param cov remove sequences that have relative abundance lower
#' than this number
#' @param design data frame containing groups for each sample
#' @param header boolean to indicate files contain headers
#' @param skip skip first line when reading files
#' @param ... arguments provided to \code{\link[SummarizedExperiment]{SummarizedExperiment}} including rowData and exptData
#' @details
#' This function parse the output of \url{http://seqcluster.readthedocs.org/mirna_annotation.html}
#' for each sample to create count matrix for isomiRs, miRNAs or isomiRs grouped in
#' types (i.e all sequences with variations at 5' but ignoring any other type). It creates
#' \code{\link[isomiRs]{IsomirDataSeq}} object to allow visualization, queries, differential
#' expression analysis and clustering.
#' To create the \code{\link[isomiRs]{IsomirDataSeq}}, it parses the isomiRs file, and generates
#' an initial matrix having all miRNAs detected among samples. As well, it creates
#' a summary for each isomiR type (trimming, addition and substitution.) to
#' visualize general isomiRs distribution omong samples.
#'
#' @rdname IsomirDataSeqFromFiles
#' @name IsomirDataSeqFromFiles
#' @return
#' \code{\link{IsomirDataSeq}} class
#' @export
IsomirDataSeqFromFiles <- function(files, design, cov=1,
                                   header=FALSE, skip=1, ...){
    listSamples <- vector("list")
    listIsomirs <- vector("list")
    idx <- 0
    for (f in files){
        idx <- idx + 1
        d <- read.table(f, header=header, skip=skip)
        if (ncol(d) < 2){
            warning(paste0("This sample has not lines: ", f))
        }else{
            d <- .filter_table(d, cov)
            out <- list(summary=0,
                        t5sum = .isomir_position(d, 6),
                        t3sum = .isomir_position(d, 7),
                        subsum = .subs_position(d, 4),
                        addsum = .isomir_position(d, 5)
            )
            listSamples[[row.names(design)[idx]]] <- d
            listIsomirs[[row.names(design)[idx]]] <- out
        }
    }
    countData <- IsoCountsFromMatrix(listSamples, design)
    se <- SummarizedExperiment(assays = SimpleList(counts=countData),
                               colData = DataFrame(design), ...)
    ids <- IsomirDataSeq(se, listSamples, listIsomirs, list())
    return(ids)
}
