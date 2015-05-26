#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq<-setClass("IsomirDataSeq",
                            contains = "SummarizedExperiment",
                            slots = list(
                                isoList="list",
                                rawList="list",
                                statsList="list"
                        ))

# setValidity( "DESeqDataSet", function( object ) {

#' The \code{IsomirDataSeq} is a subclass of \code{SummarizedExperiment},
#' used to store the input values, intermediate calculations and results of an
#' isomiR analysis.  The \code{IsomirDataSeq} class stores all raw isomiRs
#' data for each sample, processed information,
#' summary for each isomiR type,
#' raw counts, normalized counts, and data.frame with
#' columns information for each sample.
#'
#' @param se SummarizedExperiment object
#' @param expList list of samples with miraligner output
#' @param varList list of samples with summarized isomiR 
#' information for each type
#' @param sumList list of samples with general isomiR information 
#' 
#' @aliases IsomirDataSeq IsomirDataSeq-class IsomirDataSeqFromFiles
#' 
#' @docType class
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq <- function(se, expList, varList, sumList){
    new("IsomirDataSeq", se, rawList=expList, isoList=varList, statsList=sumList)
}


#' The \code{IsomirDataSeqFromFiles}
#'
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq
#' @param files files with the output of miraligner tool
#' @param cov remove sequences that have relative abundance lower
#' than this number
#' @param design data frame containing groups for each sample
#' @param header boolean to indicate files contain headers
#' @param skip skip first line when reading files
#' @param ... arguments provided to \code{SummarizedExperiment} including rowData and exptData
#' @return
#' \code{IsomirDataSeq} class
#' @export
IsomirDataSeqFromFiles <- function(files, design, cov=1, header=FALSE, skip=1, ...){
    listSamples <- vector("list")
    listIsomirs <- vector("list")
    idx <- 0
    for (f in files){
        idx <- idx + 1
        # print(idx)
        d <- read.table(f, header=header, skip=skip)

        d <- filter.table(d, cov)
        out <- list(summary=0, 
                    t5sum = .isomir_position(d, 6),
                    t3sum = .isomir_position(d, 7),
                    subsum = .subs_position(d, 4),
                    addsum = .isomir_position(d, 5)
                    )
        listSamples[[row.names(design)[idx]]] <- d
        listIsomirs[[row.names(design)[idx]]] <- out
    }
    countData <- IsoCountsFromMatrix(listSamples, design)
    se <- SummarizedExperiment(assays = SimpleList(counts=countData), colData = DataFrame(design), ...)
    IsoObj <- IsomirDataSeq(se, listSamples, listIsomirs, list())
    return(IsoObj)
}


