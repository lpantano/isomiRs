#' IsomirDataSeq object and constructors
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq<-setClass("IsomirDataSeq",
                            contains = "SummarizedExperiment",
                            slots = list(
#                                counts="matrix",
#                                normcounts="matrix",
                                varList="list",
                                expList="list",
                                sumList="list"
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
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq

IsomirDataSeq <- function(se, expList, varList, sumList){
    new("IsomirDataSeq", se, expList=expList, varList=varList, sumList=sumList)
}


#' The \code{IsomirDataSeqFromFiles}
#'
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq
#' @param files all samples
#' @param cov remove sequences that have relative abundance lower
#' than this number
#' @param design data frame containing groups for each sample
#' @param header files contain headers
#' @param skip skip first line when reading files
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
        out <- list(summary=0, t5sum=isomir.position(d, 6),
                    t3sum=isomir.position(d, 7),
                    subsum=subs.position(d, 4),
                    addsum=isomir.position(d, 5))
        listSamples[[row.names(design)[idx]]] <- d
        listIsomirs[[row.names(design)[idx]]] <- out
    }
    countData <- IsoCountsFromMatrix(listSamples, design)
    se <- SummarizedExperiment(assays = SimpleList(counts=countData), colData = design, ...)
    IsoObj <- IsomirDataSeq(se, listSamples, listIsomirs, list())
    # IsoObj@design <- design
    # IsoObj@expList <- listObj
    # IsoObj@varList <- listObjVar
    # IsoObj <- do.mir.table(IsoObj)
    return(IsoObj)
}


setGeneric("rawIso", function(x, ...) standardGeneric("rawIso"))

setGeneric("rawIso<-",
           function(x, ...) standardGeneric("rawIso<-"))

setMethod(
    f = rawIso,
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "expList")
    }
)

setReplaceMethod("rawIso", "IsomirDataSeq",
function(x, value)
    {
        slot(x, "expList") <- value
        x

    }
)
# setMethod(
#     f = "processIso",
#     signature = signature(x="IsomirDataSeq"),
#     definition = function(x){
#         x@varList
#     }
# )


# setMethod(
#     f = "summary",
#     signature = signature(x="IsomirDataSeq"),
#     definition = function(x){
#         x@sumList
#     }
# )


