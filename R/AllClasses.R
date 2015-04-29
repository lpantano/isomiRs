#' IsomirDataSeq object and constructors
#' @rdname IsomirDataSeq
#' @export
IsomirDataSeq<-setClass("IsomirDataSeq",
                            contains = "SummarizedExperiment",
                            slots = list(
#                                counts="matrix",
#                                normcounts="matrix",
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
#' @name IsomirDataSeq
#' @rdname IsomirDataSeq
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
    se <- SummarizedExperiment(assays = SimpleList(counts=countData), colData = DataFrame(design), ...)
    IsoObj <- IsomirDataSeq(se, listSamples, listIsomirs, list())
    return(IsoObj)
}


setGeneric("isoraw", function(x, ...) standardGeneric("isoraw"))

setGeneric("isoraw<-",
           function(x, ...) standardGeneric("isoraw<-"))

setMethod(
    f = isoraw,
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "rawList")
    }
)

setReplaceMethod("isoraw", "IsomirDataSeq",
function(x, value)
    {
        slot(x, "rawList") <- value
        x

    }
)

setGeneric("isoinfo", function(x, ...) standardGeneric("isoinfo"))

setGeneric("isoinfo<-",
           function(x, ...) standardGeneric("isoinfo<-"))

setMethod(
    f = "isoinfo",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "isoList")
    }
)

setReplaceMethod("isoinfo", "IsomirDataSeq",
function(x, value)
    {
        slot(x, "isoList") <- value
        x
                     
    }
)

setGeneric("isostats", function(x, ...) standardGeneric("isostats"))

setGeneric("isostats<-",
           function(x, ...) standardGeneric("isostats<-"))

setMethod(
    f = "isostats",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "statsList")
    }
)

setReplaceMethod("isostats", "IsomirDataSeq",
                 function(x, value)
                 {
                     slot(x, "statsList") <- value
                     x
                     
                 }
)

setMethod(
    f = "print",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        head(colData(x))
    }
)


#' @rdname counts                                                                                                                  
#' @export
setMethod(
    f = "counts",
    signature = signature(object="IsomirDataSeq"),
    definition = function(object){
        assays(object)[['counts']]
    }
)

#' @name counts                                                                                                                    
#' @rdname counts                                                                                                                  
#' @exportMethod "counts<-"                                                                                                        
setReplaceMethod("counts", signature(object="IsomirDataSeq", value="matrix"),
    function(object, value){
    assays(object)[["counts"]] <- value
    object
})                     
