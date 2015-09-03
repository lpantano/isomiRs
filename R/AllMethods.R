
# Accessors for the 'isoraw' slot of a IsomirDataSeq object.
setMethod(
    f = isoraw,
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "rawList")
    }
)

setReplaceMethod("isoraw", "IsomirDataSeq",
                 function(x, value){
                     slot(x, "rawList") <- value
                     x
                 }
)

setMethod(
    f = "isoinfo",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "isoList")
    }
)

setReplaceMethod("isoinfo", "IsomirDataSeq",
                 function(x, value){
                     slot(x, "isoList") <- value
                     validObject(x)
                     x
                 }
)

setMethod(
    f = "isostats",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        slot(x, "statsList")
    }
)

setReplaceMethod("isostats", "IsomirDataSeq",
                 function(x, value){
                     slot(x, "statsList") <- value
                     validObject(x)
                     x
                 }
)


#' Accessors for the 'counts' slot of a IsomirDataSeq object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each isomiR, and one
#' column for each sample. The normalized matrix
#' can be obtained by \code{normcounts} method.
#'
#' @usage
#' \S4method{counts}{IsomirDataSeq}(object, norm=FALSE)
#'
#' \S4method{counts}{IsomirDataSeq,matrix}(object)<-value
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,IsomirDataSeq-method counts<-,IsomirDataSeq,matrix-method
#'
#' @param object a \code{IsomirDataSeq} object.
#' @param value an integer matrix
#' @param norm TRUE return log2-normalized counts
#' @author Lorena Pantano
#' @examples
#' data(isomiRexp)
#' head(counts(isomiRexp))
#' @export
counts.IsomirDataSeq <- function(object, norm=FALSE) {
    if (norm){
        return(assays(object)[["norm"]])
    }
    assays(object)[["counts"]]
}

#' @rdname counts
#' @export
setMethod("counts", signature(object="IsomirDataSeq"), counts.IsomirDataSeq)

#' @name counts
#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="IsomirDataSeq", value="matrix"),
                 function(object, value){
                     assays(object)[["counts"]] <- value
                     validObject(object)
                     object
                 }
)


setMethod(
    f = "normcounts",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        assays(x)[["norm"]]
    }
)

setReplaceMethod("normcounts", "IsomirDataSeq",
                 function(x, value){
                     assays(x)[["norm"]] <- value
                     validObject(x)
                     x
                 }
)


#' Method to browse an IsomirDataSeq object.
#'
#' This method allows to select a miRNA and all its isomiRs
#' along all samples.
#'
#' @usage
#' \S4method{isoSelect}{IsomirDataSeq}(object, mirna="", norm=FALSE, minc=10)
#'
#'
#' @docType methods
#' @name isoSelect
#' @rdname isoSelect
#' @aliases isoSelect isoSelect,IsomirDataSeq-method
#'
#' @param object a \code{IsomirDataSeq} object.
#' @param mirna string of the miRNA to show
#' @param norm TRUE return log2-normalized counts
#' @param minc int minimum number of isomiR reads
#' @author Lorena Pantano
#' 
#' @examples
#' data(isomiRexp)
#' isoSelect(isomiRexp, mirna="hsa-let-7a-5p")
#' @export
isoSelect.IsomirDataSeq <- function(object, mirna="", norm=FALSE, minc=10) {
    x <- isoraw(object)
    if ( mirna == "" )
        stop("mirna parameter needs to have a value")
    l <- lapply( x, function(sample){
        sample %>% filter( mir == mirna )
    })
    IsoCountsFromMatrix(l, colData(object), ref=TRUE,iso5=TRUE,iso3=TRUE,
              add=TRUE, subs=TRUE, seed=TRUE, minc=minc)
}


#' @rdname isoSelect
#' @export
setMethod(f="isoSelect",
          signature = signature(object="IsomirDataSeq"),
          definition = isoSelect.IsomirDataSeq)
