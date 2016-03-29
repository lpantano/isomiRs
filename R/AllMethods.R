
#' Accessors for the count matrix of a IsomirDataSeq object.
#'
#' The counts slot holds the count data as a matrix of non-negative integer
#' count values, one row for each isomiR, and one
#' column for each sample. The normalized matrix
#' can be obtained by using the parameter \code{norm=TRUE}.
#'
#' @docType methods
#' @name counts
#' @rdname counts
#' @aliases counts counts,IsomirDataSeq-method counts<-,IsomirDataSeq,matrix-method
#'
#' @param object a \code{IsomirDataSeq} object
#' @param value an integer matrix
#' @param norm TRUE return log2-normalized counts
#' @return \code{\link[base]{matrix}} with raw or normalized count data.
#' @author Lorena Pantano
#' @examples
#' data(mirData)
#' head(counts(mirData))
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


#' Method to select specific miRNAs from an IsomirDataSeq object.
#'
#' This method allows to select a miRNA and all its isomiRs
#' from the count matrix.
#'
#' @docType methods
#' @name isoSelect
#' @rdname isoSelect
#' @aliases isoSelect isoSelect,IsomirDataSeq-method
#'
#' @param object a \code{IsomirDataSeq} object.
#' @param mirna string referring to the miRNA to show
#' @param minc int minimum number of isomiR reads needed
#' to be included in the table.
#' @return \code{\link[S4Vectors]{DataFrame-class}} with count
#' information. The \code{row.names}
#' show the isomiR names, and each of the columns shows the counts
#' for this isomiR in that sample. Mainly, it will return the count
#' matrix only for isomiRs belonging to the miRNA family given by
#' the \code{mirna} parameter. IsomiRs need to have counts bigger than
#' \code{minc} parameter to be included in the output.
#' 
#' @author Lorena Pantano
#' 
#' @examples
#' data(mirData)
#' # To select isomiRs from let-7a-5p miRNA
#' # and with 10000 reads or more.
#' isoSelect(mirData, mirna="hsa-let-7a-5p", minc=10000)
#' @export
isoSelect.IsomirDataSeq <- function(object, mirna,  minc=10) {
    x <- metadata(object)$rawList
    if ( mirna == "" )
        stop("mirna parameter needs to have a value")
    l <- lapply( x, function(sample){
        sample %>% filter( mir == mirna )
    })
    df <- as.matrix(IsoCountsFromMatrix(l, colData(object), ref=TRUE,iso5=TRUE,iso3=TRUE,
              add=TRUE, subs=TRUE, seed=TRUE))
    df[ df < minc ] <- 0
    DataFrame(df[ rowSums(df) > 0, , drop=FALSE])
}


#' @rdname isoSelect
#' @export
setMethod(f="isoSelect",
          signature = signature(object="IsomirDataSeq"),
          definition = isoSelect.IsomirDataSeq)
