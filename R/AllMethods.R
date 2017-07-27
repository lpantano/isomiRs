
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
counts.IsomirDataSeq <- function(object, norm=FALSE) {
    if (norm){
        return(assays(object)[["norm"]])
    }
    assays(object)[["counts"]]
}

#' @rdname counts
#' @exportMethod "counts"
setMethod("counts", signature(object="IsomirDataSeq"), counts.IsomirDataSeq)

#' @rdname counts
#' @exportMethod "counts<-"
setReplaceMethod("counts", signature(object="IsomirDataSeq", value="matrix"),
                 function(object, value){
                     assays(object)[["counts"]] <- value
                     validObject(object)
                     object
                 }
)

# normcounts
setMethod(
    f = "normcounts",
    signature = signature(x="IsomirDataSeq"),
    definition = function(x){
        assays(x)[["norm"]]
    }
)

# normcounts
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
#' \code{minc} parameter at least in one sample to be included in the output.
#' Annotation of isomiRs follows these rules:
#' 
#' * miRNA name
#' * mismatches
#' * additions
#' * 5 trimming events
#' * 3 trimming events
#'
#' @author Lorena Pantano
#'
#' @examples
#' data(mirData)
#' # To select isomiRs from let-7a-5p miRNA
#' # and with 10000 reads or more.
#' isoSelect(mirData, mirna="hsa-let-7a-5p", minc=10000)
isoSelect.IsomirDataSeq <- function(object, mirna,  minc=10) {
    x <- metadata(object)$rawList
    if ( mirna == "" )
        stop("mirna parameter needs to have a value")
    l <- lapply( x, function(sample){
        sample %>% filter( mir == mirna )
    })
    .list_samples = lapply(row.names(colData(object)), function(sample){
        d <- l[[sample]] %>%
            mutate(id=paste(mir, mism, add, t5, t3, ":", seq)) %>%
            select(id, freq) %>% mutate(sample=sample)
        d
    })
    df <- bind_rows(.list_samples) %>%
        #group_by(id, sample) %>%
        #summarise(freq=sum(freq)) %>% ungroup() %>%
        spread(key=sample, value=freq, fill=0)
    DataFrame(df[ rowSums(df[,2:ncol(df)] > minc ) > 0, , drop=FALSE])
}

design.IsomirDataSeq <- function(object) object@design

#' Accessors for the 'design' slot of a IsomirDataSeq object.
#'
#' The design holds the R \code{formula} which expresses how the
#' counts depend on the variables in \code{colData}.
#' See \code{\link{IsomirDataSeq}} for details.
#'
#' @docType methods
#' @name design
#' @rdname design
#' @aliases design design,IsomirDataSeq-method design<-,IsomirDataSeq,formula-method
#' @param object a \code{IsomirDataSeq} object
#' @param value a \code{formula} to pass to DESeq2
#' @examples
#'
#' data(mirData)
#' design(mirData) <- formula(~ 1)
#' @return design for the experiment
#' @exportMethod "design"
setMethod("design", signature(object="IsomirDataSeq"), design.IsomirDataSeq)

#' @name design
#' @rdname design
#' @exportMethod "design<-"
setReplaceMethod("design", signature(object="IsomirDataSeq", value="formula"),
                 function( object, value ) {
                     # Temporary hack for backward compatibility with "old"
                     # IsomirDataSeq objects. Remove once all serialized
                     # IsomirDataSeq objects around have been updated.
                     if (!.hasSlot(object, "rowRanges"))
                         object <- updateObject(object)
                     object@design <- value
                     validObject(object)
                     object
                 })

#' @rdname isoSelect
#' @exportMethod "isoSelect"
setMethod(f="isoSelect",
          signature = signature(object="IsomirDataSeq"),
          definition = isoSelect.IsomirDataSeq)
