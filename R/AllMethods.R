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
#' @param object A `IsomirDataSeq` object.
#' @param value An integer matrix.
#' @param norm  Boolean, return log2-normalized counts.
#' @return [base::matrix] with raw or normalized count data.
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
#' @param object A [IsomirDataSeq] object.
#' @param mirna String referring to the miRNA to show.
#' @param minc Minimum number of isomiR reads needed
#'   to be included in the table.
#' @return [S4Vectors::DataFrame-class] with count
#' information. The `row.names`
#' show the isomiR names, and each of the columns shows the counts
#' for this isomiR in that sample. Mainly, it will return the count
#' matrix only for isomiRs belonging to the miRNA family given by
#' the `mirna` parameter. IsomiRs need to have counts bigger than
#' `minc` parameter at least in one sample to be included in the output.
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
    rawData <- metadata(object)[["rawData"]] %>% 
        filter(mir == mirna)

    des <- colData(object)
    is_subs = rawData[["mism"]] != "0"
    is_add = rawData[["add"]] != "0"
    is_t5 = rawData[["t5"]] != "0"
    is_t3 = rawData[["t3"]] != "0"
    is_ref = rawData[["mism"]] != "0" & rawData[["add"]] != "0" & rawData[["t5"]] != "0" & rawData[["t3"]] != "0"
    dt <- rawData %>% 
        mutate(uid = mir) %>% 
        mutate(uid = ifelse(is_ref,
                            paste0(uid, paste0(";ref")),
                            uid)) %>% 
        mutate(uid = ifelse(is_subs,
                            paste0(uid, paste0(";iso_snv:", mism)),
                            uid)) %>% 
        mutate(uid = ifelse(is_add,
                            paste0(uid, paste0(";iso_add3p:", add)),
                            uid)) %>% 
        mutate(uid = ifelse(is_t5,
                            paste0(uid, paste0(";iso_5p:", t5)),
                            uid)) %>% 
        mutate(uid = ifelse(is_t3,
                            paste0(uid, paste0(";iso_3p:", t3)),
                            uid)) %>% 
        .[,c("uid", rownames(des))] %>% 
        group_by(!!sym("uid")) %>% 
        summarise_all(funs(sum)) %>% 
        as.data.frame() %>% 
        column_to_rownames("uid") %>% 
        as.matrix()
    
    DataFrame(dt[ rowSums(dt[,2:ncol(dt), drop=FALSE] > minc) > 0, , drop=FALSE])
}

design.IsomirDataSeq <- function(object) object@design

#' Accessors for the 'design' slot of a IsomirDataSeq object.
#'
#' The design holds the R `formula` which expresses how the
#' counts depend on the variables in `colData`.
#' See [IsomirDataSeq] for details.
#'
#' @docType methods
#' @name design
#' @rdname design
#' @aliases design design,IsomirDataSeq-method design<-,IsomirDataSeq,formula-method
#' 
#' @param object A [IsomirDataSeq] object.
#' @param value A formula to pass to DESeq2.
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
setMethod(f = "isoSelect",
          signature = signature(object="IsomirDataSeq"),
          definition = isoSelect.IsomirDataSeq)
