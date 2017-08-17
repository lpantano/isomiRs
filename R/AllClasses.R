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
#' The naming of isomiRs follows these rules:
#'
#' * miRNA name
#' * type:ref if the sequence is the same than the miRNA reference. \code{iso} if the sequence has variations.
#' * \code{t5 tag}:indicates variations at 5 position.
#' The naming contains two words: \code{direction - nucleotides},
#' where direction can be UPPER CASE NT
#' (changes upstream of the 5 reference position) or
#' LOWER CASE NT (changes downstream of the 5 reference position).
#' \code{0} indicates no variation, meaning the 5 position is
#' the same than the reference. After \code{direction},
#' it follows the nucleotide/s that are added (for upstream changes)
#'  or deleted (for downstream changes).
#' * \code{t3 tag}:indicates variations at 3 position.
#' The naming contains two words: \code{direction - nucleotides},
#' where direction can be LOWER CASE NT
#' (upstream of the 3 reference position) or
#' UPPER CASE NT (downstream of the 3 reference position).
#' \code{0} indicates no variation, meaning the 3 position is
#' the same than the reference. After \code{direction},
#' it follows the nucleotide/s that are added (for downstream changes)
#' or deleted (for upstream chanes).
#' * \code{ad tag}:indicates nucleotides additions at 3 position.
#' The naming contains two words: \code{direction - nucleotides},
#' where direction is UPPER CASE NT
#' (upstream of the 5 reference position).
#' \code{0} indicates no variation, meaning the 3 position
#' has no additions. After \code{direction},
#' it follows the nucleotide/s that are added.
#' * \code{mm tag}: indicates nucleotides substitutions along
#' the sequences. The naming contains three words:
#' \code{position-nucleotideATsequence-nucleotideATreference}.
#' * \code{seed tag}: same than \code{mm} tag,
#' but only if the change happens between nucleotide 2 and 8.
#'
#' In general nucleotides in UPPER case mean insertions respect
#' to the reference sequence, and nucleotides in LOWER case
#' mean deletions respect to the reference sequence.
#'
#' @aliases IsomirDataSeq-class
#' @examples
#' path <- system.file("extra", package="isomiRs")
#' fn_list <- list.files(path, full.names = TRUE)
#' de <- data.frame(row.names=c("f1" , "f2"), condition = c("newborn", "newborn"))
#' ids <- IsomirDataSeqFromFiles(fn_list, coldata=de)
#'
#' head(counts(ids))
#'
#' @rdname IsomirDataSeq
#' @md
#' @exportClass "IsomirDataSeq"
IsomirDataSeq <- setClass("IsomirDataSeq",
                          contains = "SummarizedExperiment",
                          representation = representation(
                              design = "formula"))

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
.IsomirDataSeq <- function(se, rawList=NULL, isoList=NULL, design=~1){
    if (!is(se, "SummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
                  se <- as(se, "SummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
                  # only to help transition from SummarizedExperiment to new
                  # RangedSummarizedExperiment objects,
                  # remove once transition is complete
                  se <- as(se, "SummarizedExperiment")
        } else {
                  stop("'se' must be a SummarizedExperiment object")
        }
    }
    ids <- new("IsomirDataSeq", se, design = design)
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
#' @param coldata data frame containing groups for each sample
#' @param design a \code{formula} to pass to \code{\link[DESeq2]{DESeqDataSet}}
#' @param rate minimum counts fraction to consider a mismatch a real mutation
#' @param canonicalAdd \code{boolean} only keep A/T non-template addition.
#' All non-template nucleotides at the 3' end will be removed if they
#' contain C/G nts.
#' @param uniqueMism \code{boolean} only keep mutations that have
#' a unique hit to one miRNA molecule
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
#' ids <- IsomirDataSeqFromFiles(fn_list, coldata=de)
#'
#' head(counts(ids))
#'
#' @export
IsomirDataSeqFromFiles <- function(files, coldata, rate=0.2,
                                   canonicalAdd=TRUE, uniqueMism=TRUE,
                                   design = ~1,
                                   header=TRUE, skip=0, quiet=TRUE, ...){
    listSamples <- vector("list")
    listIsomirs <- vector("list")
    idx <- 0
    if (header == FALSE)
      skip = 1
    for (f in files){
        idx <- idx + 1
        d <- as.data.frame(suppressMessages(read_tsv(f, skip=skip)),
                           stringsAsFactors=FALSE)
        if (quiet == FALSE)
          cat("reading file: ", f, "\n")
        if (nrow(d) < 2){
            warning(paste0("This sample hasn't any lines: ", f))
        }else{
            d <- .filter_table(d, rate=rate, canonicalAdd=canonicalAdd,
                               uniqueMism=uniqueMism)
            out <- list(summary = 0,
                        t5sum = .isomir_position(d, 6),
                        t3sum = .isomir_position(d, 7),
                        subsum = .subs_position(d, 4),
                        addsum = .isomir_position(d, 5)
            )
            listSamples[[row.names(coldata)[idx]]] <- d
            listIsomirs[[row.names(coldata)[idx]]] <- out
        }
    }
    coldata = coldata[names(listSamples),,drop=FALSE]
    countData <- IsoCountsFromMatrix(listSamples, coldata)
    se <- SummarizedExperiment(assays = SimpleList(counts=countData),
                               colData = DataFrame(coldata), ...)
    ids <- .IsomirDataSeq(se, listSamples, listIsomirs, design)
    return(ids)
}
