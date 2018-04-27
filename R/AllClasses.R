#' Class that contains all isomiRs annotation for all samples
#'
#' The [IsomirDataSeq] is a subclass of
#' [SummarizedExperiment::SummarizedExperiment]
#' used to store the raw data, intermediate calculations and results of an
#' miRNA/isomiR analysis. This class stores all raw isomiRs
#' data for each sample, processed information,
#' summary for each isomiR type,
#' raw counts, normalized counts, and table with
#' experimental information for each sample.
#'
#' [IsomirDataSeqFromFiles] creates this object using seqbuster
#' output files.
#'
#' Methods for this objects are [isomiRs::counts()] to get
#' count matrix and [isomiRs::isoSelect()]
#' for miRNA/isomiR selection. Functions
#' available for this object are [isomiRs::isoCounts()] for
#' count matrix creation,
#' [isomiRs::isoNorm()] for normalization,
#' [isomiRs::isoDE()] for
#' differential expression and [isomiRs::isoPLSDA()] for clustering.
#' [isomiRs::isoPlot()] helps with basic expression plot.
#'
#' `metadata` contains two lists: `rawList` is a list with same
#' length than number of samples and stores the input files
#' for each sample; `isoList` is a list with same length than
#' number of samples and stores information for each isomiR type summarizing
#' the different changes for the different isomiRs (trimming at 3',
#' trimming a 5', addition and substitution). For instance, you can get
#' the data stored in `isoList` for sample 1 and 5' changes
#' with this code `metadata(ids)[['isoList']][[1]]$t5sum`.
#'
#' The naming of isomiRs follows these rules:
#'
#' * miRNA name
#' * type:ref if the sequence is the same than the miRNA reference.
#' `iso` if the sequence has variations.
#' * `t5 tag`:indicates variations at 5 position.
#' The naming contains two words: `direction - nucleotides`,
#' where direction can be UPPER CASE NT
#' (changes upstream of the 5 reference position) or
#' LOWER CASE NT (changes downstream of the 5 reference position).
#' `0` indicates no variation, meaning the 5 position is
#' the same than the reference. After `direction`,
#' it follows the nucleotide/s that are added (for upstream changes)
#'  or deleted (for downstream changes).
#' * `t3 tag`:indicates variations at 3 position.
#' The naming contains two words: `direction - nucleotides`,
#' where direction can be LOWER CASE NT
#' (upstream of the 3 reference position) or
#' UPPER CASE NT (downstream of the 3 reference position).
#' `0` indicates no variation, meaning the 3 position is
#' the same than the reference. After `direction`,
#' it follows the nucleotide/s that are added (for downstream changes)
#' or deleted (for upstream chanes).
#' * `ad tag`:indicates nucleotides additions at 3 position.
#' The naming contains two words: `direction - nucleotides`,
#' where direction is UPPER CASE NT
#' (upstream of the 5 reference position).
#' `0` indicates no variation, meaning the 3 position
#' has no additions. After `direction`,
#' it follows the nucleotide/s that are added.
#' * `mm tag`: indicates nucleotides substitutions along
#' the sequences. The naming contains three words:
#' `position-nucleotideATsequence-nucleotideATreference`.
#' * `seed tag`: same than `mm` tag,
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
#' de <- data.frame(row.names=c("f1" , "f2"),
#'                  condition = c("newborn", "newborn"))
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

setValidity("IsomirDataSeq", function(object) {
    if (!("counts" %in% names(assays(object))))
        return("the assays slot must contain a matrix named 'counts'")
    if (!is.numeric(counts(object)))
        return("the count data is not numeric")
    if (any(is.na(counts(object))))
        return("NA values are not allowed in the count matrix" )
    if (any( counts(object) < 0L))
        return("the count data contains negative values")
    TRUE
})

# Constructor
.IsomirDataSeq <- function(se, rawData=NULL, design=~1L){
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
    metadata(ids) <- list(rawData = rawData)
    ids
}


#' `IsomirDataSeqFromFiles`` loads miRNA annotation from seqbuster tool
#'
#' This function parses
#' output of seqbuster tool to allow isomiRs/miRNAs analysis of samples
#' in different groups such as
#' characterization, differential expression and clustering. It creates an
#' [isomiRs::IsomirDataSeq] object.
#'
#' @param files files with the output of seqbuster tool
#' @param coldata data frame containing groups for each sample
#' @param design a `formula` to pass to [DESeq2::DESeqDataSet]
#' @param rate minimum counts fraction to consider a mismatch a real mutation
#' @param canonicalAdd `boolean` only keep A/T non-template addition.
#'   All non-template nucleotides at the 3' end will be removed if they
#'   contain C/G nts.
#' @param uniqueMism `boolean` only keep mutations that have
#'   a unique hit to one miRNA molecule. For instance, if the sequence map
#'    to two different miRNAs, then it would be removed.
#' @param uniqueHits `boolean` whether filtering ambigous sequences or not.
#' @param minHits Minimum number of reads in the sample to consider it
#'   in the final matrix.
#' @param header boolean to indicate files contain headers
#' @param skip skip first line when reading files
#' @param quiet boolean indicating to print messages
#'   while reading files. Default `FALSE`.
#' @param ... arguments provided to
#'  [SummarizedExperiment::SummarizedExperiment]
#'   including rowData.
#' @details
#' This function parses the output of
#' http://seqcluster.readthedocs.org/mirna_annotation.html
#' for each sample to create a count matrix for isomiRs, miRNAs or
#' isomiRs grouped in
#' types (i.e all sequences with variations at 5' but ignoring any other type).
#' It creates
#' [isomiRs::IsomirDataSeq] object (see link to example usage of
#' this class)
#' to allow visualization, queries, differential
#' expression analysis and clustering.
#' To create the [isomiRs::IsomirDataSeq], it parses the isomiRs
#' files, and generates
#' an initial matrix having all isomiRs detected among samples. As well,
#' it creates
#' a summary for each isomiR type (trimming, addition and substitution) to
#' visualize general isomiRs distribution.
#'
#' @rdname IsomirDataSeqFromFiles
#' @name IsomirDataSeqFromFiles
#' @return
#' [IsomirDataSeq] class object.
#' @examples
#' path <- system.file("extra", package="isomiRs")
#' fn_list <- list.files(path, full.names = TRUE)
#' de <- data.frame(row.names=c("f1" , "f2"),
#'                  condition = c("newborn", "newborn"))
#' ids <- IsomirDataSeqFromFiles(fn_list, coldata=de)
#'
#' head(counts(ids))
#'
#' @export
IsomirDataSeqFromFiles <- function(files, coldata, rate = 0.2,
                                   canonicalAdd = TRUE, uniqueMism = TRUE,
                                   uniqueHits = FALSE,
                                   design = ~1L,
                                   minHits = 1L,
                                   header = TRUE, skip = 0, quiet = TRUE, ...){
    n_filtered = 0
    idx <- 0
    if (header == FALSE)
      skip = 1
    rawData <- lapply(files, function(f) {
        s <- rownames(coldata)[files==f]
        d <- as.data.frame(suppressMessages(read_tsv(f, skip = skip)),
                           stringsAsFactors = FALSE)
        if (quiet == FALSE)
          message("reading file: ", f)
        if (nrow(d) < 2) {
            n_filtered = n_filtered + 1
            message(paste0("This sample hasn't any lines: ", f))
            return(NULL)
        }else{
            d <- .filter_table(d, rate = rate, canonicalAdd = canonicalAdd,
                               uniqueMism = uniqueMism, uniqueHits = uniqueHits)
            if (nrow(d) < minHits){
                n_filtered = n_filtered + 1
                message("Skipping sample ", f,
                        ". Low number of hits according to minHits.")
                return(NULL)
            }
        }

        d %>% 
            unite("uid", seq, mir, mism, add, t5, t3, sep = ":") %>% 
            select(uid, freq) %>%
            gather(uid, freq) %>% 
            mutate(sample = s)
    }) %>% bind_rows() %>% 
        group_by(uid, sample) %>% 
        summarise(freq = sum(freq)) %>% 
        spread(sample, freq, fill = 0) %>% 
        separate(uid,
                 into = c("seq", "mir", "mism", "add", "t5", "t3"),
                 sep = ":")
    
    if (nrow(rawData) == 0)
        stop("No samples had valids miRNA hits.")

    countData <- IsoCountsFromMatrix(rawData, coldata)
    se <- SummarizedExperiment(assays = SimpleList(counts = countData),
                               colData = DataFrame(coldata), ...)
    ids <- .IsomirDataSeq(se, rawData, design)
    message("Total samples filtered due to low number of hits: ", n_filtered)
    return(ids)
}
