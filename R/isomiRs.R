#' Do differential expression analysis with DESeq2
#'
#' This function does differential expression analysis with \code{\link[DESeq2]{DESeq2-package}}.
#' It will return a \code{\link{DESeqResults}} object.
#' 
#' @details 
#' 
#' This function collapses all isomiRs in different types.
#' Read more at \code{\link{isoCounts}}.
#' 
#' After that isomiRs, \code{\link[DESeq2]{DESeq2-package}} is used to do differential
#' expression analysis. It uses the design matrix given when the 
#' \code{\link{IsomirDataSeq}} is created
#' to construct a \code{\link[DESeq2]{DESeqDataSet}} object.
#'
#' @param ids object isomiDataSeq
#' @param formula used for DE analysis
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param subs differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @return \code{\link[DESeq2]{DESeqDataSet}} object
#' @examples
#' library(DESeq2)
#' data(isomiRexp)
#' dds <- isoDE(isomiRexp, formula=~condition)
#' @export
#' @import DESeq2
isoDE <- function(ids, formula, ref=FALSE, iso5=FALSE, iso3=FALSE,
                add=FALSE, subs=FALSE, seed=FALSE){
    if (ref | iso5 | iso3 | add | subs | seed){
        ids <- isoCounts(ids, ref, iso5, iso3, add, subs, seed)
    }
    countData <- counts(ids)
    dds <- DESeqDataSetFromMatrix(countData = countData,
                                colData = colData(ids),
                                design = formula)
    dds <- DESeq(dds, quiet=TRUE)
    dds
}

#' Heatmap of the top expressed isomiRs
#'
#' This function creates a heatmap with the top N
#' isomiRs/miRNAs. It uses the matrix under \code{counts(ids)}
#' and represent in a heatmap the raw counts for each sample.
#'
#' @param ids object isomiDataSeq
#' @param top number of isomiRs/miRNAs used
#' @examples
#' data(isomiRexp)
#' isoTop(isomiRexp)
#' @export
#' @import gplots
#' @import RColorBrewer
isoTop <- function(ids, top=20){
    select <- order(rowMeans(counts(ids)),
                    decreasing=TRUE)[1:top]
    hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
    heatmap.2(counts(ids)[select,], col = hmcol,
            scale="none", 
            dendrogram="none", trace="none")
}

#' Plot the amount of isomiRs in different samples
#'
#' This function plot different isomiRs proportion for each sample.
#' It can show trimming events at both side, additions and nucleotides
#' changes.
#'
#' @param ids object isomirDataSeq
#' @param type string (iso5, iso3, add, subs) to indicate what isomiRs
#' to use for the plot. See details for explanation.
#' @return \code{\link{IsomirDataSeq}} with new stored data to avoid
#' same calculation in the future.
#' @details 
#' There are four different values for \code{type} parameter. To plot 
#' trimming at 5' or 3' end, use \code{type="iso5"} or \code{type="iso3"}. 
#' In this case, it will plot 3 positions at both side of the reference
#' position described at miRBase. Each position refers to the number of 
#' sequences that start/end before or after the miRBase reference. The
#' color indicates the sample group. The size of the point is proportional
#' to the number of total counts. The position in \code{y} is the number of
#' different sequences. 
#' 
#' Same logic applies to \code{type="add"} and \code{type="subs"}. However,
#' when \code{type="add"}, the plot will refer to addition events from the
#' 3' end of the reference miRBase. Note that this additions doesn't match 
#' the precursor sequence. In this case, only 3 position after the 3' end
#' will appear in the plot. When \code{type="subs"}, it will appear one
#' position for each nucleotide in the reference miRNA. And the points
#' will indicate isomiRs with nucleotide changes at the given position.
#' 
#' @export
#' @import ggplot2
#' @examples
#' data(isomiRexp)
#' isoPlot(isomiRexp)
isoPlot <- function(ids, type="iso5"){
    freq <- size <- group <- abundance <- NULL
    codevn <- 2:5
    names(codevn) <- c("iso5", "iso3", "subs", "add")
    ratiov <- c(1 / 6, 1 / 6, 1 / 23, 1 / 3)
    names(ratiov) <- names(codevn)
    coden <- codevn[type]
    ratio <- ratiov[type]
    des <- colData(ids)
    table <- data.frame()
    isoList <- isoinfo(ids)
    for (sample in row.names(des)){
        uniq.dat <- as.data.frame( table(isoList[[sample]][[coden]]$size) )
        temp <- as.data.frame( isoList[[sample]][[coden]] %>%
                                group_by(size) %>%
                                summarise( freq=sum(freq) )
                              )
        total <- sum(temp$freq)
        temp <- merge(temp, uniq.dat, by=1)
        Total <- sum(temp$Freq)
        temp$abundance <- temp$freq / total
        temp$unique <- temp$Freq / Total
        table <- rbind( table,
                        data.frame( size=temp$size, abundance=temp$abundance,
                                        unique=temp$unique,
                                        sample=rep(sample, nrow(temp)),
                                        group=rep(des[sample,"condition"],
                                                nrow(temp)) ) )
    }
    isostats(ids)[[type]] <- table
    p <- ggplot(table) +
        geom_jitter(aes(x=factor(size),y=unique,colour=factor(group),
                        size=abundance)) +
        scale_colour_brewer("Groups",palette="Set1") +
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3")) +
        labs(list(title=paste(type,"distribution"),y="# of unique sequences",
                x="position respect to the reference"))
    print(p)
    ids
}

#' Create count matrix
#' 
#' This function collapses isomiRs into different groups.
#' 
#' @param ids IsomirDataSeq class
#' @param ref differenciate reference miRNA from rest
#' @param iso5 differenciate trimming at 5 miRNA from rest
#' @param iso3 differenciate trimming at 3 miRNA from rest
#' @param add differenciate additions miRNA from rest
#' @param subs differenciate nt substitution miRNA from rest
#' @param seed differenciate changes in 2-7 nt from rest
#' @param minc int minimum number of isomiR sequences
#' @details 
#' 
#' You can merge all isomiRs into miRNAs by calling the function only
#' with the first parameter \code{ids}. You can get a table with isomiRs and
#' the reference miRBase sequence by calling the function with \code{ref=TRUE}.
#' You can get a table with 5' trimming isomiRS, miRBase reference and
#' the rest by calling with \code{idx, ref=TRUE,iso5=TRUE}.
#' If you set up all parameters to TRUE, you will get a table for
#' each differnt sequence mapping to a miRNA (all isomiRs).
#' 
#' @return count table
#' @examples
#' data(isomiRexp)
#' ids <- isoCounts(isomiRexp, ref=TRUE)
#' head(counts(ids))
#' @export
isoCounts <- function(ids, ref=FALSE, iso5=FALSE, iso3=FALSE,
                      add=FALSE, subs=FALSE, seed=FALSE, minc=10){
        counts <- IsoCountsFromMatrix(isoraw(ids), colData(ids), ref,
                                      iso5, iso3,
                                      add, subs, seed, minc)
        se <- SummarizedExperiment(assays = SimpleList(counts=counts),
                                   colData = colData(ids))
        ids <- IsomirDataSeq(se, isoraw(ids), isoinfo(ids), isostats(ids))
        ids
}


#' Normalize count matrix
#'
#' This function normalizes raw count matrix using
#' \code{\link[DESeq2]{rlog}} function from \code{\link[DESeq2]{DESeq2-package}}.
#'
#' @param ids IsomirDataSeq object
#' @param formula formula that will be used for normalization
#' 
#' @return \code{\link{IsomirDataSeq}} object with the normalized
#' count matrix in a slot. The normalized matrix
#' can be access with \code{counts(ids, norm=TRUE)}.
#' 
#' @examples
#' library(DESeq2)
#' data(isomiRexp)
#' ids <- isoNorm(isomiRexp, formula=~condition)
#' head(counts(ids, norm=TRUE))
#' @export
#' @import DESeq2
isoNorm <- function(ids, formula=~condition){
    dds <- DESeqDataSetFromMatrix(countData = counts(ids),
                                colData = colData(ids),
                                design = formula)
    rld <- rlog(dds, blind=FALSE)
    normcounts(ids) <- assay(rld)
    ids
}