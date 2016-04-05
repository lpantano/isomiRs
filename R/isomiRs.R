#' Differential expression analysis with DESeq2
#'
#' This function does differential expression analysis with
#' \code{\link[DESeq2]{DESeq2-package}} using the specific formula.
#' It will return a \code{\link[DESeq2]{DESeqDataSet}} object.
#'
#' @details
#'
#' First, this function collapses all isomiRs in different types.
#' Read more at \code{\link{isoCounts}} to know the different options
#' available to collapse isomiRs.
#'
#' After that, \code{\link[DESeq2]{DESeq2-package}} is used to do differential
#' expression analysis. It uses the count matrix and design experiment
#' stored at (\code{counts(ids)} and \code{colData(ids)})
#' \code{\link{IsomirDataSeq}} object
#' to construct a \code{\link[DESeq2]{DESeqDataSet}} object.
#'
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param formula used for DE analysis
#' @param ... options to pass to \link{isoCounts} including
#' ref, iso5, iso3, add, subs and seed parameters.
#' 
#' @return \code{\link[DESeq2]{DESeqDataSet}} object. 
#' To get the differential expression isomiRs, use \link[DESeq2]{results} from
#' DESeq2 package. This allows to ask for different contrast
#' without calling again \code{isoDE}. Read \code{results} 
#' manual to know how to access all the information.
#' 
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, minc=10, mins=6)
#' dds <- isoDE(mirData, formula=~condition)
#' @export
isoDE <- function(ids, formula, ...){
    ids <- isoCounts(ids, ...)
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
#' to get the top expressed isomiRs/miRNAs using the average
#' expression value
#' and plot a heatmap with the raw counts for each sample.
#'
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param top number of isomiRs/miRNAs used
#' @examples
#' data(mirData)
#' isoTop(mirData)
#' @export
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
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param type string (iso5, iso3, add, subs) to indicate what isomiRs
#' to use for the plot. See details for explanation.
#' @param column string indicating the column in
#' \code{colData} to color samples.
#' @return \link[ggplot2]{ggplot} object showing different isomiRs changes at
#' different positions. 
#' @details
#' There are four different values for \code{type} parameter. To plot
#' trimming at 5' or 3' end, use \code{type="iso5"} or \code{type="iso3"}.
#' In this case, it will plot 3 positions at both side of the reference
#' position described at miRBase site. Each position refers to the number of
#' sequences that start/end before or after the miRBase reference. The
#' color indicates the sample group. The size of the point is proportional
#' to the number of total counts. The position at \code{y} is the number of
#' different sequences.
#'
#' Same logic applies to \code{type="add"} and \code{type="subs"}. However,
#' when \code{type="add"}, the plot will refer to addition events from the
#' 3' end of the reference position. Note that this additions don't match
#' to the precursor sequence, they are non-template additions.
#' In this case, only 3 positions after the 3' end
#' will appear in the plot. When \code{type="subs"}, it will appear one
#' position for each nucleotide in the reference miRNA. Points
#' will indicate isomiRs with nucleotide changes at the given position.
#'
#' @export
#' @examples
#' data(mirData)
#' isoPlot(mirData)
isoPlot <- function(ids, type="iso5", column="condition"){
    freq <- size <- group <- abundance <- NULL
    codevn <- 2:5
    names(codevn) <- c("iso5", "iso3", "subs", "add")
    ratiov <- c(1 / 6, 1 / 6, 1 / 23, 1 / 3)
    names(ratiov) <- names(codevn)
    coden <- codevn[type]
    ratio <- ratiov[type]
    des <- colData(ids)
    table <- data.frame()
    isoList <- metadata(ids)$isoList
    for (sample in row.names(des)){
        if (nrow(isoList[[sample]][[coden]]) > 0 ){
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
                                      group=rep(des[sample, column],
                                                nrow(temp)) ) )
        }
    }
    p <- ggplot(table) +
        geom_jitter(aes(x=factor(size),y=unique,colour=factor(group),
                        size=abundance)) +
        scale_colour_brewer("Groups",palette="Set1") +
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3")) +
        labs(list(title=paste(type,"distribution"),y="# of unique sequences",
                x="position respect to the reference"))
    print(p)
}
#' Plot nucleotides changes at a given position
#'
#' This function plot different isomiRs proportion for each sample at a given
#' position focused on the nucleotide change that happens there.
#'
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param position integer indicating the position to show
#' @param column string indicating the column in
#' colData to color samples.
#' @return \link[ggplot2]{ggplot} object showing nucleotide changes
#' at a given position. 
#' @details
#' It shows the nucleotides changes at the given position for each
#' sample in each group.
#' The color indicates the sample group. The size of the point is proportional
#' to the number of total counts of isomiRs with changes.
#' The position at \code{y} is the number of different sequences
#' supporting the change.
#'
#' @export
#' @examples
#' data(mirData)
#' isoPlotPosition(mirData)
isoPlotPosition <- function(ids, position=1, column="condition"){
    freq <- size <- group <- abundance <- NULL
    codevn <- 2:5
    type <- "subs"
    names(codevn) <- c("iso5", "iso3", "subs", "add")
    ratiov <- c(1 / 6, 1 / 6, 1 / 23, 1 / 3)
    names(ratiov) <- names(codevn)
    coden <- codevn[type]
    ratio <- ratiov[type]
    des <- colData(ids)
    table <- data.frame()
    isoList <- metadata(ids)$isoList
    for (sample in row.names(des)){
        .info <- isoList[[sample]][[coden]]
        temp <- as.data.frame( isoList[[sample]][[coden]] %>%
                                mutate(change=paste0(reference, ">", current)) %>%
                                filter(size==1) %>%
                                group_by(change) %>%
                                summarise( freq=sum(freq), times=n() )
                              )
        total <- sum(temp$freq)
        Total <- sum(temp$times)
        temp$abundance <- temp$freq / total
        temp$unique <- temp$times / Total
        table <- rbind( table,
                        data.frame( change=temp$change, abundance=temp$abundance,
                                        unique=temp$unique,
                                        sample=rep(sample, nrow(temp)),
                                        group=rep(des[sample, column],
                                                nrow(temp)) ) )
    }

    p <- ggplot(table) +
        geom_jitter(aes(x=factor(change),y=unique,colour=factor(group),
                        size=abundance)) +
        scale_colour_brewer("Groups",palette="Set1") +
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3")) +
        labs(list(title=paste(type,"distribution"),y="# of unique sequences",
                x=paste0("changes at postiion ",position," respect to the reference")))
    print(p)
}

#' Create count matrix with different summarizing options
#'
#' This function collapses isomiRs into different groups. It is a similar
#' concept than how to work with gene isoforms. With this function,
#' different changes can be put together into a single miRNA variant.
#' For instance all sequences with variants at 3' end can be
#' considered as different elements in the table
#' or analysis having the following naming
#' \code{hsa-miR-124a-5p.iso.t3:AAA}.
#'
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param ref differentiate reference miRNA from rest
#' @param iso5 differentiate trimming at 5 miRNA from rest
#' @param iso3 differentiate trimming at 3 miRNA from rest
#' @param add differentiate additions miRNA from rest
#' @param subs differentiate nt substitution miRNA from rest
#' @param seed differentiate changes in 2-7 nts from rest
#' @param minc int minimum number of isomiR sequences to be included.
#' @param mins int minimum number of samples with number of 
#' sequences bigger than \code{minc} counts.
#'
#' @details
#'
#' You can merge all isomiRs into miRNAs by calling the function only
#' with the first parameter \code{isoCounts(ids)}.
#' You can get a table with isomiRs altogether and
#' the reference miRBase sequences by calling the function with \code{ref=TRUE}.
#' You can get a table with 5' trimming isomiRS, miRBase reference and
#' the rest by calling with \code{isoCounts(ids, ref=TRUE, iso5=TRUE)}.
#' If you set up all parameters to TRUE, you will get a table for
#' each different sequence mapping to a miRNA (i.e. all isomiRs).
#' 
#' Examples for the naming used for the isomiRs are at
#' \url{http://seqcluster.readthedocs.org/mirna_annotation.html#mirna-annotation}.
#'
#' @return \code{\link{IsomirDataSeq}} object with new count table.
#' The count matrix can be access with \code{counts(ids)}.
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, ref=TRUE)
#' head(counts(ids))
#' # taking into account isomiRs and reference sequence.
#' ids <- isoCounts(mirData, ref=TRUE, minc=10, mins=6)
#' head(counts(ids))
#' @export
isoCounts <- function(ids, ref=FALSE, iso5=FALSE, iso3=FALSE,
                      add=FALSE, subs=FALSE, seed=FALSE, minc=1, mins=1){
        counts <- IsoCountsFromMatrix(metadata(ids)$rawList, colData(ids), ref,
                                      iso5, iso3,
                                      add, subs, seed)
        counts <- counts[rowSums(counts > minc) >= mins, ]
        se <- SummarizedExperiment(assays = SimpleList(counts=counts),
                                   colData = colData(ids))
        .IsomirDataSeq(se, metadata(ids)$rawList, metadata(ids)$isoList)
}


#' Normalize count matrix
#'
#' This function normalizes raw count matrix using
#' \code{\link[DESeq2]{rlog}} function from \code{\link[DESeq2]{DESeq2-package}}.
#'
#' @param ids object of class \code{\link{IsomirDataSeq}}
#' @param formula formula that will be used for normalization
#'
#' @return \code{\link{IsomirDataSeq}} object with the normalized
#' count matrix in a slot. The normalized matrix
#' can be access with \code{counts(ids, norm=TRUE)}.
#'
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, minc=10, mins=6)
#' ids <- isoNorm(mirData, formula=~condition)
#' head(counts(ids, norm=TRUE))
#' @export
isoNorm <- function(ids, formula=~condition){
    dds <- DESeqDataSetFromMatrix(countData = counts(ids),
                                colData = colData(ids),
                                design = formula)
    rld <- rlog(dds, blind=FALSE)
    normcounts(ids) <- assay(rld)
    ids
}
