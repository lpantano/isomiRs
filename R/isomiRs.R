#' Differential expression analysis with DESeq2
#'
#' This function does differential expression analysis with
#' [DESeq2::DESeq2-package] using the specific formula.
#' It will return a [DESeq2::DESeqDataSet] object.
#'
#' @details
#'
#' First, this function collapses all isomiRs in different types.
#' Read more at [isoCounts()] to know the different options
#' available to collapse isomiRs.
#'
#' After that, [DESeq2::DESeq2-package] is used to do differential
#' expression analysis. It uses the count matrix and design experiment
#' stored at (`counts(ids)` and `colData(ids)`)
#' [IsomirDataSeq] object
#' to construct a [DESeq2::DESeqDataSet] object.
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param formula Formula used for DE analysis.
#' @param ... Options to pass to [isoCounts()] including
#'   ref, iso5, iso3, add, subs and seed parameters.
#'
#' @return [DESeq2::DESeqDataSet] object.
#' To get the differential expression isomiRs, use [DESeq2::results()] from
#' DESeq2 package. This allows to ask for different contrast
#' without calling again [isoDE()]. Read `results`
#' manual to know how to access all the information.
#'
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, minc=10, mins=6)
#' dds <- isoDE(mirData, formula=~condition)
#' @export
isoDE <- function(ids, formula=NULL, ...){
    if (is.null(formula)){
        formula <- design(ids)
    }
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
#' isomiRs/miRNAs. It uses the matrix under `counts(ids)`
#' to get the top expressed isomiRs/miRNAs using the average
#' expression value
#' and plot a heatmap with the raw counts for each sample.
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param top Number of isomiRs/miRNAs used.
#' 
#' @examples
#' data(mirData)
#' isoTop(mirData)
#' @return PCA of the top expressed miRNAs
#' @export
isoTop <- function(ids, top=20, condition="condition"){
    select <- order(rowMeans(counts(ids)),
                    decreasing=TRUE)[1:top]
    degPCA(counts(ids)[select,], metadata=colData(ids), condition=condition)
}

#' Plot the amount of isomiRs in different samples
#'
#' This function plot different isomiRs proportion for each sample.
#' It can show trimming events at both side, additions and nucleotides
#' changes.
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param type String (iso5, iso3, add, snv, all) to indicate what isomiRs
#'   to use for the plot. See details for explanation.
#' @param column String indicating the column in
#'   `colData` to color samples.
#' @param use Character vector to only use these isomiRs for
#'   the plot. The id used is the rownames that comes from using
#'   isoCounts with all the arguments on TRUE.
#' @param nts Boolean to indicate whether plot positions of nucleotides
#'   changes when showing single nucleotides variants.
#' @return [ggplot2::ggplot()] Object showing different isomiRs changes at
#' different positions.
#' @details
#' There are four different values for `type` parameter. To plot
#' trimming at 5' or 3' end, use `type="iso5"` or `type="iso3"`. Get a summary of all using `type="all"`.
#' In this case, it will plot 3 positions at both side of the reference
#' position described at miRBase site. Each position refers to the % of
#' sequences that start/end before or after the miRBase reference. The
#' color indicates the sample group. The size of the point is proportional
#' to the abundance considering the total as all the sequences in the sample.
#' The position at `y` is the % of
#' different sequences considering the total as all sequences with changes
#' for the specific
#' isomiR showed.
#'
#' Same logic applies to `type="add"` and `type="subs"`. However,
#' when `type="add"`, the plot will refer to addition events from the
#' 3' end of the reference position. Note that this additions don't match
#' to the precursor sequence, they are non-template additions.
#' In this case, only 3 positions after the 3' end
#' will appear in the plot. When `type="subs"`, it will appear one
#' position for each nucleotide in the reference miRNA. Points
#' will indicate isomiRs with nucleotide changes at the given position.
#' When `type="all"` a colar coordinate map will show 
#' the abundance of each isomiR type in a single plot. 
#' Note the position is relatively to the
#' sequence not the miRNA.
#'
#' @examples
#' data(mirData)
#' isoPlot(mirData)
#' @export
isoPlot <- function(ids, type="iso5", column=NULL,
                    use = NULL, nts = FALSE){
    
    if (is.null(column)){
        column <-  names(colData(ids))[1]
    }
    supported <- c("snv", "add", "iso5", "iso3")
    if (type == "all"){return(.plot_all_iso(ids, column, use))}
    if (type=="subs"){
        type <- "snv"
    }
    stopifnot(type  %in% c("snv", "add", "iso5", "iso3"))
    
    freq <- size <- group <- abundance <- NULL
    codevn <- 3:6

    names(codevn) <- supported
    coden <- codevn[type]
    des <- colData(ids) %>% 
        as.data.frame() %>% 
        rownames_to_column("iso_sample")
    rawData <- metadata(ids)[["rawData"]]
    if (!is.null(use)){
        rawData <- .make_uid(rawData)
        rawData <- rawData[rawData[["uid"]]  %in% use,]
    }
    message("Using ", nrow(rawData), " isomiRs.")
    if (nrow(rawData) == 0)
        stop("Any of the `use` elements is in the data set.")

    if (type == "snv"){
        get_column <- "size"
        xaxis <- "position of the isomiR"
        if (nts){
            get_column <- "change"
            xaxis <- "nucleotide change"
        }
        rawData[["size"]] <- .subs_position(rawData[[coden]])[[get_column]]

    }else{
        rawData[["size"]] <- as.factor(.isomir_position(rawData[[coden]]))
        xaxis <- "position respect to the reference"
    }
    freq_data <- rawData[,c("size", des[["iso_sample"]])] %>% 
        group_by(!!sym("size")) %>% 
        summarise_all(funs(sum)) %>% 
        ungroup() %>% 
        gather("iso_sample", "sum", -!!sym("size"))
    
    n_data <- rawData[,c("size", des[["iso_sample"]])] %>%
        group_by(!!sym("size")) %>% 
        summarise_all(funs(sum(. > 0))) %>% 
        ungroup() %>% 
        gather("iso_sample", "sum", -!!sym("size"))
    
    freq_pct <- freq_data %>% 
        group_by(!!sym("iso_sample")) %>% 
        summarise(total_sum = sum(sum)) %>%
        left_join(freq_data, by = "iso_sample") %>% 
        mutate(pct_abundance = sum / total_sum * 100L,
               id = paste(iso_sample, size)) %>% 
        .[,c("id", "iso_sample", "size", "pct_abundance")]
    
    n_pct <- n_data %>% 
        group_by(!!sym("iso_sample")) %>% 
        summarise(total_sum = sum(sum)) %>%
        left_join(n_data, by = "iso_sample") %>% 
        mutate(unique = sum / total_sum *100L,
               id = paste(iso_sample, size)) %>% 
        .[,c("id", "unique")]
    
    inner_join(freq_pct, n_pct, by="id") %>%
        filter(size!="0") %>% 
        left_join(des, by ="iso_sample") %>% 
        ggplot() +
        geom_jitter(aes_string(x="size",y="unique", colour=column,
                        size="pct_abundance")) +
        scale_colour_brewer("Groups",palette="Set1") +
        theme_bw(base_size = 14, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3")) +
        labs(list(title=paste(type,"distribution"),
                  y="pct of isomiRs",
                x=xaxis))
}

#' Plot nucleotides changes at a given position
#'
#' This function plot different isomiRs proportion for each sample at a given
#' position focused on the nucleotide change that happens there.
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param position Integer indicating the position to show.
#' @param column String indicating the column in
#'   colData to color samples.
#' @return [ggplot2::ggplot()] Object showing nucleotide changes
#' at a given position.
#' @details
#' It shows the nucleotides changes at the given position for each
#' sample in each group.
#' The color indicates the sample group. The size of the point is proportional
#' to the number of total counts of isomiRs with changes.
#' The position at `y` is the % of different isomiRs
#' supporting the change. Note the position is relatively to the
#' sequence not the miRNA.
#'
#'
#' @examples
#' data(mirData)
#' isoPlotPosition(mirData)
#' @export
isoPlotPosition <- function(ids, position = 1L, column = NULL){
    if (is.null(column)){
        column <-  names(colData(ids))[1]
    }

    des <- colData(ids) %>% 
        as.data.frame() %>% 
        rownames_to_column("iso_sample")
    rawData <- metadata(ids)[["rawData"]]
    
    parsed_change <- .subs_position(rawData[["mism"]])
    rawData[["change"]] <- as.factor(parsed_change[["change"]])
    rawData[["pos"]] <- as.character(parsed_change[["size"]])

    freq_data <- rawData[,c("change", "pos", des[["iso_sample"]])] %>% 
        group_by(!!sym("change"), !!sym("pos")) %>% 
        summarise_all(funs(sum)) %>% 
        ungroup() %>% 
        gather(iso_sample, sum, -change, -pos)
    
    n_data <- rawData[,c("change", "pos", des[["iso_sample"]])] %>%
        group_by(!!sym("change"), !!sym("pos")) %>% 
        summarise_all(funs(sum(. > 0))) %>% 
        ungroup() %>% 
        gather(iso_sample, sum, -change, -pos)
    
    freq_pct <- freq_data %>% 
        group_by(!!sym("iso_sample")) %>% 
        summarise(total_sum = sum(sum)) %>%
        left_join(freq_data, by = "iso_sample", suffix = c("", "_tmp")) %>% 
        mutate(pct_abundance = sum / total_sum * 100L,
               id = paste(iso_sample, change)) %>% 
        .[.[["pos"]] == position,] %>% 
        .[,c("id", "iso_sample", "change", "pct_abundance")]
    
    n_pct <- n_data %>% 
        group_by(!!sym("iso_sample")) %>% 
        summarise(total_sum = sum(sum)) %>%
        left_join(n_data, by = "iso_sample") %>% 
        mutate(unique = sum / total_sum *100L,
               id = paste(iso_sample, change)) %>% 
        .[.[["pos"]] == position,] %>% 
        .[,c("id", "unique")]
    
    inner_join(freq_pct, n_pct, by="id") %>%
        left_join(des, by ="iso_sample") %>%
        ggplot() +
        geom_jitter(aes_string(x="change",y="unique",colour=column,
                        size="pct_abundance")) +
        scale_colour_brewer("Groups",palette="Set1") +
        theme_bw(base_size = 11, base_family = "") +
        theme(strip.background=element_rect(fill="slategray3")) +
        labs(list(title=paste("Change distribution"),
                  y="pct of isomiRs",
                x=paste0("changes at postiion ",
                         position, 
                         " respect to the reference")))
}

#' Create count matrix with different summarizing options
#'
#' This function collapses isomiRs into different groups. It is a similar
#' concept than how to work with gene isoforms. With this function,
#' different changes can be put together into a single miRNA variant.
#' For instance all sequences with variants at 3' end can be
#' considered as different elements in the table
#' or analysis having the following naming
#' `hsa-miR-124a-5p.iso.t3:AAA`.
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param ref Differentiate reference miRNA from rest.
#' @param iso5 Differentiate trimming at 5 miRNA from rest.
#' @param iso3 Differentiate trimming at 3 miRNA from rest.
#' @param add Differentiate additions miRNA from rest.
#' @param snv Differentiate nt substitution miRNA from rest.
#' @param all Differentiate all isomiRs.
#' @param seed Differentiate changes in 2-7 nts from rest.
#' @param minc Int minimum number of isomiR sequences to be included.
#' @param mins Int minimum number of samples with number of
#'   sequences bigger than `minc` counts.
#' @param merge_by Column in coldata to merge samples into a single
#'   column in counts. Useful to combine technical replicates.
#'
#' @details
#'
#' You can merge all isomiRs into miRNAs by calling the function only
#' with the first parameter `isoCounts(ids)`.
#' You can get a table with isomiRs altogether and
#' the reference miRBase sequences by calling the function with `ref=TRUE`.
#' You can get a table with 5' trimming isomiRS, miRBase reference and
#' the rest by calling with `isoCounts(ids, ref=TRUE, iso5=TRUE)`.
#' If you set up all parameters to TRUE, you will get a table for
#' each different sequence mapping to a miRNA (i.e. all isomiRs). 
#'
#' Examples for the naming used for the isomiRs are at
#' http://seqcluster.readthedocs.org/mirna_annotation.html#mirna-annotation.
#'
#' @return [IsomirDataSeq] object with new count table.
#' The count matrix can be access with `counts(ids)`.
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, ref=TRUE)
#' head(counts(ids))
#' # taking into account isomiRs and reference sequence.
#' ids <- isoCounts(mirData, ref=TRUE, minc=10, mins=6)
#' head(counts(ids))
#' @export
isoCounts <- function(ids, ref=FALSE, iso5=FALSE, iso3=FALSE,
                      add=FALSE, snv=FALSE, seed=FALSE, 
                      all = FALSE,minc=1, mins=1,
                      merge_by=NULL){
    coldata <- colData(ids)
    rawdata <- metadata(ids)[["rawData"]]
    counts <- IsoCountsFromMatrix(metadata(ids)[["rawData"]],
                                  colData(ids),
                                  ref,
                                  iso5, iso3,
                                  add, snv, seed)
    if (!is.null(merge_by)){
        stopifnot(merge_by  %in% colnames(colData(ids)))
        combined <- .merge_counts_by(counts, colData(ids), merge_by)
        counts <- combined[["counts"]]
        coldata <-  combined[["coldata"]]
        
        raw_counts <- rawdata[,7:ncol(rawdata)] %>% as.matrix()
        new <- .merge_counts_by(raw_counts, colData(ids), merge_by)
        rawdata <- as.tibble(cbind(rawdata[,1:6], new[["counts"]]))
        
    }

    counts <- counts[rowSums(counts > minc) >= mins, ]
    se <- SummarizedExperiment(assays = SimpleList(counts = counts),
                               colData = coldata)
    .IsomirDataSeq(se, rawdata)
}

#' Annotate the rawData of the [IsomirDataSeq] object
#' 
#' Get the sequence and the name information for each isomiR,
#' and the importance value (`isomir_reads/mirna_reads`) for
#' each  sample.
#' 
#' @param ids Object of class [IsomirDataSeq].
#' 
#' @details 
#' `edit_mature_position` represents the position at the mature
#'  sequence + nucleotide at reference + nucleotide at isomiR.
#' @return [data.frame] with the sequence, isomir name,
#'   and importance for each sample and isomiR.
#' 
#' @examples
#' data(mirData)
#' head(isoAnnotate(mirData))
#' @export
isoAnnotate <- function(ids){
    sample <- mir <- value <- uid <- seq <- NULL
    rawData <- metadata(ids)[["rawData"]]
    rawData <- .make_uid(rawData)
    dt <- rawData %>%  # calculate pct
        .[,c(1:2,7:ncol(.))] %>%
        gather("sample", "value", -mir, -seq, -uid) %>%
        group_by(sample, mir) %>% 
        filter(value > 0) %>%
        group_by(mir, sample, seq, uid) %>% 
        summarise(value=sum(value)) %>% 
        group_by(mir, sample) %>% 
        arrange(sample, mir, desc(value)) %>%
        mutate(pct = value / sum(value) * 100) %>%
        ungroup() %>% 
        select(-value, -mir) %>% 
        spread(sample, pct) %>% 
        as.data.frame()
    
    lift <- ((grepl("[A-Z]", rawData[["t5"]]) * -1) + (grepl("[a-z]", rawData[["t5"]]) * 1)) * nchar(rawData[["t5"]])
    pos <- as.numeric(str_extract(rawData[["mism"]], "[0-9]*")) + lift
    pos <- ifelse(pos==0, -1, pos)
    change <- reverse(str_extract(rawData[["mism"]], "[ACTG]+"))
    mature_pos <- ifelse(grepl("NA", change), "NA", paste0(pos, ":", change))
    
    dt[["edit_mature_position"]] = mature_pos
    dt
}


#' Normalize count matrix
#'
#' This function normalizes raw count matrix using
#' [DESeq2::rlog()] function from [DESeq2::DESeq2-package].
#'
#' @param ids Object of class [IsomirDataSeq].
#' @param formula Formula that will be used for normalization.
#' @param maxSamples Maximum number of samples to use with
#'   [DESeq2::rlog()], if not [limma::voom()] is used.
#' @return [IsomirDataSeq] object with the normalized
#' count matrix in a slot. The normalized matrix
#' can be access with `counts(ids, norm=TRUE)`.
#'
#' @examples
#' data(mirData)
#' ids <- isoCounts(mirData, minc=10, mins=6)
#' ids <- isoNorm(mirData, formula=~condition)
#' head(counts(ids, norm=TRUE))
#' @export
isoNorm <- function(ids, formula=NULL, maxSamples = 50){
    if (is.null(formula)){
        formula <- design(ids)
    }
    if (length(colnames(ids)) < maxSamples){
        dds <- DESeqDataSetFromMatrix(countData = counts(ids),
                                      colData = colData(ids),
                                      design = formula)
        rld <- varianceStabilizingTransformation(dds, blind=FALSE)
        normcounts(ids) <- assay(rld)
    }else{
        d <- colData(ids)
        m <- model.matrix(formula, d)
        v <- voom(counts(ids), m)
        normcounts(ids) <- v[["E"]]
    }
   ids
}
