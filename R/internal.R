.make_uid <- function(rawData){
    is_subs = rawData[["mism"]] != "0"
    is_add = rawData[["add"]] != "0"
    is_t5 = rawData[["t5"]] != "0"
    is_t3 = rawData[["t3"]] != "0"
    is_ref = rawData[["mism"]] == "0" & rawData[["add"]] == "0" & rawData[["t5"]] == "0" & rawData[["t3"]] == "0"
    rawData %>%
        mutate(uid = mir) %>%
        mutate(uid = ifelse(is_ref,
                            paste0(uid, paste0(";ref")),
                            uid)) %>%
        mutate(uid = ifelse(is_subs,
                            paste0(uid, paste0(";iso_snp:", mism)),
                            uid)) %>%
        mutate(uid = ifelse(is_add,
                            paste0(uid, paste0(";iso_add:", add)),
                            uid)) %>%
        mutate(uid = ifelse(is_t5,
                            paste0(uid, paste0(";iso_5p:", t5)),
                            uid)) %>%
        mutate(uid = ifelse(is_t3,
                            paste0(uid, paste0(";iso_3p:", t3)),
                            uid))
}


.remove_gt_n_changes <- function(iso, n=1){
    iso %>%
        mutate(changes = sapply(mism,
                                function(x) str_count(x, "[0-9]+"))) %>%
        mutate(changes = ifelse(mism == "0", 0, changes)) %>%
        filter(changes <= n)
}

.clean_noise <- function(iso, pctco=0, whitelist=NULL){
    sample <- mir <- value <- seq <- NULL
    prop <- fdr <- p.value <- NULL
    whitelist <- intersect(iso[["seq"]], whitelist)
    if (pctco==0){
        keep <- iso[["seq"]] %>% unique()
    }else{
        keep <- iso %>%  # calculate pct
            .[,c(1:2,7:ncol(.))] %>%
            gather("sample", "value", -mir, -seq) %>%
            group_by(sample, mir) %>%
            filter(value > 0) %>%
            group_by(mir, seq) %>%
            summarise(value=sum(value)) %>%
            arrange(mir, desc(value)) %>%
            mutate(rank = 1:n(),
                   total = sum(value),
                   pct = value / total * 100) %>% # statistically calculate if pct  > 10%
            # filter(pct > pctco * 100) %>%
            rowwise %>%
            mutate(prop = list(tidy(prop.test(value,
                                              total,
                                              pctco,
                                              alternative = "greater")))) %>%
            unnest(prop) %>%
            group_by(seq) %>%
            mutate(hits = n()) %>% # remove pct < 10% after p.adjust correction
            ungroup() %>%
            filter(hits == 1) %>% # only uniquely mapped reads
            mutate(fdr = p.adjust(p.value, method = "BH")) %>%
            filter(fdr < 0.05) %>% .[["seq"]] %>% unique()
    }
    iso[iso[["seq"]]  %in%  unique(c(keep, whitelist)),]
}


# put header to input files
.put_header <- function(table){
    if ( sum(colnames(table) == "seq")==0 ){
        names(table)[c(1, 3, 4, 7, 8, 9, 10, 13, ncol(table))] <- c("seq", "freq", "mir",
                                                                    "mism", "add", "t5", "t3",
                                                                    "DB", "ambiguity")
    }
    table <- table[,c("seq", "freq", "mir",
                      "mism", "add", "t5", "t3",
                      "DB", "ambiguity")]
    table[,2] <- as.numeric(table[,2])
    return(table)
}


.filter_N_add <- function(table){
    table[!grepl("N", table[["add"]]),]
}

.change_seq <- function(x, mism){
    .pos = gsub("[ATGCNU]", "", mism)
    .subs = as.vector(unlist(strsplit(gsub("[0-9]+", "", mism), "")))
    .nts = as.vector(unlist(strsplit(x, "")))
    .nts[as.numeric(.pos)] = .subs[2]
    paste0(as.vector(unlist(.nts)), collapse = "")
}

.clean_low_rate_changes <- function(tab, rate=0.20, uniqueMism=TRUE){
    if (uniqueMism){
        tab = tab  %>%
            filter(!(mism != "0" & ambiguity > 1))
    }
    tab.fil = tab %>%
        rowwise() %>%
        mutate(seq = if_else((af<rate & !is.na(af)) | grepl("N", mism),
                             .change_seq(seq, mism), seq)) %>%
        mutate(mism = if_else((af<rate & !is.na(af) | grepl("N", mism)),
                              "0", mism))
    tab.fil = tab.fil %>% ungroup() %>%
        group_by(mir, seq, mism, add, t5, t3, DB, ambiguity) %>%
        summarise(freq = sum(freq))  %>%
        dplyr::select(mir, seq, freq, mism, add, t5, t3, DB, ambiguity)
    as.data.frame(tab.fil)
}

# filter by relative abundance to reference
.filter_by_cov <- function(table, limit=0, rate=0.2,
                           canonicalAdd=TRUE, uniqueMism=TRUE){
    freq <- mir <-  NULL
    if (canonicalAdd){
        tab.fil <- table %>% filter(DB == "miRNA",
                                    !grepl("[GC]", add))
    }else{
        tab.fil <- table %>% filter(DB == "miRNA")
    }
    tab.fil.out <- as.data.frame(tab.fil %>% filter(mism==0,
                                                    nchar(add)<3) %>%
                                   group_by(mir) %>%
                                   summarise(mir_f=sum(freq)+1,
                                             mir_n=n()+1))
    tab.fil <- tab.fil %>%
        left_join(tab.fil.out, by="mir")
    tab.mism <- tab.fil %>% filter(mism!=0)
    if (nrow(tab.mism) == 0)
        return(tab.fil)

    tab.mism <- tab.mism %>%
        group_by(mir, mism, mir_n, mir_f) %>%
        summarise(mism_n=n(), mism_f=sum(freq))  %>%
        mutate(enrich=mism_n/mir_n, af=mism_f/mir_f, bias=af/enrich) %>%
        ungroup() %>%
        select(-mir_n, -mir_f)

    tab.fil <- left_join(tab.fil %>% mutate(id=paste(mir,mism)),
                     tab.mism %>% mutate(id=paste(mir,mism)) %>%
                         select(-mism, -mir),
                     by="id") %>% select(-id)
    tab.fil$score <- tab.fil$freq / tab.fil$mir_f * 100
    tab.fil <- .clean_low_rate_changes(tab.fil, rate, uniqueMism)

    tab.fil <- tab.fil %>%
        left_join(tab.fil.out, by="mir") %>%
        mutate(id=paste(mir,mism)) %>%
        left_join(tab.mism %>% mutate(id=paste(mir,mism)) %>%
                      select(-mism, -mir),
                  by="id") %>% select(-id)
    tab.fil$score <- tab.fil$freq / tab.fil$mir_f * 100
    tab.fil
}

.convert_to_new_version <- function(table){
    idx <- grepl("u-", table$t5)
    table$t5[idx] <- toupper(gsub("u-", "", table$t5[idx]))
    idx <- grepl("d-", table$t5)
    table$t5[idx] <- tolower(gsub("d-", "", table$t5[idx]))
    idx <- grepl("u-", table$t3)
    table$t3[idx] <- tolower(gsub("u-", "", table$t3[idx]))
    idx <- grepl("d-", table$t3)
    table$t3[idx] <- toupper(gsub("d-", "", table$t3[idx]))
    table$add <- toupper(gsub("u-", "", table$add))
    table
}

# Filter table reference
.filter_table <- function(table, cov=1, rate=0.2,
                          canonicalAdd=TRUE, uniqueMism=TRUE,
                          uniqueHits = FALSE){
    table <- .put_header(table)
    if (uniqueHits)
        table <- table[table[["ambiguity"]] == 1, ]
    table <- .filter_N_add(table)
    table <- .filter_by_cov(table, cov, rate, canonicalAdd, uniqueMism)
    if (sum(grepl("u-", table$add))>0)
        table <- .convert_to_new_version(table)
    table
}


# plot general information
.isomir_general_type <- function(table, colid){
    temp <- table
    temp$idfeat <- paste(table[ ,colid], table$mir)
    temp <- temp[order(temp$idfeat), ]
    temp <- temp[!duplicated(temp$idfeat), ]
    temp <- as.data.frame(summary(temp$mir))
    feat.dist <- cut(as.numeric(temp[ ,1]), breaks=c(-1, 0.5, 1.5, 2.5, Inf),
                   labels=c("0", "1", "2", ">3"))
    return ( as.data.frame( summary(feat.dist) ) )
}

# do counts table considering what isomiRs take into account
IsoCountsFromMatrix <- function(rawData, des, ref=FALSE, iso5=FALSE,
                                iso3=FALSE, add=FALSE,
                                snv=FALSE, seed=FALSE, minc=1){
    is_subs = snv & rawData[["mism"]] != "0"
    is_add = add & rawData[["add"]] != "0"
    is_t5 = iso5 & rawData[["t5"]] != "0"
    is_t3 = iso3 & rawData[["t3"]] != "0"
    is_ref = ref & rawData[["mism"]] == "0" & rawData[["add"]] == "0" & rawData[["t5"]] == "0" & rawData[["t3"]] == "0"
    dt <- rawData %>%
        mutate(uid = mir) %>%
        mutate(uid = ifelse(is_ref,
                            paste0(uid, paste0(";ref")),
                            uid)) %>%
        mutate(uid = ifelse(is_subs,
                            paste0(uid, paste0(";iso_snp:", mism)),
                            uid)) %>%
        mutate(uid = ifelse(is_add,
                            paste0(uid, paste0(";iso_add:", add)),
                            uid)) %>%
        mutate(uid = ifelse(is_t5,
                            paste0(uid, paste0(";iso_5p:", t5)),
                            uid)) %>%
        mutate(uid = ifelse(is_t3,
                            paste0(uid, paste0(";iso_3p:", t3)),
                            uid)) %>%
        .[,c("uid", rownames(des))] %>%
        group_by(!!sym("uid")) %>%
        summarise_all(sum) %>%
        as.data.frame() %>%
        remove_rownames() %>%
        column_to_rownames("uid") %>%
        as.matrix()
    if (dim(dt)[1] == 0)
        warning("No miRNA found. Make sure the third column of the file has the count value different than 0.")
    dt
}

# merge matrix by column in metadata
.merge_counts_by <- function(counts, coldata, merge_by){
    uni_counts <-lapply(unique(as.character(coldata[[merge_by]])),
                       function(s){
        one <- as.character(rownames(coldata)[coldata[[merge_by]]==s])
        round(rowMeans(counts[,one, drop=FALSE]))
    }) %>% unlist() %>%
        matrix(.,
               ncol=length(unique(as.character(coldata[[merge_by]]))))
    colnames(uni_counts) <- unique(as.character(coldata[[merge_by]]))
    rownames(uni_counts) <- rownames(counts)
    coldata = coldata %>%
        as.data.frame() %>%
        distinct(!!!sym(merge_by), .keep_all=TRUE) %>%
        DataFrame()
    rownames(coldata) <- coldata[[merge_by]]
    list(counts=uni_counts, coldata=coldata)
}

# get isomir name and sequence table
.make_isomir_naming <- function(rawData){
    is_subs = rawData[["mism"]] != "0"
    is_add = rawData[["add"]] != "0"
    is_t5 = rawData[["t5"]] != "0"
    is_t3 = rawData[["t3"]] != "0"
    is_ref = rawData[["mism"]] == "0" & rawData[["add"]] == "0" & rawData[["t5"]] == "0" & rawData[["t3"]] == "0"
    dt <- rawData %>%
        mutate(isomir = mir) %>%
        mutate(isomir = ifelse(is_ref,
                            paste0(isomir, paste0(";ref")),
                            isomir)) %>%
        mutate(isomir = ifelse(is_subs,
                            paste0(isomir, paste0(";iso_snp:", mism)),
                            isomir)) %>%
        mutate(isomir = ifelse(is_add,
                            paste0(isomir, paste0(";iso_add:", add)),
                            isomir)) %>%
        mutate(isomir = ifelse(is_t5,
                            paste0(isomir, paste0(";iso_5p:", t5)),
                            isomir)) %>%
        mutate(isomir = ifelse(is_t3,
                            paste0(isomir, paste0(";iso_3p:", t3)),
                            isomir))
    return(dt[,c("seq", "isomir")])
}

# Collapse isomiRs in miRNAs
.collapse_mirs <- function(table, ref=FALSE, iso5=FALSE, iso3=FALSE,
                        add=FALSE, snv=FALSE, seed=FALSE){
    label <- table$mir
    freq <- id <- NULL
    if (ref == TRUE){
        ref.val <- do.call(paste, table[,4:7])
        ref.val[grep("[ATGC]", ref.val, invert=TRUE,
                     ignore.case = TRUE)] <- "ref"
        ref.val[grep("[ATGC]", ref.val, ignore.case = TRUE)] <- "iso"
        label <- paste(label, ref.val, sep=".")
    }
    if (iso5 == TRUE){
        label <- paste(label, table[,"t5"], sep=".t5:")
    }
    if (seed == TRUE){
        seed.val <- as.character(table[,"mism"])
        seed.val[grep("^[2-8][ATGC]", seed.val, invert=TRUE)] <- "0"
        label <- paste(label, seed.val, sep=".seed:")
    }
    if (iso3 == TRUE){
        label <- paste(label, table[,"t3"], sep=".t3:")
    }
    if (add == TRUE){
        label <- paste(label, table[,"add"], sep=".ad:")
    }
    if (snv == TRUE){
        label <- paste(label, table[,"mism"], sep=".mm:")
    }

    table$id <- label
    table.out <- as.data.frame(table %>% group_by(id) %>%
                                 summarise(total=sum(freq)))
    table.out[is.na(table.out)] <- 0
    table.out
}

# Do summary of different isomiRs events
.isomir_position <- function(column){
    size <- sapply(column, function(x){
        if (x == "0")
            return(0)
        p <- nchar(x)
        if (grepl("[atgcn]", x))
            p <- p * -1L
        return(p)
    }) %>% unlist() %>%
        as.vector()
   size
}

# Do summary of nt substitution events
.subs_position <- function(column){
    column
    nt <- sub("[0-9]+", "", column)
    pos <- sub("[ATGC]{2}", "", column)
    pos <- data.frame(nt = as.character(nt), size = pos,
                      stringsAsFactors = FALSE) %>%
        separate(nt, sep = 1, into = c("isomir", "reference")) %>%
        unite(change, reference, isomir, sep = "->")
    pos[["change"]][pos[["size"]] == "0"] <- "0"
    pos[["size"]] <- factor(pos[["size"]],
                            levels = sort(as.numeric(unique(pos[["size"]]))))
    return(pos)
}

# function to plot isomiR summary
.plot_all_iso <- function(ids, column, use){
    des <- colData(ids) %>%
        as.data.frame() %>%
        rownames_to_column("iso_sample")

    stopifnot(column  %in% colnames(des))

    rawData <- metadata(ids)[["rawData"]]
    if (!is.null(use)){
        rawData <- .make_uid(rawData)
        rawData <- rawData[rawData[["uid"]]  %in% use,]
    }
    message("Using ", nrow(rawData), " isomiRs.")
    if (nrow(rawData) == 0)
        stop("Any of the `use` elements is in the data set.")

    is_subs <- rawData[["mism"]] != "0"
    is_add <- rawData[["add"]] != "0"
    is_t5 <- rawData[["t5"]] != "0"
    is_t3 <- rawData[["t3"]] != "0"
    is_ref <- rawData[["mism"]] == "0" & rawData[["add"]] == "0" & rawData[["t5"]] == "0" & rawData[["t3"]] == "0"

    iso_data <- rawData %>%
        mutate(uid = "") %>%
        mutate(uid = ifelse(is_ref,
                            paste0(uid, ";ref"),
                            uid)) %>%
        mutate(uid = ifelse(is_subs,
                            paste0(uid, ";iso_snp"),
                            uid)) %>%
        mutate(uid = ifelse(is_add,
                            paste0(uid, ";iso_add"),
                            uid)) %>%
        mutate(uid = ifelse(is_t5,
                            paste0(uid, ";iso_5p"),
                            uid)) %>%
        mutate(uid = ifelse(is_t3,
                            paste0(uid, ";iso_3p"),
                            uid)) %>%
        .[,c("uid", des[["iso_sample"]])] %>%
        separate_rows(!!sym("uid"), sep = ";") %>%
        filter(!!sym("uid") != "")

    freq_total <- data.frame(total_sum = colSums(rawData[, setdiff(colnames(iso_data), "uid")]),
                             iso_sample = colnames(ids),
                             stringsAsFactors = F)

    freq_data <- iso_data %>%
        group_by(!!sym("uid")) %>%
        summarise_all(sum) %>%
        ungroup() %>%
        gather("iso_sample", "sum", -!!sym("uid"))

    n_total <- data.frame(total_sum = colSums(rawData[, setdiff(colnames(iso_data), "uid")]>0),
                          iso_sample = colnames(ids),
                          stringsAsFactors = F)

    n_data <- iso_data %>%
        group_by(!!sym("uid")) %>%
        summarise_all(list(~sum(. > 0))) %>%
        ungroup() %>%
        gather("iso_sample", "sum", -!!sym("uid"))

    freq_pct <-
        left_join(freq_total, freq_data, by = "iso_sample") %>%
        mutate(pct_abundance = sum / total_sum * 100L,
               method = "isomiRs abundance") %>%
        .[,c("iso_sample", "method", "uid", "pct_abundance")]

    n_pct <-
        left_join(n_total, n_data, by = "iso_sample") %>%
        mutate(pct_abundance = sum / total_sum *100L,
               method = "isomiRs") %>%
        .[,c("iso_sample", "method", "uid", "pct_abundance")]

    bind_rows(freq_pct, n_pct) %>%
        left_join(des, by ="iso_sample") %>%
        arrange(uid) %>%
        ggplot(aes_string(x="uid", y="pct_abundance",
                          group="iso_sample", color = column)) +
        geom_polygon(fill=NA) +
        coord_polar(start=-pi) +
        scale_color_brewer(palette = "Set1") +
        facet_wrap(~method) +
        theme_bw()
}