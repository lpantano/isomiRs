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
                           .change_seq(seq, mism),seq)) %>%
        mutate(mism = if_else((af<rate & !is.na(af)), "0", mism))
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
    tab.mism = tab.fil %>% filter(mism!=0) %>% group_by(mir, mism) %>%
        summarise(mism_n=n(), mism_f=sum(freq)) %>%
        left_join(tab.fil.out, by="mir") %>%
        mutate(enrich=mism_n/mir_n, af=mism_f/mir_f, bias=af/enrich) %>%
        ungroup()
    tab.fil <- left_join(tab.fil %>% mutate(id=paste(mir,mism)),
                     tab.mism %>% mutate(id=paste(mir,mism)) %>%
                         select(-mism, -mir),
                     by="id") %>% select(-id)
    tab.fil$score <- tab.fil$freq / tab.fil$mir_f * 100
    tab.fil <- .clean_low_rate_changes(tab.fil, rate, uniqueMism)

    tab.fil <- left_join(tab.fil %>% mutate(id=paste(mir,mism)),
                         tab.mism %>% mutate(id=paste(mir,mism)) %>%
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
                          canonicalAdd=TRUE, uniqueMism=TRUE){
    table <- .put_header(table)
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
IsoCountsFromMatrix <- function(listTable, des, ref=FALSE, iso5=FALSE,
                                iso3=FALSE, add=FALSE,
                                subs=FALSE, seed=FALSE, minc=1){
    table.merge <- data.frame()
    for (sample in row.names(des)){
        d <- listTable[[sample]]
        d <- .collapse_mirs(d, ref=ref, iso5=iso5, iso3=iso3, add=add,
                            subs=subs, seed=seed)
        names(d)[ncol(d)] <- sample
        d <- d[d[,2] >= minc, ]
        if( nrow(table.merge) == 0){
            table.merge <- d
        }else{
            table.merge <- merge( table.merge, d, by=1, all=TRUE )
        }
    }

    row.names(table.merge) <- table.merge[ ,1]
    table.merge <- as.matrix(table.merge[ ,2:ncol(table.merge),drop=FALSE])
    table.merge[is.na(table.merge)] <- 0
    dt <- as.matrix(table.merge)
    if ( dim(dt)[1] == 0 )
        warning("No miRNA found.")
    dt
}

# Collapse isomiRs in miRNAs
.collapse_mirs <- function(table, ref=FALSE, iso5=FALSE, iso3=FALSE,
                        add=FALSE, subs=FALSE, seed=FALSE){
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
    if (subs == TRUE){
        label <- paste(label, table[,"mism"], sep=".mm:")
    }

    table$id <- label
    table.out <- as.data.frame(table %>% group_by(id) %>%
                                 summarise(total=sum(freq)))
    table.out[is.na(table.out)] <- 0
    table.out
}

# Do summary of different isomiRs events
.isomir_position <- function(table, colid){
    temp <- table
    temp[ ,colid] <- as.character(temp[ ,colid])
    pos <- temp[, colid, drop=FALSE]
    row.names(pos) <- 1:nrow(pos)
    pos$mir <- temp$mir
    pos$freq <- temp$freq
    pos <- pos[pos[ ,1] != 0, ]
    if (nrow(pos)==0)
      return(data.frame())
    pos$size <- apply(pos, 1, function(x){
        p <- length(unlist(strsplit(x[1], "")))
        if (grepl("[atgcn]", x[1]))
            p <- p * -1
        return(p)
    })
    pos$idfeat <- paste(pos$size, pos$mir)
    pos <- pos[order(pos$idfeat, abs(pos$size)),]
    return (pos[ ,c(2,3,4)])
}

# Do summary of nt substitution events
.subs_position <- function(table, colid){
    temp <- table
    temp[ ,colid] <- as.character(temp[ ,colid])
    nt <- sub("[0-9]+", "", temp[ ,colid])
    pos <- sub("[ATGC]{2}", "", temp[,colid])
    pos <- data.frame(nt=as.character(nt), size=pos,
                      mir=temp$mir, freq=temp$freq)
    pos$nt <- as.character(pos$nt)
    if (length(unique(pos$nt)) == 1)
      return(data.frame())
    pos <- pos[pos[,1] != "", ]
    nt.2 <- as.data.frame(t(as.data.frame((strsplit(pos$nt, "", fixed=2)))))
    names(nt.2) <- c("current","reference")
    names(nt.2$current) <- ""
    names(nt.2$reference) <- ""
    pos <- cbind(pos, nt.2)
    pos$size <- factor(pos$size, levels=1:25)
    return (pos[ ,c(3,4,2,5,6)])
}

# function to plot isomiR summary
.plot_all_iso <- function(ids, column){
    rawList <- metadata(ids)[[1]]
    metadata <- as.data.frame(colData(ids))
    metadata$sample <- as.character(row.names(metadata))
    metadata$condition <- metadata[,column]
    .l <- lapply(row.names(colData(ids)), function(sample){
        .d <- rawList[[sample]] %>%
            mutate(mismtag=ifelse(mism != "0", "Yes", "No")) %>%
            mutate(addtag=ifelse(add != "0", "Yes", "No")) %>%
            mutate(t5tag=ifelse(t5 != "0", "Yes", "No")) %>%
            mutate(t3tag=ifelse(t3 != "0", "Yes", "No"))
        .s <- data.frame(sample=sample,
                   type=c("mism","add","t5","t3","ref"),
                   method="freq",
                   freq=c(sum(.d$freq[.d$mismtag=="Yes"], na.rm=TRUE),
                     sum(.d$freq[.d$addtag=="Yes"], na.rm=TRUE),
                     sum(.d$freq[.d$t5tag=="Yes"], na.rm=TRUE),
                     sum(.d$freq[.d$t3tag=="Yes"], na.rm=TRUE),
                     sum(.d$freq[.d$mismtag!="Yes" & .d$addtag!="Yes"  & .d$t5tag!="Yes" & .d$t3tag!="Yes"])),
                     stringsAsFactors = FALSE
        )
        .s$freq <- .s$freq/sum(.s$freq) * 100
        .u <- data.frame(sample=sample,
                         type=c("mism","add","t5","t3","ref"),
                         method="unique",
                         freq=c(sum(.d$mismtag=="Yes"),
                                sum(.d$addtag=="Yes"),
                                sum(.d$t5tag=="Yes"),
                                sum(.d$t3tag=="Yes"),
                                sum(.d$mismtag!="Yes" & .d$addtag!="Yes"  & .d$t5tag!="Yes" & .d$t3tag!="Yes")),
                         stringsAsFactors = FALSE
        )
        .u$freq <- .u$freq/sum(.u$freq) * 100
        rbind(.s, .u)
    })
    do.call(rbind, .l) %>% arrange(type) %>%
        left_join(metadata, by="sample") %>%
        ggplot(aes_string(x="type", y="freq", group="sample", color="condition")) +
        geom_polygon(fill=NA) +
        coord_polar(start=-pi) +
        scale_color_brewer(palette = "Set1") +
        facet_wrap(~method) +
        theme_bw()
}