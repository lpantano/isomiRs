# put header to input files
.put_header <- function(table){
    names(table)[c(1, 3, 4, 7, 8, 9, 10, 13, 14)] <- c("seq", "freq", "mir",
                                               "subs", "add", "t5", "t3",
                                               "DB", "ambiguity")
    table <- table[,c(1, 3, 4, 7, 8, 9, 10, 13, 14)]
    table[,2] <- as.numeric(table[,2])
    return(table)
}

# filter by relative abundance to reference
.filter_by_cov <- function(table, limit=0){
    freq <- mir <-  NULL
    tab.fil <- table[table$DB == "miRNA",]
    tab.fil.out <- as.data.frame(tab.fil %>% group_by(mir) %>%
                                   summarise(total=sum(freq)))
    tab.fil <- merge(tab.fil[ ,c(3, 1:2, 4:ncol(tab.fil)) ],
                     tab.fil.out,
                     by=1)
    tab.fil$score <- tab.fil$freq / tab.fil$total * 100
    tab.fil[tab.fil$score >= limit,]
}

.convert_to_new_version <- function(table){
    idx <- grepl("u-", table$t5)
    table$t5[idx] <- toupper(gsub("u-", "", table$t5[idx]))
    idx <- grepl("d-", table$t5)
    table$t5[idx] <- tolower(gsub("d-", "", table$t5[idx]))
    idx <- grepl("u-", table$t3)
    table$t3[idx] <- toupper(gsub("u-", "", table$t3[idx]))
    idx <- grepl("d-", table$t3)
    table$t3[idx] <- tolower(gsub("d-", "", table$t3[idx]))
    table$add <- toupper(gsub("u-", "", table$add))
    table
}

# Filter table reference
.filter_table <- function(table, cov=1){
    table <- .put_header(table)
    table <- .filter_by_cov(table,cov)
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
        d <- d[d[,2] > minc, ]
        if( nrow(table.merge) == 0){
            table.merge <- d
        }else{
            table.merge <- merge( table.merge, d, by=1, all=TRUE )
        }
    }
    row.names(table.merge) <- table.merge[ ,1]
    table.merge <- as.matrix(table.merge[ ,2:ncol(table.merge)])
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
        ref.val <- do.call(paste,table[,4:7])
        ref.val[grep("[ATGC]", ref.val, invert=TRUE)] <- "ref"
        ref.val[grep("[ATGC]", ref.val)] <- "iso"
        label <- paste(label, ref.val, sep=".")
    }
    if (iso5 == TRUE){
        label <- paste(label, table[,6], sep=".t5:")
    }
    if (seed == TRUE){
        seed.val <- as.character(table[,4])
        seed.val[grep("^[2-8][ATGC]", seed.val, invert=TRUE)] <- "0"
        label <- paste(label, seed.val, sep=".seed:")
    }
    if (iso3 == TRUE){
        label <- paste(label, table[,7], sep=".t3:")
    }
    if (add == TRUE){
        label <- paste(label, table[,5], sep=".ad:")
    }
    if (subs == TRUE){
        label <- paste(label, table[,4], sep=".mm:")
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
    pos <- pos[pos[,1] != "", ]
    nt.2 <- as.data.frame(t(as.data.frame((strsplit(pos$nt, "", fixed=2)))))
    names(nt.2) <- c("current","reference")
    names(nt.2$current) <- ""
    names(nt.2$reference) <- ""
    pos <- cbind(pos, nt.2)
    pos$size <- factor(pos$size, levels=1:25)
    return (pos[ ,c(3,4,2,5,6)])
}
