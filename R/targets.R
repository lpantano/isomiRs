.no_na <- function(v){
    return(v[!is.na(v)])
}

pairsMatrix <- function(df){
    if (!is.data.frame(df))
        stop("Need a data.frame with 2 columns: gene, mir.")
    if (sum(!(names(df) %in% c("gene","mir"))) > 0){
        message("Assuming columns names as gene and mir.")
        names(df) <- c("gene", "mir")
    }
    df[["value"]] <- 1
    ma = df %>% .[c("gene", "mir", "value")] %>%
        dplyr::distinct() %>% spread(!!sym("mir"),
                                     !!sym("value"),
                                     fill = 0) %>%
        filter(!is.na(gene))
    row.names(ma) = ma[["gene"]]
    ma
}

.cor_matrix <- function(mirna, gene, target, min_val = -.6){
    stopifnot(is.matrix(mirna), is.matrix(gene), is.matrix(target))
    cor = cor(t(mirna), t(gene), method = "kendall")
    cor[cor > min_val] <- 0
    cor_target <- cor
    common_gene <- intersect(rownames(gene), rownames(target))
    common_mirna <- intersect(rownames(mirna), colnames(target))
    cor_target <- cor_target[common_mirna, common_gene]
    cor_target[ t(target[common_gene, common_mirna]) == 0 ] = 0
    cat("Dimmension of cor matrix:", dim(cor_target), "\n")
    stopifnot(nrow(cor_target) > 1 & ncol(cor_target) > 1)
    t(cor_target) # cor matrix with values only if mirna-gene are in target
}

.scale <- function(ma){
    ma_scaled = t(apply(ma, 1, function(e){
        scale(e)
    }))
    colnames(ma_scaled) = colnames(ma)
    ma_scaled
}

.cluster_exp <- function(ma){
    meta <- data.frame(row.names = colnames(ma),
                      xaxis = colnames(ma))
    r <- degPatterns(ma,
                     meta, time = "xaxis", minc = 0, plot = FALSE)
    groups <- r$df$cluster
    names(groups) <- r$df$genes
    groups
}


.median_by_group <- function(e, g){
    sapply(levels(g), function(i){
        idx = which(g == i)
        median(e[idx], na.rm = TRUE)
    })
}

.apply_median <- function(ma, group, minfc = 0.5){
    new_exp = t(apply(ma, 1, function(x, g){
        .median_by_group(x, g)
    }, group))
    new_exp[abs(rowMax(new_exp) - rowMin(new_exp)) > minfc, ]
}

#' Find miRNAs target using mRNA/miRNA expression
#'
#' This function creates a matrix with rows (genes) and
#' columns (mirnas) with values indicating if miRNA-gene
#' pair is target according putative targets and negative
#' correlation of the expression of both molecules.
#'
#' @param mirna_rse \code{SummarizedExperiment} with miRNA
#' information. See details.
#' @param gene_rse \code{SummarizedExperiment} with gene
#' information. See details.
#' @param target Data.frame with two columns: gene and miRNA.
#' @param summarize Character column name in colData(rse) to use to group
#'   samples and compare betweem miRNA/gene expression.
#' @param min_cor Numeric cutoff for correlation value that
#'   will be use to consider a miRNA-gene pair as valid.
#' @examples
#'
#' data(isoExample)
#' mirna_ma <- data.frame(gene = names(gene_ex_rse)[1:20],
#'                        mir = names(mirna_ex_rse))
#' corMat <- findTargets(mirna_ex_rse, gene_ex_rse, mirna_ma)
#' @return mirna-gene matrix
#' @export
findTargets <- function(mirna_rse, gene_rse, target,
                         summarize = "group", min_cor= -.6){
    target = pairsMatrix(target)
    mirna = assay(mirna_rse,"norm")[, order(colData(mirna_rse)[, summarize])]
    mirna = mirna[intersect(colnames(target), rownames(mirna)),]
    message("Number of mirnas ", nrow(mirna))
    gene = assay(gene_rse, "norm")[, order(colData(gene_rse)[, summarize])]
    gene = gene[intersect(rownames(target), rownames(gene)),]
    message("Number of genes ", nrow(gene))
    mirna_group = droplevels(colData(mirna_rse)[order(colData(mirna_rse)[, summarize]), summarize])
    gene_group = droplevels(colData(gene_rse)[order(colData(gene_rse)[, summarize]), summarize])
    lvs_common_group = intersect(levels(mirna_group), levels(gene_group))
    mirna_group = droplevels(mirna_group[mirna_group  %in%  lvs_common_group])
    gene_group = droplevels(gene_group[gene_group  %in%  lvs_common_group])
    mirna = mirna[,mirna_group  %in%  lvs_common_group]
    gene = gene[,gene_group  %in%  lvs_common_group]
    stopifnot(is.factor(mirna_group))
    stopifnot(is.factor(gene_group))
    message("Factors genes", paste(levels(gene_group)))
    message("Factors mirnas", paste(levels(mirna_group)))
    message("Order genes", paste(gene_group))
    message("Order mirnas", paste(mirna_group))

    mirna_norm <- .apply_median(as.matrix(mirna), mirna_group, minfc = 0.5)
    gene_norm <- .apply_median(as.matrix(gene), gene_group, minfc = 0.5)
    
    message("Calculating correlation matrix")
    cor_target <- .cor_matrix(mirna_norm, gene_norm, as.matrix(target), min_cor)
    return(cor_target)
}

.predict_correlate_mirna_targets <- function(mirna_rse,
                                             gene_rse, summarize,
                                             org, gene_id){
    message("Predict miRNA targets with mirna2targetscan")
    stopifnot(is(org, "OrgDb"))
    mirna_de = as.character(metadata(mirna_rse)$sign)
    gene_de = as.character(metadata(gene_rse)$sign)
    
    targets = mirna2targetscan(mirna_de, org, gene_id)
    imp_targets = targets[,c(gene_id, "mir")] %>% 
        filter(!!sym(gene_id)  %in% gene_de)
    findTargets(mirna_rse, gene_rse, imp_targets,
                summarize = summarize)
}

.detect_gene_symbol <- function(names){
    names <- as.character(names)
    if (grepl("^ENS", names[1]))
        return("ENSEMBL")
    if (grepl("^[0-9]*$", names[1]))
        return("ENTREZID")
    return("SYMBOL")
}

.is_mapping_needed <- function(names1, names2){
    names1 <- as.character(names1)
    names2 <- as.character(names2)
    id1 <- .detect_gene_symbol(names1)
    id2 <- .detect_gene_symbol(names2)
    if (id1 != id2)
        return(list(from = id2, to = id1))
    return(FALSE)
}

.convert_names <- function(org, names, keytype, columns){
    AnnotationDbi::select(org, keys = names,
                          keytype = keytype, columns = columns)
}
#' Clustering miRNAs-genes pairs in similar pattern expression
#'
#' Clustering miRNAs-genes pairs
#'
#' @param mirna_rse \code{SummarizedExperiment} with miRNA
#'   information. See details.
#' @param gene_rse \code{SummarizedExperiment} with gene
#'   information. See details.
#' @param target Matrix with miRNAs (columns) and genes (rows)
#'   target prediction (1 if it is a target, 0 if not).
#' @param org \code{AnnotationDb} obejct. For example:(org.Mm.eg.db)
#' @param enrich The output of clusterProfiler of similar functions.
#' @param summarize Character column name in `colData(rse)` to use to group
#'   samples and compare betweem miRNA/gene expression.
#' @param genename Character keytype of the gene
#'   names in gene_rse object.
#' @param min_cor Numeric cutoff to consider a miRNA to regulate a target.
#' @param min_fc Numeric cutoff to consider as the minimum log2FoldChange
#'   between groups to be considered in the analysis.
#' @details
#'
#' This function will correlate miRNA and gene expression data using
#' a specific metadata variable to group samples and detect pattern
#' of expression that will be annotated with GO terms.
#' mirna_rse and gene_rse can be created using the following code:
#'
#' `mi_rse = SummarizedExperiment(assays=SimpleList(norm=mirna_matrix),
#'                                colData, metadata=list(sign=mirna_keep))`
#'
#' where, `mirna_matrix` is the normalized counts expression,
#' `colData` is the metadata information and `mirna_keep`
#' the list of miRNAs to be used by this function.
#' 
#' @examples
#' # library(org.Mm.eg.db)
#' # library(clusterProfiler)
#' data(isoExample)
#' # ego <- enrichGO(row.names(assay(gene_ex_rse, "norm")),
#' #                 org.Mm.eg.db, "ENSEMBL", ont = "BP")
#' data <- isoNetwork(mirna_ex_rse, gene_ex_rse, 
#'                    summarize = "group", target = ma_ex,
#'                    enrich = ego)
#' isoPlotNet(data, minGenes = 5)
#' @return list with network information
#' @export
isoNetwork <- function(mirna_rse, gene_rse,
                       summarize = NULL,
                       target = NULL, org = NULL,
                       enrich = NULL,
                       genename = "ENSEMBL",
                       min_cor = -.6, min_fc = 0.5){
    stopifnot(summarize  %in% names(colData(mirna_rse)))
    stopifnot(is(gene_rse, "SummarizedExperiment"))
    stopifnot(is(mirna_rse, "SummarizedExperiment"))
    stopifnot("sign"  %in% names(metadata(gene_rse)))
    stopifnot("sign"  %in% names(metadata(mirna_rse)))
    mirna_de = as.character(metadata(mirna_rse)$sign)
    gene_de = as.character(metadata(gene_rse)$sign)
    if (is.null(target))
        target <- .predict_correlate_mirna_targets(mirna_rse,
                                                   gene_rse,
                                                   summarize,
                                                   org, genename)
    target = as.matrix(target)
    stopifnot(is.data.frame(target) | is.matrix(target))
    
    mirna = assay(mirna_rse,"norm")[,order(colData(mirna_rse)[,summarize])]
    message("Number of mirnas ", nrow(mirna), " with these columns:", paste(colnames(mirna)))
    gene = assay(gene_rse, "norm")[,order(colData(gene_rse)[,summarize])]
    message("Number of genes ", nrow(gene), " with these columns:", paste(colnames(gene)))
    mirna_group = colData(mirna_rse)[order(colData(mirna_rse)[,summarize]),summarize]
    mirna_group = droplevels(mirna_group)
    gene_group = colData(gene_rse)[order(colData(gene_rse)[,summarize]),summarize]
    gene_group = droplevels(gene_group)

    if (!all(levels(mirna_group) == levels(gene_group))){
        message("levels in mirna and gene data are not the same.")
        message("Reducing data to common levels.")
        common = intersect(as.character(gene_group), 
                           as.character(mirna_group))
        mirna_group = droplevels(mirna_group[mirna_group  %in% common])
        gene_group = droplevels(gene_group[gene_group  %in% common])
    }
    
    mirna_norm <- .apply_median(mirna[mirna_de,], mirna_group, min_fc)
    message("Number of mirnas ", nrow(mirna_norm),
            " with these columns:", paste(colnames(mirna_norm)))
    gene_norm <- .apply_median(gene[gene_de,], gene_group, min_fc)
    message("Number of mirnas ", nrow(gene_norm),
            " with these columns:", paste(colnames(gene_norm)))
    
    cor_target <- .cor_matrix(mirna_norm, gene_norm, target, min_cor)

    is_target_and_de <- rownames(cor_target)[apply(cor_target, 1, min) != 0]
    profile <- rownames(gene)[rowSums(gene > 0) > 3]
    stopifnot(length(is_target_and_de) > 10)
    
    res <- NULL
    if (is.null(enrich))
        stop("Run enrich method, please. See clusterProfiler or ReactomePA.")
    if (is(enrich, "enrichResult")){
        res <- slot(enrich, "result")
    }else if (is(enrich, "data.frame")){
        res <- enrich
    }else{
        stop("No significant genes.")
    }

    cor_long <- cor_target %>%
        as.data.frame() %>%
        rownames_to_column("gene") %>% 
        gather("mir", "value", -gene) %>%
        filter(value != 0)
    net <- res %>% separate_rows("geneID") %>% 
        .[,c("ID", "Description", "geneID")] 
    # browser()
    if (!are_intersecting_sets(net[["geneID"]], names(gene_rse))){
        message("No matching genes between enrich and gene_rse.")
        mapping <- .is_mapping_needed(names(gene_rse), net$geneID)
        message(" Converting from ", mapping[[1]], " to ", mapping[[2]])
        df_map <- .convert_names(org, net[["geneID"]],
                                 mapping[["from"]],
                                 mapping[["to"]])
        net <- left_join(net, df_map,
                         by = c("geneID" = mapping[["from"]]))
        net[["geneID"]] <- net[[mapping[["to"]]]]
    }
    net <- net %>%
        inner_join(cor_long, by = c("geneID" = "gene")) %>% 
        mutate(go = Description, gene = geneID)
    

    
    res_by_mir <- net %>% group_by(!!sym("mir"), !!sym("ID")) %>% 
        dplyr::summarise(n()) %>%
        group_by(!!sym("ID")) %>%
        dplyr::summarise(nmir = n())
    res <- res[res$ID %in% res_by_mir$ID,]
    res[match(res_by_mir$ID, res$ID), "nmir"] <- res_by_mir$nmir
    
    obj = .viz_mirna_gene_enrichment(list(network = net, 
                                          summary = res),
                                     mirna_norm, gene_norm, 
                                     mirna_group, org, plot = FALSE)
    list(network = net, summary = res, analysis = obj)
}

.plot_profiles = function(groups, ma){
    lapply(unique(groups), function(g){
        .ma = ma[names(groups)[groups == g],]
        ma_long = suppressMessages(reshape::melt.array(as.matrix(.ma)))
        names(ma_long) = c("gene", "group", "average")
        ma_long$group = factor(ma_long$group, 
                               gtools::mixedsort(levels(ma_long$group)))
        n = round(length(unique(ma_long$group)) / 2L)
        ggplot(ma_long, aes_string(x = "group", y = "average")) +
            stat_smooth(data = ma_long, size = 1,
                        aes_string(x = "group", y = "average",
                                   group = 1L),
                        color = "black",
                        method = "lm",formula = y~poly(x,n)) +
            theme_bw(base_size = 11L) + xlab(NULL) + ylab(NULL) +
            theme(axis.text.x = element_blank()) +
            theme(axis.text.y = element_blank(),
                  axis.ticks = element_blank())
    })
}

.viz_mirna_gene_enrichment <- function(obj, mirna_norm, mrna_norm, group, org, plot=FALSE){
    net <- as.data.frame(obj$network)
    summary <- obj$summary
    ma_g = .scale(as.data.frame(mrna_norm)[as.character(unique(net$gene)),])
    ma_m = .scale(as.data.frame(mirna_norm)[as.character(unique(net$mir)),])
    groups = .cluster_exp(mrna_norm[as.character(unique(net$gene)),])
    final_df = data.frame()
    for (x in summary$Description) {
        message(x)
        if (plot)
            cat("## GO term:", x)
        .m = as.character(unique(unlist(net[net$go == x, "mir"])))
        .g = as.character(unique(unlist(net[net$go == x, "gene"])))
        .df = reshape::melt.array(rbind(ma_g[.g,, drop = FALSE],
                                        ma_m[.m,, drop = FALSE]))
        .groups = unique(groups[.g])
        for (c in unique(.groups)) {
            .in_g = intersect(names(groups[groups == c]), .g)
            if (length(.in_g) < 2)
                next
            .reg = intersect(unlist(net[net$gene %in% .in_g, "mir"]), .m)
            .net = .df %>% dplyr::filter(X1 %in% c(.in_g, .reg)) %>%
                mutate(type = "gene")
            .net[.net$X1 %in% .reg, "type"] = "miR"
            .net$X2 = factor(.net$X2, levels = levels(group))
            p = ggplot(.net, aes_string(x = "X2", y = "value",
                                        colour = "type", group = "X1")) +
                geom_line(size = 0.1) +
                stat_smooth(aes_string(x = "X2", y = "value",
                                       group = "type"), se = TRUE,
                            method = "lm",
                            formula = y~poly(x,3)) +
                ggtitle(paste("Group:", x, "(", length(.in_g), " genes )")) +
                theme_bw(base_size = 11) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                xlab("") + ylab("scaled normalized expression.")
            if (plot)
                print(p)
            out = data.frame(gene = paste(.in_g,
                                          collapse = ","),
                             mir = paste(.reg, collapse = ","),
                             ngene = length(.in_g),
                             group = c,
                             term = x)
            final_df = rbind(final_df, out)
        }
    }
    list(table = final_df, profiles = .plot_profiles(groups, ma_g))
}

.reduce_mirna = function(mirs){
    sapply(mirs, function(x){
        v = unlist(strsplit(x, split = "-"))[2:3]
        v[2] = gsub("[a-z]", "", v[2])
        paste(v, collapse = "-")
    })
}

.summary_mirna = function(df){
    do.call(rbind, apply(df, 1, function(r){
        exp = unlist(strsplit(r[2], split = ","))
        exp = .reduce_mirna(exp)
        data.frame(term = r[5], group = r[4],
                   mir = exp, row.names = NULL)
    }))
}

.fix_term_name <- function(term){
    if (nchar(term) > 50)
        term <- paste0(substr(term, 1, 50), "...")
    words <- strsplit(term, " ")[[1]]
    splits <- ceiling(cumsum(nchar(words)) / 20) - 1
    temp = splits[1]
    short = words[1]
    if (length(words) == 1)
        return(term)
    for (s in 2:length(words)){
        if (splits[s] > temp)
            short = paste(short, words[s], sep = "\n")
        else short = paste(short, words[s], sep = " ")
        temp = splits[s]
    }
    short
}

#' Functional miRNA / gene expression profile plot
#' 
#' Plot analysis from [isoNetwork()]. See that function
#' for an example of the figure.
#' 
#' @param obj Output from [isoNetwork()].
#' @param minGenes Minimum number of genes per term to be kept.
#' @return Network ggplot.
#' @export
isoPlotNet = function(obj, minGenes = 2){
    # browser()
    df = obj$analysis$table
    df = df[df[["ngene"]] >= minGenes,]
    ma = as.matrix(df %>% .[,c("term", "group", "ngene")] %>% 
                       spread(group, ngene, fill=0))
    

    df$term_short = sapply(as.character(df$term), .fix_term_name)
    df$group = factor(df$group)
    df$term_short = factor(df$term_short,
                           levels = unique(df$term_short[order(df$term)]))

    terms_vs_profile =
        ggplot(df, aes_string(x = "group",
                              y = "term_short",
                              size = "ngene")) +
        geom_point(color = "grey75") +
        geom_text(aes_string(label = "ngene"), size = 3) +
        theme_bw() + xlab("Gene expression profiles") + ylab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle("Number genes in each term\n over expression profile") +
        theme(plot.title = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 7, angle = 0)) +
        theme(legend.position="none")

    mirna_df = .summary_mirna(df)
    ma = mirna_df %>% group_by(mir) %>% dplyr::summarise(total=n())
    mirna_df = mirna_df[mirna_df$mir %in% ma$mir[ma$total > 1], ]
    mirna_ma = mirna_df %>% dplyr::distinct() %>%
        mutate(value = 1) %>% spread(mir, value, fill=0)

    ma_temp = as.matrix(mirna_ma[,2:ncol(mirna_ma)])
    d = dist(t(ma_temp), method = "binary")
    hc = hclust(d, method = "ward.D")
    
    mirna_df$mir = factor(mirna_df$mir, levels = hc$labels[hc$order])
    mirna_df$term = factor(mirna_df$term, levels = sort(unique(df$term)))
    # browser()
    terms_vs_mirnas =
        ggplot(mirna_df, aes_string(x = "mir", y = "term")) +
        geom_text(aes_string(label = "group"), size = 3) +
        ggtitle("miRNAs targeting that term") +
        theme_bw() + ylab("") + xlab("") +
        theme(axis.text.x = element_text(size = 7, angle = 90,
                                         hjust = 1, vjust = 0.5)) +
        theme(axis.text.y = element_blank(), axis.ticks = element_blank()) +
        theme(plot.title = element_text(size = 10))

    profiles = obj$analysis$profiles
    opts = theme(plot.margin = rep(unit(0,"null"),4),
                 panel.spacing = unit(0,"null"),
                 axis.ticks.length = unit(0,"null")
                 )
    pp.profiles <- arrangeGrob(padding = unit(0, "line"),
                               grobs = lapply(sort(unique(df$group)),
                                              function(x){
                                   dt = data.frame(x = 3, y = 0.9, p = x)
                                   ggplotGrob(profiles[[x]] + opts +
                                                  geom_text(data = dt,
                                                            size = 4,
                                                            aes_string(x = "x",
                                                                       y = "y",
                                                                       label = "p")))
                                   }), ncol = 4)


    p1 <- ggplot_gtable(ggplot_build(terms_vs_profile))
    p2 <- ggplot_gtable(ggplot_build(terms_vs_mirnas))
    maxWidth = unit.pmax(p1$heights[c(2, 7)], p2$heights[c(2, 7)])
    p1$heights[c(2, 7)] <- maxWidth
    p2$heights[c(2, 7)] <- maxWidth
    p <- plot_grid(plot_grid(p1,p2, align = "h"),
              pp.profiles, rel_widths = c(1,2),
              nrow=2, rel_heights = c(3,1))
    show(p)
    invisible(p)
}


#' Find targets in targetscan database
#' 
#' From a list of miRNA names, find their targets
#' in targetscan.Hs.eg.db annotation package.
#' 
#' @param mirna Character vector with miRNA names as in
#'   miRBase 21.
#' @param species hsa or mmu supported right now.
#' @param org \code{AnnotationDb} obejct. For example:(org.Mm.eg.db)
#' @param keytype Character mentioning the gene id to use.
#'   For example, `ENSEMBL`.
#'   
#' @return [data.frame] with 4 columns:
#'   * miRFamily
#'   * Seedmatch
#'   * PCT
#'   * entrezGene
#' @examples 
#' library(targetscan.Hs.eg.db)
#' mirna2targetscan(c("hsa-miR-34c-5p"))
#' @export
mirna2targetscan <- function(mirna, species = "hsa", org = NULL, keytype = NULL){
    mirna_map <- data.frame(mir = mirna,
                            ext = gsub("-[53]p$", "", mirna),
                            num = gsub("-[1-9]$", "",
                                       gsub("-[53]p$", "",mirna)),
                            letter = gsub("[a-z]$", "",
                                          gsub("-[1-9]$", "",
                                               gsub("-[53]p$", "",mirna))),
                            stringsAsFactors = FALSE)
    if (species == "hsa"){
        egMIRNA <- keys(targetscan.Hs.egMIRNA)
        egMIRBASE2FAMILY <- targetscan.Hs.egMIRBASE2FAMILY
        egTARGETS <- revmap(targetscan.Hs.egTARGETS)
        egTARGETSFULL <- targetscan.Hs.egTARGETSFULL
    }else if (species == "mmu"){
        egMIRNA <- keys(targetscan.Mm.egMIRNA)
        egMIRBASE2FAMILY <- targetscan.Mm.egMIRBASE2FAMILY
        egTARGETS <- revmap(targetscan.Mm.egTARGETS)
        egTARGETSFULL <- targetscan.Mm.egTARGETSFULL
    }else{
        stop("Species not supported: ", species)
    }

    mirna_map[["targetscan"]]  <- "None"
    mirna_map[["targetscan"]][mirna_map[["mir"]] %in% egMIRNA] <- mirna_map[["mir"]][mirna_map[["mir"]]  %in% egMIRNA]
    mirna_map[["targetscan"]][mirna_map[["ext"]]  %in% egMIRNA] <- mirna_map[["ext"]][mirna_map[["ext"]]  %in% egMIRNA]
    mirna_map[["targetscan"]][mirna_map[["num"]]  %in% egMIRNA] <- mirna_map[["num"]][mirna_map[["num"]]  %in% egMIRNA]
    mirna_map[["targetscan"]][mirna_map[["letter"]]  %in% egMIRNA] <- mirna_map[["letter"]][mirna_map[["letter"]]  %in% egMIRNA]

    mir <- unique(mirna_map[["targetscan"]])
    
    if (sum(mirna_map[["targetscan"]] == "None") > 0)
        message("Missing miRNAs: ",
                mirna_map[["mir"]][mirna_map[["targetscan"]] == "None"])
    
    fam <- mget(intersect(mir, egMIRNA),
                egMIRBASE2FAMILY)
    genes = mget(unlist(fam), egTARGETS)
    full = mget(unlist(genes)[!is.na(unlist(genes))], egTARGETSFULL) 
    df <- lapply(names(full), function(x){
        data.frame(miRFamily = full[[x]]@miRFamily, 
                   Seedmatch = full[[x]]@Seedmatch,
                   PCT = full[[x]]@PCT,
                   gene = x,
                   stringsAsFactors = FALSE)
    }) %>% bind_rows() %>%
        .[.[["miRFamily"]] %in% unlist(fam),] %>% 
        dplyr::distinct()
    if (!is.null(org)){
        map <- AnnotationDbi::select(org, unique(df[["gene"]]),
                                     columns = keytype,
                                     keytype = "ENTREZID")
        df <- df %>% left_join(map, by = c("gene" = "ENTREZID")) 
    }
    
    df %>% 
        left_join(data.frame(targetscan = names(fam),
                             miRFamily = unlist(fam),
                             stringsAsFactors = FALSE) %>% 
                      left_join(mirna_map, by = "targetscan") %>% 
                      filter(!is.na(!!sym("mir"))),
                  by = "miRFamily")
}
