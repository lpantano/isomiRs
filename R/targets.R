.no_na <- function(v){
    return(v[!is.na(v)])
}

from_pairs_to_matrix <- function(df){
    if (!is.data.frame(df))
        error("Need a data.frame with 2 columns: gene, mir.")
    if ( !(names(df) %in% c("gene","mir")) )
        error("The columns should be names as gene and mir.")
    df$value <- 1
    ma = pairs %>% dplyr::select(gene, mir, value) %>%
        dplyr::distinct() %>% spread(mir, value, fill=0) %>% filter(!is.na(gene))
    row.names(ma) = ma$gene
    ma
}

.cor_matrix <- function(mirna, gene, target, min_val = -.6){
    stopifnot(is.matrix(mirna), is.matrix(gene), is.matrix(target))
    cor = cor( t(mirna), t(gene), method="kendall" )
    cor[cor > min_val] <- 0
    cor_target <- cor
    common_gene <- intersect(rownames(gene), rownames(target))
    common_mirna <- intersect(rownames(mirna), colnames(target))
    cor_target <- cor_target[common_mirna, common_gene]
    cor_target[ t(target[common_gene, common_mirna]) == 0 ] = 0 
    cat("Dimmension of cor matrix:", dim(cor_target), "\n")
    stopifnot(nrow(cor_target)>1 & ncol(cor_target)>1)
    t(cor_target) # cor matrix with values only if mirna-gene are in target
}

.scale <- function(ma){
    ma_scaled = t(apply(ma, 1, function(e){
        # (e - min(e))/(max(e) - min(e))
        scale(e)
    }))
    colnames(ma_scaled) = colnames(ma)
    ma_scaled
}

.cluster_exp <- function(ma){
    m = (1-cor(t(ma), method = "kendall"))
    m[m<0] = 0
    d = as.dist(m^2)
    c = cluster::diana(d, diss = TRUE, stand = FALSE)
    
    cutree(as.hclust(c), h = c$dc)
}

.run_enricher <- function(target, universe, org, genename, max_group=30){
    cat("GO enrichment with ", length(universe), " as universe", universe[1], "\n")
    cat("and ", length(target), " as query:", target[1], "\n")
    if (length(universe) > length(sel_genes)*3){
        ego <- enrichGO(sel_genes, org, ont = "BP", keytype = genename, universe = universe)
    }else{
        ego <- enrichGO(sel_genes, org, ont = "BP", keytype = genename)
    }
    if (is.null(ego))
        return(NULL)
    if (nrow(ego@result)==0)
        return(NULL)
    tab <- ego@result %>% filter(Count<max_group)
    ego_s <- ego
    ego_s@result <- tab
    ego_s <- clusterProfiler::simplify(ego_s)
    summary(ego_s)
}

.median_by_group <- function(e, g){
    sapply(levels(g), function(i){
        idx = which(g == i)
        mean(e[idx], na.rm=TRUE)
    })
}

.apply_median <- function(ma, group, minfc=0.5){
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
#' @param mirna_rse \link[SummarizedExperiment]{SummarizedExperiment} with miRNA
#' information. See details.
#' @param gene_rse \link[SummarizedExperiment]{SummarizedExperiment} with gene
#' information. See details.
#' @param target matrix with miRNAs (columns) and genes (rows)
#' target prediction values (1 if it is a target, 0 if not).
#' @param summarize character column name in colData(rse) to use to group
#' samples and compare betweem miRNA/gene expression.
#' @param min_cor numeric cutoff for correlation value that
#' will be use to consider a miRNA-gene pair as valid.
#' @examples 
#' 
#' pairs <- as.matrix(data.frame(row.names=c("gene1", "gene2"),
#'                             mirna1=c(0,1), mirna2=c(1,0)))
#' mirna_matrix <- as.matrix(data.frame(row.names=c("mirna1", "mirna2"),
#'                                     time0_1=c(1,1),time0_2=c(1.2,0.9),
#'                                     time1_1=c(8,8),time1_2=c(8.2,7.9)))
#' gene_matrix <- as.matrix(data.frame(row.names=c("gene1", "gene2"),
#'                                     time0_1=c(8,8),time0_2=c(8.2,7.9),
#'                                     time1_1=c(1,1),time1_2=c(1.2,0.9)))
#' mirna_col <- data.frame(row.names=c("time0_1","time0_2","time1_1","time1_2"),
#'                        group=c("t0","t0","t1","t1"))
#' gene_col <- data.frame(row.names=c("time0_1","time0_2","time1_1","time1_2"),
#'                        group=c("t0","t0","t1","t1"))
#'
#' mirna <- SummarizedExperiment(assays=SimpleList(norm=as.matrix(mirna_matrix)), 
#'                             colData= mirna_col)
#' gene <- SummarizedExperiment(assays=SimpleList(norm=as.matrix(gene_matrix)), 
#'                              colData= gene_col)
#' find_targets(mirna, gene, pairs)
find_targets <- function(mirna_rse, gene_rse, target, summarize="group", min_cor= -.6){
    mirna = assay(mirna_rse,"norm")[,order(colData(mirna_rse)[,summarize])]
    message("Number of mirnas ", nrow(mirna), " with these columns:", paste(colnames(mirna)))
    gene = assay(gene_rse, "norm")[,order(colData(gene_rse)[,summarize])]
    message("Number of genes ", nrow(gene), " with these columns:", paste(colnames(gene)))
    mirna_group = colData(mirna_rse)[order(colData(mirna_rse)[,summarize]),summarize]
    gene_group = colData(gene_rse)[order(colData(gene_rse)[,summarize]),summarize]
    message("Factors genes", levels(gene_group))
    message("Factors mirnas", levels(mirna_group))
    message("Order genes", gene_group)
    message("Order mirnas", mirna_group)
    
    mirna_norm <- .apply_median(as.matrix(mirna), mirna_group, minfc=0.5)
    gene_norm <- .apply_median(as.matrix(gene), gene_group, minfc=0.5)
    cor_target <- .cor_matrix(mirna_norm, gene_norm, as.matrix(target), min_cor)
    return(cor_target)
}

#' Clustering miRNAs-genes pairs in similar pattern expression
#' 
#' @param mirna_rse \link[SummarizedExperiment]{SummarizedExperiment} with miRNA
#' information. See details.
#' @param gene_rse \link[SummarizedExperiment]{SummarizedExperiment} with gene
#' information. See details.
#' @param target matrix with miRNAs (columns) and genes (rows)
#' target prediction (1 if it is a target, 0 if not).
#' @param summarize character column name in colData(rse) to use to group
#' samples and compare betweem miRNA/gene expression.
#' @param org \link[AnnotationDbi]{AnnotationDb}. (org.Mm.eg.db)
#' @param genename character keytype of the gene
#' names in gene_rse object.
#' @param min_cor numeric cutoff to consider a miRNA to regulate a target
#' @details 
#' 
#' This function will correlate miRNA and gene expression data using 
#' a specific metadata variable to group samples and detect pattern
#' of expression that will be annotated with GO terms.
#' mirna_rse and gene_rse can be created using the following code:
#' 
#' \code{mi_rse = SummarizedExperiment(assays=SimpleList(norm=mirna_matrix), colData, metadata=list(sign=mirna_keep))}
#' 
#' where, \code{mirna_matrix} is the normalized counts expression, 
#' \code{colData} is the metadata information and \code{mirna_keep}
#' the list of miRNAs to be used by this function.
#' @examples
#' 
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' data(isoExample)
#' ego <- enrichGO(row.names(assay(gene_ex_rse, "norm")), org.Mm.eg.db, ont = "BP", keytype="ENSEMBL")
#' data = isoNetwork(mirna_ex_rse, gene_ex_rse, ma_ex, org=ego@result)
#' isoPlotNet(data)
isoNetwork <- function(mirna_rse, gene_rse, target, 
                       summarize="group", org, genename="ENSEMBL", min_cor = -.6){
    stopifnot(class(gene_rse)=="SummarizedExperiment")
    stopifnot(class(mirna_rse)=="SummarizedExperiment")
    stopifnot(is.data.frame(target) | is.matrix(target))
    
    target = as.matrix(target)
    mirna = assay(mirna_rse,"norm")[,order(colData(mirna_rse)[,summarize])]
    message("Number of mirnas ", nrow(mirna), " with these columns:", paste(colnames(mirna)))
    gene = assay(gene_rse, "norm")[,order(colData(gene_rse)[,summarize])]
    message("Number of genes ", nrow(gene), " with these columns:", paste(colnames(gene)))
    mirna_de = as.character(metadata(mirna_rse)$sign)
    gene_de = as.character(metadata(gene_rse)$sign)
    mirna_group = colData(mirna_rse)[order(colData(mirna_rse)[,summarize]),summarize]
    gene_group = colData(gene_rse)[order(colData(gene_rse)[,summarize]),summarize]
    
    if (!all(levels(mirna_group) == levels(gene_group)))
        stop("levels in mirna and gene data are not the same")
    
    mirna_norm <- .apply_median(mirna[mirna_de,], mirna_group, minfc=0.5)
    message("Number of mirnas ", nrow(mirna_norm), " with these columns:", paste(colnames(mirna_norm)))
    gene_norm <- .apply_median(gene[gene_de,], gene_group, minfc=0.5)
    message("Number of mirnas ", nrow(gene_norm), " with these columns:", paste(colnames(gene_norm)))
    cor_target <- .cor_matrix(mirna_norm, gene_norm, target, min_cor)
    
    is_target_and_de <- rownames(cor_target)[apply(cor_target, 1, min) != 0]
    profile <- rownames(gene)[rowSums(gene>0)>3]
    stopifnot(length(is_target_and_de)>10)
    if (class(org)=="character")
        res <- .run_enricher(is_target_and_de, profile, org, genename=genename)
    if (class(org)=="data.frame")
        res <- org
    if ( is.null(res) )
        stop("No significant genes.")
    
    net <- do.call(rbind, apply(res, 1, function(x){
        .genes <- unlist(strsplit(x[8], split = "/"))
        .idx <- which(colSums(cor_target[.genes,,drop=FALSE] != 0 ) > 0)
        if ( length(.idx) > 0 & length(.genes) > 2 ){
            .tab = data.frame(reshape::melt.array(as.matrix(cor_target[.genes, .idx]))) %>%
            dplyr::filter(value != 0) %>%
            dplyr::select(gene=X1, mir=X2) %>%
            mutate(goid=x[1], go=x[2])
            return(.tab)
        }
        return(data.frame())
    }))
    
    res_by_mir <- net %>% group_by(mir, go) %>% dplyr::summarise(n()) %>%
        group_by(go) %>% dplyr::summarise(nmir=n())
    res <- res[res$Description %in% res_by_mir$go,]
    res[match(res_by_mir$go, res$Description), "nmir"] <- res_by_mir$nmir
    obj = .viz_mirna_gene_enrichment(list(network = net, summary = res),
                                     mirna_norm, gene_norm, mirna_group, org)
    list(network = net, summary = res, analysis=obj)
}

.plot_profiles = function(groups, ma){
    lapply(unique(groups), function(g){
        .ma = ma[names(groups)[groups==g],]
        ma_long = suppressMessages(reshape::melt.array(as.matrix(.ma)))
        names(ma_long) = c("gene", "group", "average")
        ma_long$group = factor(ma_long$group, gtools::mixedsort(levels(ma_long$group)))
        n = round(length(unique(ma_long$group)) / 2)
        ggplot(ma_long, aes(x=group, y=average)) + 
            # geom_boxplot(outlier.size = 0.2, size=0.5) +
            stat_smooth(data=ma_long, size=0.5,
                        aes(x=group, y=average, group=1),method = "lm",formula = y~poly(x,n)) +
            theme_bw(base_size = 11) + xlab(NULL) + ylab(NULL) +
            theme(axis.text.x=element_blank()) +
            theme(axis.text.y=element_blank(),axis.ticks=element_blank())
    })
}

.viz_mirna_gene_enrichment <- function(obj, mirna_norm, mrna_norm, group, org, plot=FALSE){
    net <- obj$network
    summary <- obj$summary
    
    ma_g = .scale(as.data.frame(mrna_norm)[as.character(unique(net$gene)),])
    ma_m = .scale(as.data.frame(mirna_norm)[as.character(unique(net$mir)),])
    groups = .cluster_exp(mrna_norm[as.character(unique(net$gene)),])
    final_df = data.frame()
    for (x in summary$Description){
        if (plot)
            cat("## GO term:", x)
        .m = as.character(unique(net[net$go==x, "mir"]))
        .g = as.character(unique(net[net$go==x, "gene"]))
        .df = reshape::melt.array(rbind(ma_g[.g,], ma_m[.m,]))
        .groups = unique(groups[.g])
        for (c in unique(.groups)){
            .in_g = intersect(names(groups[groups==c]), .g)
            .reg = intersect(net[net$gene %in% .in_g, "mir"], .m)
            .net = .df %>% dplyr::filter(X1 %in% c(.in_g, .reg)) %>% mutate(type="gene")
            .net[.net$X1 %in% .reg, "type"] = "miR"
            .net$X2 = factor(.net$X2, levels = levels(group))
            if (length(.in_g)<4)
                next
            p = ggplot(.net, aes(x=X2, y=value, colour=type, group=X1)) + 
                geom_line(size=0.1) +
                stat_smooth(aes(x=X2, y=value, group=type),se=T,method = "lm",
                            formula = y~poly(x,3)) +
                ggtitle(paste("Group:", x, "(", length(.in_g), " genes )")) +
                theme_bw(base_size = 11) +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                xlab("") + ylab("scaled normalized expression.")
            if (plot)
                print(p)
            out = data.frame(gene = paste(.in_g,
                                          collapse=","),
                             mir = paste(.reg, collapse=","),
                             ngene = length(.in_g),
                             group = c,
                             term = x)
            final_df = rbind(final_df, out)
        }
    }
    list(table=final_df, profiles=.plot_profiles(groups, ma_g))
}

.reduce_mirna = function(mirs){
    sapply(mirs, function(x){
        v=unlist(strsplit(x, split = "-"))[2:3]
        v[2] = gsub("[a-z]", "", v[2])
        paste(v, collapse = "-")
    })
}

.summary_mirna = function(df){
    do.call(rbind, apply(df, 1, function(r){
        exp = unlist(strsplit(r[2], split = ","))
        exp = .reduce_mirna(exp)
        data.frame(term=r[5], mir=exp, row.names = NULL)
    }))
}

#' Functional miRNA / gene expression profile plot
#' 
#' @param obj output from \link{isoNetwork}
isoPlotNet = function(obj){
    df = obj$analysis$table
    ma = as.matrix(df %>% dplyr::select(term, group, ngene) %>% spread(group, ngene, fill=0))
    df = df[rowSums(ma>4)>1,]
    
    df$term_short = sapply(df$term, function(x){
        paste0(substr(x, 1, 50), "...")})
    
    df$term_short = factor(df$term_short, 
                           levels=unique(df$term_short[order(df$term)]))
    
    terms_vs_profile = 
        ggplot(df, aes(x=as.factor(group), y=term_short, size=ngene)) +
        geom_point(color="grey75") + 
        geom_text(aes(label=ngene), size=5) + scale_size("", guide = FALSE) +
        theme_bw() + xlab("profiles") + ylab("") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        ggtitle("Number genes in each term and expression profile") +
        theme(plot.title = element_text(size = 10)) +
        theme(axis.text.y = element_text(size = 11))
    
    mirna_df = .summary_mirna(df)
    ma = mirna_df %>% group_by(mir) %>% dplyr::summarise(total=n())
    mirna_df = mirna_df[mirna_df$mir %in% ma$mir[ma$total>1], ]
    mirna_ma = mirna_df %>% distinct() %>% mutate(value=1) %>% spread(mir, value, fill=0)
    
    ma_temp = as.matrix(mirna_ma[,2:ncol(mirna_ma)])
    d = dist(t(ma_temp), method = "binary")
    hc = hclust(d, method="ward.D")
    mirna_df$mir = factor(mirna_df$mir, levels=hc$labels[hc$order])
    
    mirna_df$term = factor(mirna_df$term, levels=sort(unique(df$term)))
    
    terms_vs_mirnas = 
        ggplot(mirna_df, aes(x=mir, y=term)) +
        geom_point() + ggtitle("miRNAs targeting that term") +
        theme_bw() + ylab("") + xlab("") +
        theme(axis.text.x = element_text(size = 10, angle = 90, hjust = 1, vjust=0.5)) +
        theme(axis.text.y=element_blank(),axis.ticks=element_blank()) +
        theme(plot.title = element_text(size = 10))
              
    profiles = obj$analysis$profiles
    opts = theme(plot.margin = rep(unit(0,"null"),4),
                 panel.margin = unit(0,"null"),
                 axis.ticks.length = unit(0,"null")
                 )
    pp.profiles <- arrangeGrob(padding = unit(0, "line"),
                               grobs=lapply(sort(unique(df$group)), function(x){
                                   dt = data.frame(x=3, y=0.9, p=x)
                                   ggplotGrob(profiles[[x]] + opts +
                                                  geom_text(data=dt,size=4,aes(x=x,y=y,label=p)))
                                   }), ncol=4)
    

    p1 <- ggplot_gtable(ggplot_build(terms_vs_profile))
    p2 <- ggplot_gtable(ggplot_build(terms_vs_mirnas))
    maxWidth = unit.pmax(p1$heights[4:5], p2$heights[4:5])
     
    p1$heights[4:5] <- maxWidth
    p2$heights[4:5] <- maxWidth
    p = grid.arrange(top=textGrob("miRNA-gene-term interaction"),
                     arrangeGrob(
                 arrangeGrob(p1,p2, widths=c(2,3), ncol=2),
                 arrangeGrob(pp.profiles), heights = c(2,1)
                 ))
    
    #p = grid.arrange(top=textGrob("summary"),
    #                 arrangeGrob(pp.terms,pp.mirs, widths=c(2,1), ncol=2),
    #                 arrangeGrob(profiles, ncol=3),
    #                 nrow=2, heights=c(2,2))
    invisible(p)
}