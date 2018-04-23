#' Deprecated functions in package isomiRs
#' 
#' 
#' @rdname isomiRs-deprecated
#' 
#' @details
#' The following functions are deprecated and will be made defunct; 
#' use the replacement indicated below:
#' * isoCorrection dLQNO qLQNO pLQNO rLQNOLQNO isoLQNO

dLQNO <- function()
{
    .Deprecated("scounts::dLQNO")
}

qLQNO <- function()
{
    .Deprecated("scounts::qLQNO")
}

pLQNO <- function()
{
    .Deprecated("scounts::pLQNO")
}

rLQNO <- function()
{
    .Deprecated("scounts::rLQNO")
}

LQNO <- function()
{
    .Deprecated("scounts::LQNO")
}

isoLQN <- function()
{
    .Deprecated("scounts::isoLQN")
}

isoCorrection <- function()
{
    .Deprecated("scounts::isoCorrection")
}

#' isomiRs
#'
#' @import BiocGenerics
#' @import S4Vectors
#' @import IRanges
#' @import GenomicRanges
#' @import SummarizedExperiment
#' @import methods
#' @import DESeq2
#' @importFrom limma voom
#' @import GGally
#' @import gtools
#' @import grDevices
#' @import Biobase
#' @import assertive.sets
#' @import cluster
#' @importFrom AnnotationDbi keys mget revmap
#' @importFrom reshape melt melt.array melt.data.frame melt.list
#' @importFrom tidyr spread gather separate_rows
#' @importFrom readr read_tsv
#' @importFrom rlang sym
#' @importFrom dplyr select arrange summarise rowwise mutate filter 
#'             if_else group_by "%>%" distinct n left_join inner_join
#'             bind_rows ungroup
#' @importFrom tibble rownames_to_column
#' @importFrom DiscriMiner plsDA
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @importFrom ggplot2 aes aes_string element_rect geom_jitter 
#'             ggplot element_text
#'             labs ggtitle xlab ylab scale_size facet_wrap
#'             scale_color_brewer scale_colour_brewer theme theme_bw
#'             stat_smooth coord_polar element_blank
#'             ggplot_gtable ggplot_build
#'             geom_text geom_line geom_point ggplotGrob geom_polygon
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob unit unit.pmax
#' @importFrom graphics pairs
#' @importFrom stats as.dist as.hclust cutree dist hclust
#'             dnorm pnorm predict qnorm rnorm deviance
#'             model.matrix nlminb p.adjust

#' 
"_PACKAGE"