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
#' @import Biobase
#' @import assertive.sets
#' @import cluster
#' @importFrom AnnotationDbi keys mget revmap
#' @importFrom reshape melt melt.array melt.data.frame melt.list
#' @importFrom tidyr spread gather separate_rows unite separate unite unnest
#' @importFrom broom tidy
#' @importFrom readr read_tsv
#' @importFrom rlang sym
#' @importFrom DEGreport degPatterns degPCA
#' @importFrom dplyr select arrange summarise rowwise mutate filter 
#'             if_else group_by "%>%" distinct n left_join inner_join
#'             bind_rows ungroup summarise_all funs
#' @importFrom tibble rownames_to_column column_to_rownames
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
#' @importFrom cowplot plot_grid
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid textGrob unit unit.pmax
#' @importFrom graphics pairs
#' @importFrom stats as.dist as.hclust cutree dist hclust
#'             dnorm pnorm predict qnorm rnorm deviance
#'             model.matrix nlminb p.adjust prop.test
#' @importFrom stringr str_extract

#' 
"_PACKAGE"