##' @title sc_dim_fraction
##'
##' @description
##' Cell-type composition as a bar plot, to place *beside* a dimension
##' reduction plot (e.g. with [aplot::plot_list()]) instead of around it.
##'
##' This is the companion view of [sc_dim_geom_ring()]: the rings show the same
##' metadata composition wrapped around the embedding, while `sc_dim_fraction()`
##' shows it as an ordinary bar plot that can be assembled with \pkg{aplot}.
##' The fill colours are taken from the colour scale of `sc_dim_plot`, so a cell
##' type has the same colour in the embedding and in the bar plot.
##'
##' @param sc_dim_plot a dimension reduction plot produced by [sc_dim()].
##' @param group.by the metadata column defining the bars (one bar per level),
##' e.g. the sample, treatment or condition.
##' @param celltype the metadata column defining the cell types to stack.
##' Default is `NULL`, which uses the variable mapped to `colour` in the plot.
##' @param position `"fill"` (the default) normalises each bar to proportions
##' summing to 1; `"stack"` shows raw cell counts.
##' @param ... additional parameters passed to [geom_col][ggplot2::geom_col].
##' @return a `ggplot` bar plot of the cell-type composition
##' @importFrom ggplot2 ggplot aes geom_col scale_y_continuous labs theme_minimal
##'   theme element_blank scale_fill_manual waiver
##' @seealso [sc_dim()], [sc_dim_geom_ring()], [sc_dim_count()]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- suppressWarnings(logNormCounts(sce))
##' colLabels(sce) <- clusterCells(sce, assay.type = 'logcounts')
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p <- sc_dim(sce, reduction = 'UMAP')
##' sc_dim_fraction(p, group.by = 'Treatment')
##'
##' ## assemble beside the embedding with aplot.  The bar plot has its own
##' ## discrete x axis, so it is stacked as a separate panel (`plot_list`)
##' ## rather than axis-aligned (`insert_bottom`, which would drop the bars).
##' if (requireNamespace('aplot', quietly = TRUE)) {
##'   pf <- sc_dim_fraction(p, group.by = 'Treatment') +
##'     theme(legend.position = 'none')   # the embedding already has the legend
##'   aplot::plot_list(gglist = list(p, pf), ncol = 1, heights = c(3, 1.2))
##' }
##' @export
sc_dim_fraction <- function(sc_dim_plot, group.by, celltype = NULL,
                            position = c("fill", "stack"), ...) {
    position <- match.arg(position)
    d <- sc_dim_plot$data

    ct <- celltype %||% .dim_colour_var(sc_dim_plot)
    .fraction_check_col(d, ct, "celltype")
    .fraction_check_col(d, group.by, "group.by")
    if (is.numeric(d[[ct]])) {
        cli::cli_abort(c(
            "`celltype` must be a categorical column.",
            "x" = "The {.var {ct}} column is numeric."
        ))
    }

    tab <- as.data.frame(table(d[[ct]], d[[group.by]]))
    names(tab) <- c("celltype", "group", "n")
    tab <- tab[tab$n > 0, , drop = FALSE]

    ## Keep the embedding's cell-type order and colours so the two views agree.
    cmap <- if (identical(ct, .dim_colour_var(sc_dim_plot))) {
        .dim_colour_map(sc_dim_plot)
    } else {
        NULL
    }
    ct_levels <- if (!is.null(cmap)) names(cmap) else .ring_levels(d[[ct]])
    tab$celltype <- factor(tab$celltype, levels = ct_levels)
    tab$group <- .as_group_factor(tab$group)

    p <- ggplot(tab, aes(x = .data$group, y = .data$n, fill = .data$celltype)) +
        geom_col(position = position, ...) +
        scale_y_continuous(
            expand = c(0, 0),
            labels = if (position == "fill") scales::label_percent() else waiver()) +
        labs(x = NULL,
             y = if (position == "fill") "Fraction of cells" else "Cell count",
             fill = ct) +
        theme_minimal() +
        theme(panel.grid.major.x = element_blank())

    if (!is.null(cmap)) {
        p <- p + scale_fill_manual(values = cmap, name = ct)
    }
    p
}

## Validate that a fraction column exists and is usable.
.fraction_check_col <- function(d, col, arg) {
    if (length(col) != 1L || is.na(col) || !col %in% colnames(d)) {
        cli::cli_abort(c(
            "The {.var {arg}} column could not be found in the plot data.",
            "i" = "It must be a column of the metadata stored in the {.cls sc_dim} plot."
        ))
    }
    invisible(TRUE)
}

## A discrete x axis: keep factor levels, otherwise use the order of appearance.
.as_group_factor <- function(x) {
    if (is.factor(x)) {
        return(x)
    }
    factor(x, levels = unique(x))
}

## The trained cell-type -> colour mapping of the plot's colour scale, so a
## companion plot can reuse exactly the same colours.  Returns NULL when the
## plot has no discrete colour scale.
.dim_colour_map <- function(plot) {
    var <- .dim_colour_var(plot)
    if (is.null(var)) {
        return(NULL)
    }
    sc <- ggplot_build(plot)$plot$scales$get_scales("colour")
    if (is.null(sc) || !sc$is_discrete()) {
        return(NULL)
    }
    lv <- sc$get_limits()
    stats::setNames(sc$map(lv), lv)
}
