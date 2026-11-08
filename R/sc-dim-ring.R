##' @title sc_dim_geom_ring
##'
##' @description
##' Add concentric annotation rings around a dimension reduction plot.
##'
##' The rings are drawn *outside* the point cloud and outside the current
##' panel range; `ggplot2` expands the panel on its own, so the embedding is
##' never rescaled.  Each cell is assigned to a sector (by default the
##' clusters used to colour `sc_dim()`), the sectors are laid out around a
##' circle and the ring values are drawn as coloured arcs inside them.
##'
##' A ring is a *composition* of the cells in a sector, not a cell-by-cell
##' position: the arc lengths are proportional to the number of cells of each
##' level.  Consecutive rings share the same sectors, so they line up.
##'
##' @param mapping aesthetic mapping, only `colour` (or `fill`) is used and it
##' must be a column of the metadata.  Default is the sector variable, which
##' gives the ring that matches the point colours.
##' @param sector the metadata column that defines the sectors.  Default is
##' `NULL`, which uses the variable mapped to `colour` in the plot.
##' @param width how the sector width is computed: `"n"` is proportional to the
##' number of cells (the honest scale), `"sqrt"` and `"log10"` compress the
##' large sectors so that small clusters stay visible.
##' @param gap,gap_last the gap (in degrees) between sectors, and the gap
##' between the last and the first sector (defaults to `gap`).
##' @param palette colours used for the ring levels; a named vector is matched
##' by name, otherwise the values are taken in the order of the levels.  The
##' default is the qualitative palette of `ggplot2`.
##' @param r the inner and outer radius of the ring; by default the ring is
##' placed just outside the points and outside any ring already added.
##' @param thickness the ring width; by default 4% of the radius of the point
##' cloud.
##' @param legend whether to draw a legend for this ring.  `NULL` (the default)
##' hides the legend when the ring shows the same variable that already colours
##' the points (the point legend would be duplicated) and shows it otherwise.
##' `TRUE`/`FALSE` force the behaviour.  Each ring with a legend gets its own
##' fill scale, so more than one legend-bearing ring requires \pkg{ggnewscale}.
##' @param legend.title the legend title; defaults to the name of the variable
##' the ring displays.
##' @param ... additional parameters passed to
##' [geom_arc_bar][ggforce::geom_arc_bar], e.g. `linewidth`.
##' @return layer of arcs to add to a dimension reduction plot
##' @seealso [sc_dim()], [sc_dim_geom_ring_label()]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- suppressWarnings(logNormCounts(sce))
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p <- sc_dim(sce, reduction = 'UMAP')
##' p +
##'   sc_dim_geom_ring() +
##'   sc_dim_geom_ring(mapping = aes(colour = Treatment)) +
##'   sc_dim_geom_ring_label()
##' @export
sc_dim_geom_ring <- function(mapping = NULL, sector = NULL,
                             width = c("n", "sqrt", "log10"),
                             gap = 2, gap_last = gap,
                             palette = NULL, r = NULL, thickness = NULL,
                             legend = NULL, legend.title = NULL,
                             ...) {
    structure(
        list(mapping = mapping,
             sector = sector,
             width = match.arg(width),
             gap = gap,
             gap_last = gap_last,
             palette = palette,
             r = r,
             thickness = thickness,
             legend = legend,
             legend.title = legend.title,
             params = list(...)),
        class = 'sc_dim_geom_ring'
    )
}

##' @method ggplot_add sc_dim_geom_ring
##' @importFrom ggplot2 coord_fixed
##' @importFrom rlang .data
##' @export
ggplot_add.sc_dim_geom_ring <- function(object, plot, object_name, ...) {
    rlang::check_installed('ggforce', 'for `sc_dim_geom_ring()`.')

    if (!inherits(plot$facet, 'FacetNull')) {
        cli::cli_abort(c(
            "`sc_dim_geom_ring()` does not support faceted plots.",
            "x" = "The ring would be drawn in every panel.",
            "i" = "Add the rings to a single-panel dimension reduction plot."
        ))
    }

    xy <- .dim_xy_vars(plot)
    state <- plot$ggsc_ring
    if (is.null(state)) {
        state <- .ring_layout(plot, xy, object)
    }

    var <- .ring_var(object, plot, state)
    lv <- .ring_levels(plot$data[[var]])
    rr <- object$r
    if (is.null(rr)) {
        r0 <- state$r_next
        r1 <- r0 + (object$thickness %||% state$thickness)
    } else {
        if (length(rr) == 1L) rr <- c(rr, rr + (object$thickness %||% state$thickness))
        if (length(rr) != 2L || anyNA(rr) || rr[2L] <= rr[1L]) {
            cli::cli_abort("`r` must be a single radius or an increasing pair.")
        }
        r0 <- rr[1L]
        r1 <- rr[2L]
    }

    arcs <- .ring_arcs(state, plot$data, var, lv, r0, r1)

    cols <- object$palette
    if (is.null(cols)) {
        ## the rings are categorical, so a qualitative palette; with a single
        ## sector variable this also matches the colours `sc_dim()` gives the
        ## points.
        cols <- scales::hue_pal()(length(lv))
    }
    if (!is.null(names(cols)) && any(lv %in% names(cols))) {
        level_cols <- cols[lv]
    } else {
        level_cols <- rep_len(cols, length(lv))
    }
    names(level_cols) <- lv

    ## Whether this ring gets a legend.  By default a ring that repeats the
    ## variable already colouring the points stays silent (the point legend is
    ## enough); any other metadata track gets its own legend.
    show_legend <- object$legend
    if (is.null(show_legend)) {
        show_legend <- !identical(var, .dim_colour_var(plot))
    }

    params <- object$params
    params$data <- arcs
    params$inherit.aes <- FALSE
    if (is.null(params$colour)) params$colour <- NA
    params$show.legend <- show_legend

    p <- plot
    if (show_legend) {
        ## Map fill to the level and let a discrete scale produce the legend.
        ## Each legend-bearing ring needs its own fill scale; ggnewscale frees
        ## the previous one so the rings do not overwrite each other's legend.
        if (state$n_fill_scales > 0L) {
            rlang::check_installed('ggnewscale',
                'for more than one ring legend in `sc_dim_geom_ring()`.')
            p <- p + ggnewscale::new_scale_fill()
        }
        arcs$group <- factor(arcs$group, levels = lv)
        params$data <- arcs
        params$mapping <- aes(x0 = .data$x0, y0 = .data$y0, r0 = .data$r0,
                              r = .data$r1,
                              start = .data$start * pi / 180,
                              end = .data$end * pi / 180,
                              fill = .data$group)
        p <- ggplot_add(do.call(ggforce::geom_arc_bar, params), p, object_name, ...)
        p <- p + ggplot2::scale_fill_manual(
            values = level_cols,
            name = object$legend.title %||% var)
        state$n_fill_scales <- state$n_fill_scales + 1L
    } else {
        ## `I()` keeps the colours as they are, so the ring adds no fill scale
        ## and no legend, and never fights the scale `sc_dim()` uses for points.
        arcs$fill <- unname(level_cols[arcs$group])
        params$data <- arcs
        params$mapping <- aes(x0 = .data$x0, y0 = .data$y0, r0 = .data$r0,
                              r = .data$r1,
                              start = .data$start * pi / 180,
                              end = .data$end * pi / 180,
                              fill = I(.data$fill))
        p <- ggplot_add(do.call(ggforce::geom_arc_bar, params), p, object_name, ...)
    }

    if (!.has_fixed_coord(p$coordinates)) {
        p <- p + coord_fixed()
    }

    state$r_max <- r1
    state$r_next <- r1 + state$thickness * 0.35
    state$rings <- c(state$rings, var)
    p$ggsc_ring <- state
    p
}

##' @title sc_dim_geom_ring_label
##'
##' @description
##' Label the sectors of the rings added by [sc_dim_geom_ring()].  The labels
##' are placed outside the outermost ring and rotated tangentially.
##'
##' @param size,colour,fontface appearance of the labels
##' @param padding distance between the outermost ring and the labels, in
##' units of the ring thickness
##' @param ... additional parameters passed to [geom_text][ggplot2::geom_text]
##' @return layer of labels to add to a dimension reduction plot
##' @seealso [sc_dim_geom_ring()]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- suppressWarnings(logNormCounts(sce))
##' colLabels(sce) <- clusterCells(sce, assay.type = 'logcounts')
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' sc_dim(sce, reduction = 'UMAP') +
##'   sc_dim_geom_ring() +
##'   sc_dim_geom_ring_label()
##' @export
sc_dim_geom_ring_label <- function(size = 3, colour = 'black',
                                   fontface = 'bold', padding = 2, ...) {
    structure(
        list(size = size, colour = colour, fontface = fontface,
             padding = padding, params = list(...)),
        class = 'sc_dim_geom_ring_label'
    )
}

##' @method ggplot_add sc_dim_geom_ring_label
##' @importFrom ggplot2 geom_text
##' @importFrom rlang .data
##' @export
ggplot_add.sc_dim_geom_ring_label <- function(object, plot, object_name, ...) {
    state <- plot$ggsc_ring
    if (is.null(state)) {
        cli::cli_abort(c(
            "There is no ring to label.",
            "x" = "`sc_dim_geom_ring_label()` needs the layout of a ring.",
            "i" = "Add `sc_dim_geom_ring()` before the labels."
        ))
    }

    th <- state$tab$start + state$tab$span / 2
    r <- state$r_max + state$thickness * object$padding
    angle <- 90 - th
    angle <- ifelse(angle < -90 | angle > 90, angle + 180, angle)

    d <- data.frame(x = state$cx + r * cos(th * pi / 180),
                    y = state$cy + r * sin(th * pi / 180),
                    label = state$tab$sector,
                    angle = angle)

    params <- object$params
    params$data <- d
    params$mapping <- aes(x = .data$x, y = .data$y,
                          label = .data$label, angle = .data$angle)
    params$inherit.aes <- FALSE
    params$size <- object$size
    params$colour <- object$colour
    if (is.null(params$fontface)) params$fontface <- object$fontface

    ggplot_add(do.call(geom_text, params), plot, object_name, ...)
}

## `coord_fixed()` is an alias of `coord_cartesian(ratio = 1)` in ggplot2 >= 4.0,
## where the resulting object keeps the `CoordCartesian` class, so the class
## alone does not tell us whether the aspect ratio is already fixed.
.has_fixed_coord <- function(coord) {
    inherits(coord, 'CoordFixed') || isTRUE(coord$ratio == 1)
}

## `%||%` without adding a rlang import for one operator.
`%||%` <- function(x, y) if (is.null(x)) y else x

## Build the sector layout from the plot data: the centre and the radius of
## the point cloud, plus one row per sector (its cells, start angle and span).
.ring_layout <- function(plot, xy, object) {
    d <- plot$data
    sector <- object$sector
    if (is.null(sector)) {
        sector <- .dim_colour_var(plot)
    }
    ok <- length(sector) == 1L && !is.na(sector) &&
        sector %in% colnames(d) && !is.numeric(d[[sector]])
    if (!ok) {
        cli::cli_abort(c(
            "The sector column could not be determined.",
            "x" = "The plot has no non-numeric column to sector the rings by.",
            "i" = "Map a non-numeric column to `colour`, or supply `sector = ...`."
        ))
    }

    x <- d[[xy[1]]]
    y <- d[[xy[2]]]
    cx <- mean(range(x, na.rm = TRUE))
    cy <- mean(range(y, na.rm = TRUE))
    rc <- max(sqrt((x - cx)^2 + (y - cy)^2), na.rm = TRUE)

    lv <- .ring_levels(d[[sector]])
    n <- as.integer(table(factor(d[[sector]], levels = lv)))
    w <- switch(object$width,
                n = n,
                sqrt = sqrt(n),
                log10 = log10(n) + 1)
    span <- (360 - object$gap * (length(lv) - 1L) - object$gap_last) * w / sum(w)
    if (any(span <= 0)) {
        cli::cli_abort(c(
            "The sectors do not fit in 360 degrees.",
            "i" = "Use a smaller `gap`/`gap_last` or fewer sectors."
        ))
    }
    start <- c(0, cumsum(utils::head(span, -1L) + object$gap))

    list(xy = xy, cx = cx, cy = cy, rc = rc, sector = sector,
         tab = data.frame(sector = lv, n = n, start = start, span = span),
         thickness = rc * 0.04,
         r_max = rc * 1.02,
         r_next = rc * 1.02,
         rings = character(0),
         n_fill_scales = 0L)
}

## The variable mapped to `colour` in the plot.  `plot$labels` is only filled
## in by `ggplot_build()`, while the global mapping is already available, so
## read that first and only build when the colour comes from a layer.
.dim_colour_var <- function(plot) {
    if (!is.null(plot$mapping$colour)) {
        return(ggfun::get_aes_var(plot$mapping, 'colour'))
    }
    ggplot_build(plot)$plot$labels$colour
}

## The variable a ring displays: the mapped aesthetic, or the sector itself.
.ring_var <- function(object, plot, state) {
    m <- object$mapping
    if (is.null(m)) {
        return(state$sector)
    }
    key <- if (!is.null(m$colour)) {
        'colour'
    } else if (!is.null(m$fill)) {
        'fill'
    } else {
        NULL
    }
    if (is.null(key)) {
        return(state$sector)
    }
    var <- ggfun::get_aes_var(m, key)
    if (!var %in% colnames(plot$data)) {
        cli::cli_abort(c(
            "The {.var {var}} column is not in the plot data.",
            "i" = "Map a column of the metadata to `colour`."
        ))
    }
    var
}

## Levels of a ring variable, keeping the order of a factor and the order of
## appearance otherwise.
.ring_levels <- function(x) {
    if (is.factor(x)) {
        return(levels(droplevels(x)))
    }
    unique(as.character(x))
}

## Split every sector into one arc per level of `var`, proportional to the
## number of cells.  The arcs of all sectors tile the sector span.
.ring_arcs <- function(state, d, var, lv, r0, r1) {
    sector <- state$sector
    out <- lapply(seq_len(nrow(state$tab)), function(i) {
        si <- state$tab$sector[i]
        idx <- as.character(d[[sector]]) == si
        cnt <- as.integer(table(factor(d[[var]][idx], levels = lv)))
        keep <- cnt > 0L
        if (!any(keep)) {
            return(NULL)
        }
        cnt <- cnt[keep]
        edge <- state$tab$start[i] +
            c(0, cumsum(cnt)) / sum(cnt) * state$tab$span[i]
        data.frame(sector = si, group = lv[keep],
                   start = utils::head(edge, -1L), end = edge[-1L],
                   stringsAsFactors = FALSE)
    })
    arcs <- do.call(rbind, out)
    arcs <- arcs[arcs$end > arcs$start, , drop = FALSE]
    arcs$x0 <- state$cx
    arcs$y0 <- state$cy
    arcs$r0 <- r0
    arcs$r1 <- r1
    arcs
}
