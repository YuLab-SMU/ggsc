## Coverage for the circular ring layers.  Like the other `sc_dim_geom_*`
## tests, a hand-built plot keeps the tests fast and free of scater/scuttle.

ring_plot <- function(n = 20) {
    set.seed(1)
    d <- data.frame(
        .BarcodeID = sprintf("c%03d", seq_len(3 * n)),
        UMAP1 = rnorm(3 * n),
        UMAP2 = rnorm(3 * n),
        cluster = rep(c("a", "b", "c"), each = n),
        phase = rep(c("G1", "S"), length.out = 3 * n),
        stringsAsFactors = FALSE
    )
    ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = cluster)) +
        ggplot2::geom_point()
}

test_that(".ring_levels() keeps the factor order and the order of appearance", {
    expect_identical(.ring_levels(factor(c("b", "a"), levels = c("b", "a", "c"))),
                     c("b", "a"))
    expect_identical(.ring_levels(c("b", "a", "b")), c("b", "a"))
})

test_that("the ring layout fills 360 degrees except for the gaps", {
    p <- ring_plot() + sc_dim_geom_ring()
    tab <- p$ggsc_ring$tab

    expect_equal(nrow(tab), 3L)
    expect_equal(tab$n, rep(20L, 3))
    expect_equal(sum(tab$span) + 2 * 2 + 2, 360)
    ## sectors are laid out in the order of the levels, without overlap
    expect_equal(tab$start, c(0, 2 + tab$span[1], 4 + sum(tab$span[1:2])))
})

test_that("`width` controls how much room each sector gets", {
    d <- ring_plot()$data
    d$cluster <- c(rep("a", 40), rep("b", 15), rep("c", 5))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = cluster)) +
        ggplot2::geom_point()

    n <- (p + sc_dim_geom_ring(width = "n"))$ggsc_ring$tab$span
    expect_equal(n[1] / n[3], 8)

    s <- (p + sc_dim_geom_ring(width = "sqrt"))$ggsc_ring$tab$span
    expect_equal(s[1] / s[3], sqrt(8))

    l <- (p + sc_dim_geom_ring(width = "log10"))$ggsc_ring$tab$span
    expect_true(l[1] / l[3] < sqrt(8))
})

test_that("each ring draws one arc per sector and level", {
    p <- ring_plot() + sc_dim_geom_ring() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase))
    arcs <- p$layers[[3]]$data

    expect_equal(nrow(arcs), 6L)
    expect_setequal(arcs$group, c("G1", "S"))
    ## the arcs of one sector tile it without leaving a slice out
    by_sector <- split(arcs, arcs$sector)
    expect_true(all(vapply(by_sector,
                           function(x) isTRUE(all.equal(
                               sum(x$end - x$start),
                               p$ggsc_ring$tab$span[p$ggsc_ring$tab$sector == x$sector[1]]
                           )),
                           logical(1))))
})

test_that("rings stack outwards instead of overlapping", {
    p <- ring_plot() + sc_dim_geom_ring() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase))
    inner <- p$layers[[2]]$data
    outer <- p$layers[[3]]$data

    expect_true(all(outer$r0 > inner$r1))
    expect_equal(p$ggsc_ring$r_max, unique(outer$r1))
})

test_that("the embedding is not rescaled when a ring is added", {
    p <- ring_plot()
    b0 <- ggplot2::ggplot_build(p)$data[[1]]
    b1 <- ggplot2::ggplot_build(p + sc_dim_geom_ring())$data[[1]]

    expect_equal(b0$x, b1$x)
    expect_equal(b0$y, b1$y)
})

test_that("the aspect ratio is fixed so the rings stay circular", {
    p <- ring_plot() + sc_dim_geom_ring()
    expect_true(isTRUE(p$coordinates$ratio == 1) ||
                inherits(p$coordinates, "CoordFixed"))
})

test_that("`r` overrides the radii and rejects a decreasing pair", {
    p <- ring_plot() + sc_dim_geom_ring(r = c(1, 1.5))
    expect_equal(unique(p$layers[[2]]$data$r1), 1.5)
    expect_error(ring_plot() + sc_dim_geom_ring(r = c(2, 1)), "increasing")
})

test_that("ring labels sit outside the outermost ring", {
    p <- ring_plot() + sc_dim_geom_ring() + sc_dim_geom_ring_label()
    lab <- p$layers[[3]]$data
    st <- p$ggsc_ring

    expect_setequal(lab$label, c("a", "b", "c"))
    r <- sqrt((lab$x - st$cx)^2 + (lab$y - st$cy)^2)
    expect_true(all(r > st$r_max))
})

test_that("the ring layers report what cannot be done", {
    p <- ring_plot()

    expect_error(p + sc_dim_geom_ring_label(), "no ring to label")
    expect_error(p + sc_dim_geom_ring(sector = "nope"), "sector column")
    expect_error(p + sc_dim_geom_ring(mapping = ggplot2::aes(colour = nope)),
                 "not in the plot data")
    expect_error(p + ggplot2::facet_wrap(~phase) + sc_dim_geom_ring(),
                 "faceted plots")

    ## a numeric column cannot sector the cells
    d <- p$data
    d$v <- seq_len(nrow(d))
    pv <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = v)) +
        ggplot2::geom_point()
    expect_error(pv + sc_dim_geom_ring(), "sector column")
})

## ---- legends -------------------------------------------------------------

## Names of the fill scales of a built plot (ggnewscale renames all but the
## last to `fill_ggnewscale_N`, so collect them from the scale list).
fill_scale_names <- function(p) {
    sc <- ggplot2::ggplot_build(p)$plot$scales$scales
    is_fill <- vapply(sc, function(s) any(grepl("^fill", s$aesthetics)), logical(1))
    vapply(sc[is_fill], function(s) s$name %||% "", character(1))
}

test_that("a ring repeating the point colour adds no legend by default", {
    p <- ring_plot() + sc_dim_geom_ring()
    expect_identical(p$ggsc_ring$n_fill_scales, 0L)
    expect_length(fill_scale_names(p), 0L)
})

test_that("a metadata ring gets its own fill scale and legend", {
    skip_if_not_installed("ggforce")
    p <- ring_plot() + sc_dim_geom_ring() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase))
    expect_identical(p$ggsc_ring$n_fill_scales, 1L)
    expect_identical(fill_scale_names(p), "phase")
})

test_that("legend.title overrides the legend name", {
    p <- ring_plot() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase),
                         legend.title = "Cell phase")
    expect_identical(fill_scale_names(p), "Cell phase")
})

test_that("legend = TRUE forces a legend even for the point colour", {
    p <- ring_plot() + sc_dim_geom_ring(legend = TRUE)
    expect_identical(p$ggsc_ring$n_fill_scales, 1L)
    expect_identical(fill_scale_names(p), "cluster")
})

test_that("legend = FALSE hides a metadata ring legend", {
    p <- ring_plot() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase), legend = FALSE)
    expect_identical(p$ggsc_ring$n_fill_scales, 0L)
    expect_length(fill_scale_names(p), 0L)
})

test_that("two legend-bearing rings keep two independent fill scales", {
    skip_if_not_installed("ggnewscale")
    d <- ring_plot()$data
    d$batch <- rep(c("b1", "b2", "b3"), length.out = nrow(d))
    p <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = cluster)) +
        ggplot2::geom_point()
    p <- p + sc_dim_geom_ring() +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = phase)) +
        sc_dim_geom_ring(mapping = ggplot2::aes(colour = batch))
    expect_identical(p$ggsc_ring$n_fill_scales, 2L)
    expect_setequal(fill_scale_names(p), c("phase", "batch"))
    expect_silent(ggplot2::ggplot_build(p))
})
