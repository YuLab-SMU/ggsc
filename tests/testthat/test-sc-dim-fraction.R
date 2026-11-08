## Coverage for `sc_dim_fraction()`, the bar-plot companion of the rings.
## A hand-built plot keeps the tests fast and free of scater/scuttle.

fraction_plot <- function(n = 20) {
    set.seed(1)
    d <- data.frame(
        .BarcodeID = sprintf("c%03d", seq_len(3 * n)),
        UMAP1 = rnorm(3 * n),
        UMAP2 = rnorm(3 * n),
        cluster = factor(rep(c("a", "b", "c"), each = n), levels = c("c", "a", "b")),
        treat = rep(c("ctrl", "drug"), length.out = 3 * n),
        stringsAsFactors = FALSE
    )
    ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = cluster)) +
        ggplot2::geom_point()
}

test_that("sc_dim_fraction() returns a bar plot of the composition", {
    p <- fraction_plot()
    pf <- sc_dim_fraction(p, group.by = "treat")

    expect_s3_class(pf, "ggplot")
    expect_setequal(names(pf$data), c("celltype", "group", "n"))
    expect_equal(sum(pf$data$n), nrow(p$data))
    expect_setequal(as.character(pf$data$group), c("ctrl", "drug"))
})

test_that("position = 'fill' normalises each bar to proportions", {
    p <- fraction_plot()
    pf <- sc_dim_fraction(p, group.by = "treat", position = "fill")
    built <- ggplot2::ggplot_build(pf)$data[[1]]

    top <- tapply(built$ymax, built$x, max)
    expect_true(all(abs(top - 1) < 1e-8))
    expect_identical(pf$labels$y, "Fraction of cells")
})

test_that("position = 'stack' keeps raw counts", {
    p <- fraction_plot()
    pf <- sc_dim_fraction(p, group.by = "treat", position = "stack")
    built <- ggplot2::ggplot_build(pf)$data[[1]]

    top <- tapply(built$ymax, built$x, max)
    expect_equal(as.numeric(top), rep(nrow(p$data) / 2, 2))
    expect_identical(pf$labels$y, "Cell count")
})

test_that("the fill colours are taken from the embedding's colour scale", {
    p <- fraction_plot()
    pf <- sc_dim_fraction(p, group.by = "treat")

    cp <- ggplot2::ggplot_build(p)$plot$scales$get_scales("colour")
    cf <- ggplot2::ggplot_build(pf)$plot$scales$get_scales("fill")
    lv <- cp$get_limits()
    expect_identical(unname(cp$map(lv)), unname(cf$map(lv)))
    ## the cell-type order follows the embedding, not alphabetical order
    expect_identical(levels(pf$data$celltype), lv)
})

test_that("celltype can be a different metadata column", {
    p <- fraction_plot()
    pf <- sc_dim_fraction(p, group.by = "cluster", celltype = "treat")
    expect_setequal(as.character(pf$data$celltype), c("ctrl", "drug"))
    expect_identical(pf$labels$fill, "treat")
})

test_that("sc_dim_fraction() reports unusable columns", {
    p <- fraction_plot()

    expect_error(sc_dim_fraction(p, group.by = "nope"), "group.by")
    expect_error(sc_dim_fraction(p, group.by = "treat", celltype = "nope"), "celltype")

    d <- p$data
    d$num <- seq_len(nrow(d))
    pn <- ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2, colour = cluster)) +
        ggplot2::geom_point()
    expect_error(sc_dim_fraction(pn, group.by = "treat", celltype = "num"),
                 "categorical")
})
