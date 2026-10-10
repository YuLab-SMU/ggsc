## Regression coverage for `.dim_group_col()`, the helper that resolves the
## grouping column of `sc_dim_geom_ellipse()` / `sc_dim_geom_label()`.
##
## A hand-built plot is used instead of a real `sc_dim()` figure so the tests
## stay fast and do not need the scater/scuttle toolchain; a single
## end-to-end check with `scater::mockSCE()` (skipped when unavailable) covers
## the original report as well.

## `sc_dim()` with no `colLabels()`: the barcode column plus the two reduction
## coordinates and no non-numeric column mapped to `colour`.
dim_group_plot <- function(n = 30, extra = NULL) {
    set.seed(1)
    d <- data.frame(
        .BarcodeID = sprintf("cell%03d", seq_len(2 * n)),
        UMAP1 = c(rnorm(n), rnorm(n) + 4),
        UMAP2 = rnorm(2 * n),
        stringsAsFactors = FALSE
    )
    if (!is.null(extra)) {
        d[names(extra)] <- extra
    }
    ggplot2::ggplot(d, ggplot2::aes(x = UMAP1, y = UMAP2)) +
        ggplot2::geom_point()
}

dim_group_layer <- function(key = "group") {
    structure(list(mapping = NULL), class = "sc_dim_geom_ellipse")
}

test_that(".dim_group_col() reports a missing colour mapping instead of failing", {
    p <- dim_group_plot()
    grp <- .dim_group_col(dim_group_layer(), p, "group")

    ## `plot$data[[NULL]]` used to raise
    ## "attempt to select less than one element in get1index"
    expect_false(grp$ok)
    expect_null(grp$lab.text)
})

test_that(".dim_group_col() reports a label that is not a column of the data", {
    p <- dim_group_plot()
    obj <- dim_group_layer()
    obj$mapping <- ggplot2::aes(group = no_such_column)
    grp <- .dim_group_col(obj, p, "group")

    expect_false(grp$ok)
    expect_null(grp$lab.text)
})

test_that(".dim_group_col() reports a numeric grouping column", {
    p <- dim_group_plot(extra = list(v = seq_len(60)))
    p <- p + ggplot2::aes(colour = v)

    grp <- .dim_group_col(dim_group_layer(), p, "group")
    expect_false(grp$ok)
    expect_null(grp$lab.text)
})

test_that(".dim_group_col() accepts a non-numeric column and drops the mapping", {
    p <- dim_group_plot(extra = list(cluster = rep(c("a", "b"), each = 30)))
    obj <- dim_group_layer()
    obj$mapping <- ggplot2::aes(group = cluster)

    grp <- .dim_group_col(obj, p, "group")
    expect_true(grp$ok)
    expect_identical(grp$lab.text, "cluster")

    ## the aesthetic is taken over by the layer's own default mapping
    expect_null(grp$object$mapping$group)
})

test_that("sc_dim_geom_ellipse() explains a plot with no grouping column", {
    p <- dim_group_plot()
    expect_error(p + sc_dim_geom_ellipse(),
                 "The group column could not be determined")
})

test_that("sc_dim_geom_label() explains a plot with no grouping column", {
    p <- dim_group_plot()
    expect_error(p + sc_dim_geom_label(),
                 "The label column could not be determined")
})

test_that("sc_dim_geom_label() explains a numeric label column", {
    p <- dim_group_plot(extra = list(v = seq_len(60)))
    p <- p + ggplot2::aes(colour = v)
    expect_error(p + sc_dim_geom_label(),
                 "The label column could not be determined")
})

test_that("an explicit `mapping = aes(group = )` still builds", {
    p <- dim_group_plot(extra = list(cluster = rep(c("a", "b"), each = 30)))
    built <- ggplot2::ggplot_build(
        p + sc_dim_geom_ellipse(mapping = ggplot2::aes(group = cluster))
    )
    expect_equal(length(unique(built$data[[2]]$group)), 2L)
})

test_that("an explicit `mapping = aes(label = )` still builds", {
    p <- dim_group_plot(extra = list(cluster = rep(c("a", "b"), each = 30)))
    built <- ggplot2::ggplot_build(
        p + sc_dim_geom_label(mapping = ggplot2::aes(label = cluster))
    )
    expect_setequal(built$data[[2]]$label, c("a", "b"))
})

test_that("the plot's non-numeric colour column still groups both layers", {
    p <- dim_group_plot(extra = list(cluster = rep(c("a", "b"), each = 30)))
    p <- p + ggplot2::aes(colour = cluster)

    expect_equal(
        length(unique(ggplot2::ggplot_build(p + sc_dim_geom_ellipse())$data[[2]]$group)),
        2L
    )
    expect_setequal(
        ggplot2::ggplot_build(p + sc_dim_geom_label())$data[[2]]$label,
        c("a", "b")
    )
})

test_that("the reported `sc_dim()` + no `colLabels()` figure aborts cleanly", {
    skip_if_not_installed("scater")
    skip_if_not_installed("scuttle")

    set.seed(1)
    sce <- scuttle::mockSCE()
    sce <- suppressWarnings(scuttle::logNormCounts(sce))
    sce <- suppressWarnings(scater::runUMAP(sce, assay.type = "logcounts"))
    p <- suppressWarnings(sc_dim(sce, reduction = "UMAP"))

    expect_error(p + sc_dim_geom_ellipse(),
                 "The group column could not be determined")
    expect_error(p + sc_dim_geom_label(),
                 "The label column could not be determined")
})
