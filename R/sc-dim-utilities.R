##' @title sc_dim_count
##' @rdname sc-dim-count
##' @param sc_dim_plot dimension reduction plot of single cell data
##' @return a bar plot to present the cell numbers of different clusters
##' @importFrom ggplot2 coord_flip
##' @importFrom ggplot2 layer_data
##' @importFrom ggplot2 geom_col
##' @importFrom ggplot2 scale_y_continuous
##' @importFrom ggplot2 theme_minimal
##' @importFrom stats setNames
##' @seealso [sc_dim()]
##' @export
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p <- sc_dim(sce, reduction = 'UMAP')
##' p1 <- sc_dim_count(p)
sc_dim_count <- function(sc_dim_plot) {
  x <- layer_data(sc_dim_plot)


  if (length(sc_dim_plot$layers) >= 2) {
    ## assume label exists
    dd <- layer_data(sc_dim_plot, 2)

    pos <- dplyr::group_by(x, .data$colour) |>
      dplyr::summarize(x=mean(.data$x), y=mean(.data$y))
    dd <- unique(dd[, c("x", "y", "label")])
    # idx <- match(paste(dd$x, dd$y), paste(pos$x, pos$y)) # find closest one is better
    idx <- vapply(seq_len(nrow(dd)), function(i) {
                which.min((dd$x[i] - pos$x)^2 + (dd$y[i] - pos$y)^2)
            }, FUN.VALUE = numeric(1)
        )
    dd$colour <- pos$colour[idx]
    y <- setNames(dd$label, dd$colour)
  } else {
    d2 <- unique(x[, c("colour", "group")])
    y <- setNames(d2$group, d2$colour)
  }

  d <- as.data.frame(sort(table(x$colour)))
  d$group <- y[as.character(d$Var1)] |> as.factor()
  
  rlang::check_installed("forcats", "for sc_dim_count()")

  ggplot(d, 
    aes(forcats::fct_rev(.data$group), 
        .data$Freq, 
        fill = I(as.character(.data$Var1)))) + 
    geom_col() + coord_flip() +
    scale_y_continuous(expand=c(0,0)) +
    theme_minimal() +
    xlab(NULL) +
    ylab(NULL)
}


##' @title sc_dim_geom_feature
##' @rdname sc-dim-geom-feature
##' @param object Seurat or SingleCellExperiment object
##' @param features selected features (i.e., genes)
##' @param dims selected dimensions (must be a two-length vector) that are used 
##' in visualization
##' @param ncol number of facet columns if 'length(features) > 1'
##' @param .fun user defined function that will be applied to selected features 
##' (default is to filter out genes with no expression values)
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return layer of points for selected features
##' @export
##' @seealso [sc_feature()]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p1 <- sc_dim(sce, reduction = 'UMAP')
##' set.seed(123)
##' genes <- rownames(sce) |> sample(6)
##' f1 <- p1 + 
##'       sc_dim_geom_feature(
##'         object = sce, 
##'         features = genes
##'       )
sc_dim_geom_feature <- function(object, features, dims = c(1,2), ncol=3, ..., 
            .fun=function(.data) dplyr::filter(.data, .data$value > 0)) {
    params <- list(...)
    structure(
      list(data=object, features=features, 
           dims=dims, ncol=ncol, params=params, .fun=.fun), 
       class = 'sc_dim_geom_feature'
    )
    #d <- get_dim_data(object, dims=dims, features=features)
    #d <- tidyr::pivot_longer(d, 4:ncol(d), names_to = "features")
    #d$features <- factor(d$features, levels = features)
    #p <- sc_geom_point(data = .fun(d), ...)
    #list(p, 
    #    .feature_setting(features=features, ncol=ncol)
    #)
}

##' @method ggplot_add sc_dim_geom_feature
##' @importFrom tibble as_tibble
##' @export
ggplot_add.sc_dim_geom_feature <- function(object, plot, object_name, ...){
    if (inherits(object$data, 'Seurat')){
        d <- get_dim_data(object$data, 
                          dims=object$dims, 
                          features=object$features)
    }else{
        d <- .extract_sce_data(object$data, 
                               dims = object$dims, 
                               features = object$features)
    }
    #d <- as_tibble(d, rownames='.ID.NAME')

    d <- tidyr::pivot_longer(
           d, 
           seq(ncol(d) - length(object$features) + 1, ncol(d)), 
           names_to = "features"
         ) |> 
         dplyr::select(-c(2, 3)) |>
         dplyr::left_join(plot$data[,seq_len(3)],
                          by='.BarcodeID'
         )
    if (is.numeric(object$features)){
        object$features <- rownames(object$data)[object$features]
    }
    d$features <- factor(d$features, levels = object$features)
    d <- object$.fun(d)
    sc.point.params <- object$params
    sc.point.params$data <- d
    geomfun <- .extract_geom_name(plot)
    if (geomfun == 'geom_scattermore2'){
        geomfun <- "sc_geom_point"
    }else{
        sc.point.params$pixels <- NULL
    }    
    p <- do.call(geomfun, sc.point.params)
    ly <- list(p,
        .feature_setting(features=object$features, ncol=object$ncol)
    )
    ggplot_add(ly, plot, object_name, ...)
}


##' @title sc_dim_geom_label
##' @rdname sc-dim-geom-label
##' @param geom geometric layer (default: geom_text) to display the lables
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to the geom
##' @return layer of labels
##' @export
##' @seealso [sc_dim_geom_label()]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p1 <- sc_dim(sce, reduction = 'UMAP', mapping = aes(colour = Cell_Cycle))
##' p2 <- sc_dim(sce, reduction = 'UMAP')
##' f1 <- p1 + sc_dim_geom_label()
sc_dim_geom_label <- function(geom = ggplot2::geom_text, mapping=NULL, ...) {
    structure(list(geom = geom, mapping = mapping, ...),
        class = "sc_dim_geom_label")
}

##' @importFrom ggplot2 ggplot_add
##' @importFrom rlang .data
##' @method ggplot_add sc_dim_geom_label
##' @export
ggplot_add.sc_dim_geom_label <- function(object, plot, object_name, ...) {
    dims <- names(plot$data)[seq_len(3)]
    if (!is.null(object$mapping$label)){
        lab.text <- ggfun::get_aes_var(object$mapping, 'label')
        object$mapping$label <- NULL
    }else{
        lab.text <- ggplot_build(plot)$plot$labels$colour
    }
    flag1 <- lab.text %in% colnames(plot$data) && !is.numeric(plot$data[[lab.text]]) 
    if (is.null(object$data) && flag1) {
        object$data <- split(plot$data, plot$data[[lab.text]]) |> 
            lapply(function(x).calculate_ellipse(x, vars = dims[c(2, 3)], level=object$level)) |>
            dplyr::bind_rows(.id=lab.text) 
        object$level <- NULL
        object$data <- .set_label_levels(object$data, plot, lab.text)
    }else{
        cli::cli_abort("The `label` in mapping should be specified, and the data should not be numeric type!")
    }

    geom <- object$geom
    object$geom <- NULL
    default_mapping <- aes(x=!!rlang::sym(dims[2]), y = !!rlang::sym(dims[3]), label = !!rlang::sym(lab.text))
    if (is.null(object$mapping)) {
        object$mapping <- default_mapping
    } else {
        object$mapping <- utils::modifyList(default_mapping, object$mapping)
    }

    flag2 <- .check_colour(plot, object)
    if (flag2){
        object$colour <- 'black'
    }
    
    object <- .set_inherit.aes(object)

    ly <- do.call(geom, object)    
    ggplot_add(ly, plot, object_name, ...)
}



##' @title sc_dim_geom_ellipse
##' @rdname sc-dim-geom-ellipse
##' @param geom the layer function, default is \code{stat_ellipse},
##' other option is \code{geom_mark_hull} of \code{ggforce}.
##' @param mapping aesthetic mapping
##' @param level the level at which to draw an ellipse
##' @param ... additional parameters pass to the stat_ellipse
##' @return layer of ellipse
##' @seealso
##'  [stat_ellipse][ggplot2::stat_ellipse]; 
##' @export
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p1 <- sc_dim(sce, reduction = 'UMAP', mapping = aes(colour = Cell_Cycle))
##' p2 <- sc_dim(sce, reduction = 'UMAP')
##' f1 <- p1 + sc_dim_geom_ellipse()
sc_dim_geom_ellipse <- function(geom = stat_ellipse, mapping = NULL, level = 0.95, ...) {
    structure(list(geom = geom, mapping = mapping, level = level, ...), 
              class = "sc_dim_geom_ellipse")
}

##' @importFrom ggplot2 ggplot_add
##' @importFrom rlang .data
##' @method ggplot_add sc_dim_geom_ellipse
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 stat_ellipse
##' @export
ggplot_add.sc_dim_geom_ellipse <- function(object, plot, object_name, ...) {
    dims <- names(plot$data)[seq_len(3)]
    if (!is.null(object$mapping$group)){
        lab.text <- ggfun::get_aes_var(object$mapping, 'group')
        object$mapping$group <- NULL
    }else{
        lab.text <- ggplot_build(plot)$plot$labels$colour
    }
    flag1 <- lab.text %in% colnames(plot$data) && !is.numeric(plot$data[[lab.text]])
    if (!flag1){
        cli::cli_abort("The `group` in mapping should be specified, and the data should not be numeric type!")
    }
    default_mapping <- aes(x = !!rlang::sym(dims[2]), 
                           y = !!rlang::sym(dims[3]), 
                           group = !!rlang::sym(lab.text))
    if (is.null(object$mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, object$mapping)
    }
    object$mapping <- mapping
    
    flag2 <- .check_colour(plot, object)
    if (flag2){
        object$colour <- 'black'
    }
    object <- .set_inherit.aes(object)    
    geomfun <- object$geom
    object$geom <- NULL

    ly <- do.call(geomfun, object)
    ggplot_add(ly, plot, object_name, ...)
}

##' @title sc_dim_geom_subset
##' @rdname sc-dim-geom-subset
##' @param mapping aesthetic mapping
##' @param subset subset of clusters to be displayed
##' @param .column which column represents cluster (e.g., 'ident')
##' @param ... additional parameters pass to sc_geom_point
##' @return plot with a layer of specified clusters
##' @export
##' @seealso [sc_dim_geom_sub]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p1 <- sc_dim(sce, reduction = 'UMAP')
##' f1 <- p1 + sc_dim_geom_sub(subset = c(1, 2), .column = 'label', bg_colour='black')
sc_dim_geom_sub <- function(mapping = NULL, subset, .column = "ident", ...) {
  structure(list(mapping = mapping, 
        subset = subset, 
        .column = .column,
        ...), 
    class = "dim_geom_sub")
}

##' @method ggplot_add dim_geom_sub
##' @export     
ggplot_add.dim_geom_sub <- function(object, plot, object_name, ...) {
  ii <- plot$data[[object$.column]] %in% object$subset
  object$data <- plot$data[ii, ]
  default_mapping <- aes(color = .data[[object$.column]])
  if (is.null(object$mapping)) {
    mapping <- default_mapping
  } else {
    mapping <- modifyList(default_mapping, object$mapping)
  }
  object$mapping <- mapping
  object$subset <- NULL
  object$.column <- NULL
  geomfun <- .extract_geom_name(plot)
  if (geomfun == 'geom_scattermore2'){
      geomfun <- "sc_geom_point"
  }else{
      object$pixels <- NULL
  }
  ly <- do.call(geomfun, object)
  ggplot_add(ly, plot, object_name, ...)
}

##' @title sc_dim_sub
##' @rdname sc-dim-sub
##' @param subset subset of clusters to be displayed
##' @param .column which column represents cluster (e.g., 'ident')
##' @return update plot with only subset displayed
##' @export
##' @seealso [sc_dim]
##' @examples
##' library(scuttle)
##' library(scater)
##' library(scran)
##' library(ggplot2)
##' sce <- mockSCE()
##' sce <- logNormCounts(sce)
##' clusters <- clusterCells(sce, assay.type = 'logcounts')
##' colLabels(sce) <- clusters
##' sce <- runUMAP(sce, assay.type = 'logcounts')
##' p1 <- sc_dim(sce, reduction = 'UMAP')
##' f1 <- p1 + sc_dim_sub(subset = c(1, 2), .column = 'label')
sc_dim_sub <- function(subset, .column = "ident") {
  structure(list(subset = subset, .column = .column), class = "dim_sub")
}

##' @method ggplot_add dim_sub
##' @export     
ggplot_add.dim_sub <- function(object, plot, object_name, ...) {
  ii <- plot$data[[object$.column]] %in% object$subset
  plot$data <- plot$data[ii, ]
  plot
}
