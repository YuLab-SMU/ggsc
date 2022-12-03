##' @title sc_dim_geom_feature
##' @rdname sc-dim-geom-feature
##' @param object Seurat object
##' @param features selected features (i.e., genes)
##' @param dims selected dimensions (must be a two-length vector) that are used in visualization
##' @param ncol number of facet columns if 'length(features) > 1'
##' @param .fun user defined function that will be applied to selected features (default is to filter out genes with no expression values)
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return layer of points for selected features
##' @export
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
ggplot_add.sc_dim_geom_feature <- function(object, plot, object_name){
    d <- get_dim_data(object$data, dims=object$dims, features=object$features)
    d <- as_tibble(d, rownames='.ID.NAME')
    d <- tidyr::pivot_longer(d, 5:ncol(d), names_to = "features") |> 
         dplyr::select(-c(2, 3, 4)) |>
         dplyr::left_join(plot$data[,seq_len(3)] |>
                          tibble::as_tibble(rownames='.ID.NAME'), 
                          by='.ID.NAME'
         )
    d$features <- factor(d$features, levels = object$features)
    d <- object$.fun(d)
    sc.point.params <- object$params
    sc.point.params$data <- d
    p <- do.call(sc_geom_point, sc.point.params)
    ly <- list(p,
        .feature_setting(features=object$features, ncol=object$ncol)
    )
    ggplot_add(ly, plot, object_name)
}


##' @title sc_dim_geom_label
##' @rdname sc-dim-geom-label
##' @param geom geometric layer (default: geom_text) to display the lables
##' @param ... additional parameters pass to the geom
##' @return layer of labels
##' @export
sc_dim_geom_label <- function(geom = ggplot2::geom_text, ...) {
    structure(list(geom = geom, ...),
        class = "sc_dim_geom_label")
}

##' @importFrom ggplot2 ggplot_add
##' @importFrom rlang .data
##' @method ggplot_add sc_dim_geom_label
##' @export
ggplot_add.sc_dim_geom_label <- function(object, plot, object_name) {
    dims <- names(plot$data)[1:2]
    if (is.null(object$data)) {
        object$data <- dplyr::group_by(plot$data, .data$ident) |> 
            dplyr::summarize(x=mean(get(dims[1])), y=mean(get(dims[2])))
    }

    geom <- object$geom
    object$geom <- NULL
    default_mapping <- aes(x=.data$x, y = .data$y, label = .data$ident)
    if (is.null(object$mapping)) {
        object$mapping <- default_mapping
    } else {
        object$mapping <- utils::modifyList(default_mapping, object$mapping)
    }
    ly <- do.call(geom, object)    
    ggplot_add(ly, plot, object_name)
}



##' @title sc_dim_geom_ellipse
##' @rdname sc-dim-geom-ellipse
##' @param mapping aesthetic mapping
##' @param level the level at which to draw an ellipse
##' @param ... additional parameters pass to the stat_ellipse
##' @return layer of ellipse
##' @seealso
##'  [stat_ellipse][ggplot2::stat_ellipse]; 
##' @export
sc_dim_geom_ellipse <- function(mapping = NULL, level = 0.95, ...) {
    structure(list(mapping = mapping, level = level, ...), class = "sc_dim_geom_ellipse")
}

##' @importFrom ggplot2 ggplot_add
##' @importFrom rlang .data
##' @method ggplot_add sc_dim_geom_ellipse
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 stat_ellipse
##' @export                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ggplot_add.sc_dim_geom_ellipse <- function(object, plot, object_name) {
    dims <- names(plot$data)[1:2]
    default_mapping <- aes(x = .data[[dims[1]]], y = .data[[dims[2]]], group = .data$ident)
    if (is.null(object$mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, object$mapping)
    }
    object$mapping <- mapping

    ly <- do.call(stat_ellipse, object)
    ggplot_add(ly, plot, object_name)
}

##' @title sc_dim_geom_subset
##' @rdname sc-dim-geom-subset
##' @param mapping aesthetic mapping
##' @param subset subset of clusters to be displayed
##' @param .column which column represents cluster (e.g., 'ident')
##' @param ... additional parameters pass to sc_geom_point
##' @return plot with a layer of specified clusters
##' @export
sc_dim_geom_sub <- function(mapping = NULL, subset, .column = "ident", ...) {
  structure(list(mapping = mapping, 
        subset = subset, 
        .column = .column,
        ...), 
    class = "dim_geom_sub")
}

##' @method ggplot_add dim_geom_sub
##' @export     
ggplot_add.dim_geom_sub <- function(object, plot, object_name) {
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
  ly <- do.call(sc_geom_point, object)
  ggplot_add(ly, plot, object_name)
}

##' @title sc_dim_sub
##' @rdname sc-dim-sub
##' @param subset subset of clusters to be displayed
##' @param .column which column represents cluster (e.g., 'ident')
##' @return update plot with only subset displayed
##' @export
sc_dim_sub <- function(subset, .column = "ident") {
  structure(list(subset = subset, .column = .column), class = "dim_sub")
}

##' @method ggplot_add dim_sub
##' @export     
ggplot_add.dim_sub <- function(object, plot, object_name) {
  ii <- plot$data[[object$.column]] %in% object$subset
  plot$data <- plot$data[ii, ]
  plot
}
