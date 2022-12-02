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
    d <- get_dim_data(object, dims=dims, features=features)
    d <- tidyr::pivot_longer(d, 4:ncol(d), names_to = "features")
    d$features <- factor(d$features, levels = features)
    p <- sc_geom_point(data = .fun(d), ...)
    list(p, 
        .feature_setting(features=features, ncol=ncol)
    )
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
##' @param mapping, aesthetic mapping
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

