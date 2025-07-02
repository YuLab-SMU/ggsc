#' @title add the annotation layer for ggsc object
#' @param data The data to be displayed in this layer. There are three
#' options:
#'  If \code{NULL}, the default, the data is inherited from the plot
#'  data as specified in the call to \code{ggplot()}.
#'  A \code{data.frame}, will override the plot data. the \code{data.frame}
#'  should have a barcode id or features column.
#'  A \code{function} will be called with a single argument, the plot
#'  data. The return value must be a ‘data.frame’, and will be
#'  used as the layer data. A \code{function} can be created from a
#'  ‘formula’ (e.g. ‘~ head(.x, 10)’).
#' @param mapping Set of aesthetic mappings created by \code{aes()}. If specified
#' and \code{inherit.aes = TRUE} (the default), it is combined with the default 
#' mapping at the top level of the plot. You must supply \code{mapping} if there is no plot mapping.
#' @inheritParams geom_scattermore2
#' @export
#' @return layer object
sc_geom_annot <- function(
     data=NULL, 
     mapping=NULL,
     pointsize = 2, 
     pixels = c(512, 512),
     gap_colour = "white",
     gap_alpha = 1,
     bg_line_width = 0.3,
     gap_line_width = 0.1,
     show.legend = NA,
     na.rm = FALSE,
     ...
){
    params <- list(...)
    x <- structure(
      list(data = data, 
           mapping = mapping, 
	   pointsize = pointsize,
           pixels = pixels, 
           gap_colour = gap_colour,
           gap_alpha = gap_alpha,
           bg_line_width = bg_line_width,
           gap_line_width = gap_line_width,
           show.legend = show.legend,
           na.rm = na.rm, 
           params=params),
       class = 'sc_geom_annot'
    )
    return(x)
}

#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add sc_geom_annot
#' @export
ggplot_add.sc_geom_annot <- function(object, plot, object_name, ...){
    object <- .check_layer_data(object, plot)
    params <- object$params
    object$params <- NULL
    geomfun <- .extract_geom_name(plot)
    if (geomfun == 'geom_scattermore2'){
        geomfun <- "sc_geom_point"
    }else{
        object$pixels <- NULL
    }
    ly <- do.call(geomfun, c(object, params))
    ggplot_add(ly, plot, object_name, ...)
}


.check_layer_data <- function(object, plot){
    if (is.data.frame(object$data)){
        object$data <- plot$data |> 
             dplyr::left_join(object$data, suffix = c("", ".y")) |> 
             suppressMessages()
    }
    return(object)
}

.extract_geom_name <- function(plot){
    ind <- 1
    if (length(plot$layers)>1){
        ind <- 2
    }
    lay <- plot$layers[[ind]]
    x <- snakeize(class(lay$geom))[[1]]
    x <- gsub("new_","", x)
    return(x)
}

# this is from the internal function of ggplot2
snakeize <- function(x){
    x <- gsub("([A-Za-z])([A-Z])([a-z])", "\\1_\\2\\3", x)
    x <- gsub(".", "_", x, fixed = TRUE)
    x <- gsub("([a-z])([A-Z])", "\\1_\\2", x)
    to_lower_ascii(x)
}

upper_ascii <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
lower_ascii <- "abcdefghijklmnopqrstuvwxyz"
to_lower_ascii <- function(x){
   chartr(upper_ascii, lower_ascii, x)
}
