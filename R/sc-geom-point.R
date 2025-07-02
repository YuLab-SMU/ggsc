##' @title sc_geom_point
##' @rdname sc-geom-point
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return layer of points
##' @importFrom scattermore geom_scattermore
##' @seealso [sc_dim()] and [sc_feature()]
##' @export
##' @examples
##' library(ggplot2)
##' ggplot(iris,
##'  aes(x= Sepal.Length, y = Petal.Width, color=Species)
##' ) + 
##' sc_geom_point()
sc_geom_point <- function(mapping=NULL, ...){
    default_params <- list(mapping = mapping, 
                        pointsize = 3,
                        pixels = c(700, 700)
                    )
    params <- modifyList(default_params, list(...))
    do.call(geom_scattermore2, params)
}

#' @title geom_bgpoint
#' @description
#' this add the background color for \code{\link[ggplot2]{geom_point}}
#' @inheritParams ggplot2::layer
#' @param na.rm If \code{FALSE}, the default, missing values are removed 
#' with a warning, if \code{TRUE}, missing values are silently removed.
#' @param gap_colour colour of gap background between the bottom background
#' and top point point layer, default is \code{white}.
#' @param gap_alpha numeric the transparency of gap background colour, default is 1.
#' @param bg_line_width numeric the line width of background point layer,
#' default is \code{0.3}.
#' @param gap_line_width numeric the line width of gap between the background and
#' top point point layer, default is \code{.1}.
#' @param pointsize numeric the size of point, default is NULL, will use the 
#' internal size aesthetics of \code{geom_bgpoint}
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}.
#' @details Aesthetics \code{geom_bgpoint} understands the following aesthetics. Required
#' aesthetics are displayed in bold and defaults are displayed for optional aesthetics:
#'  \itemize{
#'     \item \strong{\code{x}}.
#'     \item \strong{\code{y}}.
#'     \item \code{colour} the colour of point, default is \code{black}.
#'     \item \code{bg_colour} the colour of background point, default is \code{NA}.
#'     \item \code{alpha} the transparency of colour, default is 1.
#'     \item \code{subset} subset the data frame which meet conditions to display.
#'  }
#' @importFrom rlang list2
#' @author Shuangbin Xu
#' @export
#' @examples
##' library(ggplot2)
##' ggplot(iris,
##'  aes(x= Sepal.Length, y = Petal.Width, color=Species, bg_colour=Species)
##' ) +
##' geom_bgpoint(pointsize=4, gap_line_width = .1, bg_line_width = .3)
geom_bgpoint <- function(
    mapping = NULL, 
    data = NULL,
    stat = "identity", 
    position = "identity",
    ...,
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    gap_colour = 'white', 
    gap_alpha = 1, 
    bg_line_width = .3, 
    gap_line_width = .1,
    pointsize = NULL
){

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomBgpoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list2(
      na.rm = na.rm,
      gap_colour = gap_colour,
      gap_alpha = gap_alpha,
      bg_line_width = bg_line_width,
      gap_line_width = gap_line_width,  
      pointsize = pointsize,
      ...
    )
  )
}


#' @title geom_scattermore2
#' @description
#' this add the background colour for the \code{\link[scattermore]{geom_scattermore}} 
#' @inheritParams ggplot2::layer
#' @param na.rm If `FALSE`, the default, missing values are removed with
#'     a warning. If `TRUE`, missing values are silently removed.
#' @param interpolate A logical value indicating whether to linearly interpolate
#' the image (the alternative is to use nearest-neighbour interpolation, 
#' which gives a more blocky result). Default \code{FALSE}, 
#' passed to \code{\link[grid]{rasterGrob}}.
#' @param pointsize Radius of rasterized point. Use ‘0’ for single pixels (fastest).
#' @param pixels Vector with X and Y resolution of the raster, default \code{c(512,512)}.
#' @param gap_colour colour of gap background between the bottom background 
#' and top point point layer, default is \code{white}. 
#' @param gap_alpha numeric the transparency of gap background colour, default is 1.
#' @param bg_line_width numeric the line width of background point layer, 
#' default is \code{0.3}.
#' @param gap_line_width numeric the line width of gap between the background and 
#' top point point layer, default is \code{.1}.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}.
#' @details Aesthetics \code{geom_scattermore2} understands the following aesthetics. Required
#' aesthetics are displayed in bold and defaults are displayed for optional aesthetics:
#'  \itemize{
#'     \item \strong{\code{x}}.
#'     \item \strong{\code{y}}.
#'     \item \code{colour} the colour of point, default is \code{black}.
#'     \item \code{bg_colour} the colour of background point, default is \code{NA}.
#'     \item \code{alpha} the transparency of colour, default is 1.
#'     \item \code{subset} subset the data frame which meet conditions to display.
#'  }
#' @return polygonal point layer
#' @importFrom ggplot2 layer
#' @author Shuangbin Xu
#' @export
#' @examples
##' library(ggplot2)
##' ggplot(iris,
##'  aes(x= Sepal.Length, y = Petal.Width, color=Species, bg_colour=Species)
##' ) + 
##' geom_scattermore2(pointsize=4, gap_line_width = .1, bg_line_width = .3)
geom_scattermore2 <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", ...,
                             na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
                             interpolate = FALSE, pointsize = 0, pixels = c(512, 512),
                             gap_colour = 'white', gap_alpha = 1, bg_line_width = .3, gap_line_width = .1){
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    position = position,
    geom = GeomScattermore2,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list2(
      na.rm = na.rm,
      interpolate = interpolate,
      pointsize = pointsize,
      pixels = pixels,
      gap_colour = gap_colour,
      gap_alpha = gap_alpha,
      bg_line_width = bg_line_width,
      gap_line_width = gap_line_width,
      ...
    )
  )
}


#' @importFrom scattermore scattermore
GeomScattermore2 <- ggplot2::ggproto("GeomScattermore2", ggplot2::Geom,
  required_aes = c("x", "y"),
  non_missing_aes = c("alpha", "colour"),
  optional_aes = c("subset"),
  default_aes = ggplot2::aes(
    shape = 19, colour = "black", size = 1.5, fill = NA,
    alpha = 1, stroke = 0.5, bg_colour = NA
  ),
  setup_data = function(data, params){
      if (is.null(data$subset))
          return(data)
      data[which(data$subset),]
  },                                     
  draw_panel = function(data, pp, coord,
                        pointsize = 0,
                        interpolate = FALSE,
                        na.rm = FALSE,
                        pixels = c(512, 512),
                        gap_colour = 'white',
                        gap_alpha = 1,
                        bg_line_width = .3,
                        gap_line_width = .1){
    coords <- coord$transform(data, pp)

    upperimage <- scattermore(cbind(coords$x, coords$y),
                              rgba = grDevices::col2rgb(alpha = TRUE, scales::alpha(coords$colour, coords$alpha)),
                              cex = pointsize,
                              xlim = c(0, 1),
                              ylim = c(0, 1),
                              size = pixels)

    if (!all(is.null(coords$bg_colour)) && !all(is.na(coords$bg_colour))){
        tmpsize <- sqrt(pointsize)
        gapsize <- (tmpsize + tmpsize * gap_line_width * 2)^2
        bgsize <- gapsize + (sqrt(bg_line_width) + tmpsize * bg_line_width * 2) ^2
        bgimage <- scattermore(cbind(coords$x, coords$y),
                               rgba = grDevices::col2rgb(alpha = TRUE, scales::alpha(coords$bg_colour, 1)),
                               cex = bgsize,
                               xlim = c(0, 1),
                               ylim = c(0, 1),
                               size = pixels)
        gapimage <- scattermore(cbind(coords$x, coords$y),
                                rgba = grDevices::col2rgb(alpha = TRUE, scales::alpha(gap_colour, gap_alpha)),
                                cex = gapsize,
                                xlim = c(0, 1),
                                ylim = c(0, 1),
                                size = pixels)

    }else{
        bgimage <- gapimage <- NULL
    }

    ggplot2:::ggname(
      "geom_scattermore2",
      crasterGrob(
        upperimage, bgimage, gapimage,
        0, 0, 1, 1,
        default.units = "native",
        just = c("left", "bottom"),
        interpolate = interpolate
      )
    )
  },
  draw_key = draw_key_bgpoint
)

#' @importFrom ggplot2 translate_shape_string draw_key_point
GeomBgpoint <- ggplot2::ggproto("GeomBgpoint", ggplot2::Geom,
  required_aes = c("x", "y"),
  non_missing_aes = c("alpha", "colour"),
  optional_aes = c("subset"),
  default_aes = ggplot2::aes(
    shape = 19, colour = "black", size = 1.5, fill = NA,
    alpha = 1, stroke = 0.5, bg_colour = NA
  ),
  setup_data = function(data, params){
      if (is.null(data$subset))
          return(data)
      data[which(data$subset),]
  },
  draw_panel = function(self, 
                        data,
                        panel_params,
                       	coord,
                       	na.rm = FALSE,
                        gap_colour = 'white',
                        gap_alpha = 1,
                        bg_line_width = .3,
                        gap_line_width = .1,
                        pointsize = NULL 
			){
    if (is.character(data$shape)) {
      data$shape <- translate_shape_string(data$shape)
    }

    coords <- coord$transform(data, panel_params)
    stroke_size <- coords$stroke
    stroke_size[is.na(stroke_size)] <- 0
    if (is.null(pointsize)){
       pointsize <- coords$size * .pt + stroke_size * .stroke / 2
    }else{
       pointsize <- rep(pointsize, nrow(coords)) * .pt + stroke_size * .stroke/2
    }
    ggname("geom_bgpoint",
      cpointsGrob(
        coords$x, coords$y,
        pch = coords$shape,
        bg_line_width = bg_line_width,
        gap_line_width = gap_line_width,
        bg_colour = coords$bg_colour,
        gap_colour = alpha(gap_colour, gap_alpha),
        gp = gpar(
          col = alpha(coords$colour, coords$alpha),
          fill = fill_alpha(coords$fill, coords$alpha),
          fontsize = pointsize,
          lwd = coords$stroke * .stroke / 2
        )
      )
    )
  },
  draw_key = draw_key_bgpoint
)


#' @importFrom grid grobName
ggname <- function (prefix, grob){
   grob$name <- grobName(grob, prefix)
   grob
}


