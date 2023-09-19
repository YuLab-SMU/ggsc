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
sc_geom_point <- function(mapping=NULL, ...) {
    default_params <- list(mapping = mapping, 
                        pointsize = 2, 
                        pixels = c(700, 700)
                    )
    params <- modifyList(default_params, list(...))
    do.call(scattermore::geom_scattermore, params)
}
