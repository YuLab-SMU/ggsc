##' @title sc_dim
##' @rdname sc-dim-methods
##' @param object Seurat object or SingleCellExperiment object
##' @param dims selected dimensions (must be a two-length vector) that 
##' are used in visualization
##' @param reduction reduction method, default is NULL and will use the 
##' default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping, the \code{x} and \code{y} is set internally,
##' other character of geometric layer, such as \code{color}, \code{size},
##' \code{alpha} or (\code{shape} when geom = geom_point) can be set manually.
##' @param geom the function of geometric layer, default is sc_geom_point,
##' other geometric layer, such as \code{geom_bgpoint} or \code{geom_point} also works.
##' @param ... additional parameters pass to \code{geom_scattermore2()}.
##' \itemize{
##'     \item \code{bg_colour} the colour of background point, default is \code{NA}.
##'      this character also can be set in \code{mappint}.
##'     \item \code{gap_colour} the colour of gap background, default is \code{'white'}.
##'     \item \code{bg_line_width} the line width of background point, 
##'      default is \code{.3}.
##'     \item \code{gap_line_width} the gap line width of background point, 
##'      default is \code{.1}.
##'     \item \code{alpha} the transparency of colour, default is 1.
##'     \item \code{subset} subset the data frame which meet conditions to display.
##'      this should be set in \code{mapping}.
#'  }
##' @return dimension reduction plot
##' @seealso
##'  [geom_scattermore][scattermore::geom_scattermore]; 
##' @export
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
##' p1 <- sc_dim(sce, reduction = 'UMAP', mapping = aes(colour = Cell_Cycle))
##' p2 <- sc_dim(sce, reduction = 'UMAP')
##' f1 <- p1 + sc_dim_geom_label()
##' f2 <- p2 + 
##'       sc_dim_geom_label(
##'         geom = shadowtext::geom_shadowtext,
##'         color='black',
##'         bg.color='white'
##'       )
setGeneric('sc_dim', 
    function(object, 
       dims=c(1,2), 
       reduction=NULL, 
       cells=NULL, 
       slot = "data", 
       mapping = NULL, 
       geom = sc_geom_point,
       ...)
   standardGeneric('sc_dim')
)

#' @importFrom methods setMethod
#' @rdname sc-dim-methods
#' @aliases sc_dim,Seurat
#' @exportMethod sc_dim
setMethod("sc_dim", 
  'Seurat', 
  function(
   object, 
   dims=c(1,2), 
   reduction=NULL, 
   cells=NULL, 
   slot = "data", 
   mapping = NULL,
   geom = sc_geom_point,
   ...) 
{
    d <- get_dim_data(object = object, features = NULL,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot)
    mapping <- .check_aes_mapping(object, mapping, data = d, prefix = 'ident')

    p <- sc_dim_internal(d, mapping, geom = geom, ...)
    p <- .add_class(p, "ggsc")
    return(p)
})

#' @rdname sc-dim-methods
#' @aliases sc_dim,SingleCellExperiment
#' @exportMethod sc_dim
setMethod('sc_dim', 'SingleCellExperiment', 
    function(
       object, 
       dims = c(1, 2), 
       reduction = NULL, 
       cells = NULL, 
       slot = 'data', 
       mapping = NULL, 
       geom = sc_geom_point,
       ...)
  {
    d <- .extract_sce_data(object = object, features = NULL, dims = dims, 
                      reduction = reduction, cells = cells, slot = slot)

    mapping <- .check_aes_mapping(object, mapping, data = d, prefix = 'label')

    p <- sc_dim_internal(d, mapping, geom = geom, ...)
    p <- .add_class(p, "ggsc")
    return(p)
})

##' @importFrom tidydr theme_dr
sc_dim_internal <- function(data, mapping, geom = sc_geom_point, ...) {
    default_mapping <- .set_default_mapping(data)
    mapping <- modifyList(default_mapping, mapping)
    p <- ggplot(data, mapping)
    params <- list(...)
    layers <- do.call(geom, params)
    p <- p + layers + theme_dr()
    return(p)
} 


.set_default_mapping <- function(data){
    dims <- names(data)[seq_len(3)]
    x <- aes(x=!!rlang::sym(dims[2]), y=!!rlang::sym(dims[3]))
    return(x)
}
