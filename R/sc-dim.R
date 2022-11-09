##' @title sc_dim
##' @rdname sc-dim
##' @param object Seurat object
##' @param dims selected dimensions (must be a two-length vector) that are used in visualization
##' @param reduction reduction method, default is NULL and will use the default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return dimension reduction plot
##' @seealso
##'  [geom_scattermore][scattermore::geom_scattermore]; 
##' @export
sc_dim <- function(object, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", mapping = NULL, ...) {
    d <- get_dim_data(object = object, features = NULL,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot)

    default_mapping <- aes_string(color="ident")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }               
    p <- sc_dim_internal(d, mapping, ...)
    return(p)
}

##' @importFrom tidydr theme_dr
sc_dim_internal <- function(data, mapping, ...) {
    dims <- names(data)[1:2]
    ggplot(data, aes_string(dims[1], dims[2])) + 
        sc_geom_point(mapping, ...) + 
        theme_dr()
} 

##' @importFrom SeuratObject FetchData
get_dim_data <- function(object, features = NULL, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data") {
    if (is.null(cells)) {
        cells <- colnames(object)
    }
    if (!is.null(dims)) {
        if (is.null(reduction)) {
            reduction <- SeuratObject::DefaultDimReduc(object)
        }
        dims <- paste0(SeuratObject::Key(object = object[[reduction]]), dims)
    }
    SeuratObject::FetchData(object, vars = c(dims, "ident", 
        features), cells = cells, slot = slot)
}


