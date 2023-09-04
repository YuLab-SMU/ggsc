##' @title sc_dim
##' @rdname sc-dim-methods
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
setGeneric('sc_dim', 
           function(object, dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", mapping = NULL, 
                    ...)
               standardGeneric('sc_dim')
)

#' @importFrom methods setMethod
#' @rdname sc-dim-methods
#' @aliases sc_dim,Seurat
#' @exportMethod sc_dim
setMethod("sc_dim", 'Seurat', function(object, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", mapping = NULL, ...) {
    d <- get_dim_data(object = object, features = NULL,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot)

    default_mapping <- aes(color=.data$ident)
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }               
    p <- sc_dim_internal(d, mapping, ...)
    return(p)
})

#' @rdname sc-dim-methods
#' @aliases sc_dim,SingleCellExperiment
#' @exportMethod sc_dim
setMethod('sc_dim', 'SingleCellExperiment', function(object, dims = c(1, 2), reduction = NULL, 
                                                     cells = NULL, slot = 'data', mapping = NULL, ...){
    d <- .extract_sce_data(object = object, features = NULL, dims = dims, 
                      reduction = reduction, cells = cells, slot = slot)
    
    default_mapping <- aes(color = .data$label)
    if (is.null(mapping)){
        mapping <- default_mapping
    }else{
        mapping <- modifyList(default_mapping, mapping)
    }
    p <- sc_dim_internal(d, mapping, ...)
    return(p)

})

##' @importFrom methods as
##' @importFrom SingleCellExperiment reducedDims reducedDimNames
##' @importFrom SummarizedExperiment assay colData
##' @importFrom cli cli_abort
.extract_sce_data <- function(object, features = NULL, dims = c(1, 2), reduction = NULL, cells = NULL, slot = 1){
    if (!is.null(cells)){
        object <- object[, cells]
    }
    
    xx <- colData(object) |> as.data.frame(check.names = FALSE) |> suppressWarnings()

    if (!is.null(dims)){
        if (length(reducedDimNames(object)) == 0){
            cli::cli_abort(c("The {.cls {class(object)}} didn't contain the results of reduction."))
        }
        if (is.null(reduction)){
            reduction <- 1
        }
        tmp.reduced <- reducedDims(object)[[reduction]][,dims] |> as.data.frame(check.names = FALSE)
        xx <- merge(xx, tmp.reduced, by = 0)
        rownames(xx) <- xx$Row.names
        xx$Row.names <- NULL
    }

    if (!is.null(features)){
        if (slot == 'data'){
            slot <- 1
        }
        tmp <- assay(object, slot)
        tmp <- tmp[features, ,drop=FALSE] |> 
               as('matrix') |> 
               t() |> 
               as.data.frame(check.names=FALSE)

        xx <- merge(xx, tmp, by = 0)
        rownames(xx) <- xx$Row.names
        xx$Row.names <- NULL
    }
    return(xx)
}

##' @importFrom tidydr theme_dr
sc_dim_internal <- function(data, mapping, ...) {
    dims <- names(data)[1:2]
    ggplot(data, aes(.data[[dims[1]]], .data[[dims[2]]])) + 
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


