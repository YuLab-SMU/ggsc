##' @title sc_dim
##' @rdname sc-dim-methods
##' @param object Seurat object
##' @param dims selected dimensions (must be a two-length vector) that 
##' are used in visualization
##' @param reduction reduction method, default is NULL and will use the 
##' default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
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
##' sce <- logNormCounts(sce)
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
setMethod('sc_dim', 'SingleCellExperiment', 
          function(object, dims = c(1, 2), reduction = NULL, 
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
##' @importFrom SummarizedExperiment assay colData assayNames
##' @importFrom cli cli_abort
.extract_sce_data <- function(object, features = NULL, dims = c(1, 2), 
                              reduction = NULL, cells = NULL, slot = 1, 
                              plot.pie = FALSE, density=FALSE, grid.n = 400, 
                              joint = FALSE, joint.fun = prod, sp.coords=NULL){
    if (!is.null(cells)){
        object <- object[, cells]
    }
    
    xx <- colData(object) |> data.frame()
    reduced.dat <- NULL
    if (!is.null(dims)){
        if (length(reducedDimNames(object)) == 0){
            cli::cli_abort(c("The {.cls {class(object)}} didn't contain the results of reduction."))
        }
        if (is.null(reduction)){
            reduction <- 1
        }
        reduced.dat <- reducedDims(object)[[reduction]][,dims] |> 
            as.data.frame(check.names = FALSE)
        xx <- merge(reduced.dat, xx, by = 0)
        rownames(xx) <- xx$Row.names
        xx$Row.names <- NULL
    }

    if (!is.null(features)){
        if (slot == 'data'){
            if ('logcounts' %in% assayNames(object)){
                slot <- 'logcounts'
            }else{
                slot <- 1
            }
        }
        
        tmp <- assay(object, slot)
        tmp <- tmp[features, ,drop=FALSE] 

        if (density && !is.null(reduced.dat) && !plot.pie){
          tmp <- .buildWkde(w = tmp, coords = reduced.dat, n = grid.n, 
                            joint = joint, joint.fun = joint.fun) 
        }else if (density && !is.null(sp.coords) && !plot.pie){
          tmp <- .buildWkde(w = tmp, coords = sp.coords, n = grid.n,
                            joint = joint, joint.fun = joint.fun)
        }else{
          tmp <- tmp |> 
                 as('matrix') |> 
                 t() |> 
                 as.data.frame(check.names=FALSE)
        }

        xx <- merge(xx, tmp, by = 0)
        rownames(xx) <- xx$Row.names
        xx$Row.names <- NULL
    }
    return(xx)
}

##' @importFrom tidydr theme_dr
sc_dim_internal <- function(data, mapping, ...) {
    dims <- names(data)[seq_len(2)]
    ggplot(data, aes(.data[[dims[1]]], .data[[dims[2]]])) + 
        sc_geom_point(mapping, ...) + 
        theme_dr()
} 

get_dim_data <- function(object, features = NULL, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", 
                    plot.pie=FALSE, density = FALSE, 
                    grid.n = 400, joint = FALSE, 
                    joint.fun = prod, sp.coords = NULL
                    ) {
    rlang::check_installed('SeuratObject', 'for the internal function `get_dim_data()`.')
    reduced.dat <- NULL
    
    if (is.null(cells)) {
        cells <- colnames(object)
    }
    #xx <- data.frame(ident=SeuratObject::Idents(object)[cells])
    xx <- cbind(data.frame(ident = SeuratObject::Idents(object)[cells]), object@meta.data[cells,,drop=FALSE])
    
    if (!is.null(dims)) {
        if (is.null(reduction)) {
            reduction <- SeuratObject::DefaultDimReduc(object)
        }
        dims <- paste0(SeuratObject::Key(object = object[[reduction]]), dims)
        reduced.dat <- as.data.frame(SeuratObject::Embeddings(object[[reduction]])[cells, dims])
    }

    if (!is.null(features)){
        tmp <- SeuratObject::FetchData(object, vars = features, cells = cells, slot = slot)
        if (density && !is.null(reduced.dat) && !plot.pie){
            tmp <- .buildWkde(t(tmp), reduced.dat, grid.n, joint, joint.fun)
            xx <- cbind(reduced.dat, xx, tmp)
        }else if(density && !is.null(sp.coords) && !plot.pie){
            tmp <- .buildWkde(t(tmp), sp.coords, grid.n, joint, joint.fun)
            xx <- cbind(xx, tmp)
        }else if (!is.null(reduced.dat) && !density){
            xx <- cbind(reduced.dat, xx, tmp)
        }else{
            xx <- cbind(xx, tmp)
        }
    }else{
        if (!is.null(reduced.dat)){
            xx <- cbind(reduced.dat, xx)
        }
    }

    return(xx)
}


