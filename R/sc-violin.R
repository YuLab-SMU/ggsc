##' @title sc_violin
##' @rdname sc-violin-methods
##' @param object Seurat object
##' @param features selected features
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param .fun user defined function that will be applied to selected features (default is NULL and there is no data operation)
##' @param mapping aesthetic mapping
##' @param ncol number of facet columns if 'length(features) > 1'
##' @param ... additional parameters pass to 'ggplot2::geom_geom_violin()'
##' @return violin plot to visualize feature expression distribution
##' @seealso
##'  [geom_violin][ggplot2::geom_violin]; 
##' @importFrom utils modifyList
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_violin
##' @importFrom ggplot2 facet_wrap
##' @importFrom tidyr pivot_longer
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
##' set.seed(123)
##' genes <- rownames(sce) |> sample(6) 
##' sc_violin(sce, genes[1], slot = 'logcounts')
##' sc_violin(sce, genes[1], slot = 'logcounts',
##'      .fun=function(d) dplyr::filter(d, value > 0)
##'      ) +
##'      ggforce::geom_sina(size=.1)
##' sc_violin(sce, genes, slot = 'logcounts') +
##'   theme(axis.text.x = element_text(angle=45, hjust=1))
setGeneric('sc_violin', function(object, features, cells=NULL, 
                                 slot = "data", .fun = NULL, mapping = NULL, ncol=3, ...)
    standardGeneric('sc_violin')
)

##' @rdname sc-violin-methods
##' @aliases sc_violin,Seurat
##' @exportMethod sc_violin
setMethod("sc_violin", 'Seurat', function(object, features, 
                    cells=NULL, slot = "data", .fun = NULL, 
                    mapping = NULL, ncol=3, ...) {
    d <- get_dim_data(object, dims=NULL, features=features)
    indx.f <- seq(ncol(d) - length(features) + 1, ncol(d))
    d <- tidyr::pivot_longer(d, cols=indx.f, names_to = "features") 
    d$features <- factor(d$features, levels = features)
    if (!is.null(.fun)) {
        d <- .fun(d)
    }
    default_mapping <- aes_string(fill="ident")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
    p <- ggplot(d, aes_string("ident", "value")) + 
        geom_violin(mapping, ...) #+ 
        #ggforce::geom_sina(size=.1)
    
    if (length(features) > 1) {
        p <- p + facet_wrap(~features, ncol=ncol)
    }
    return(p)
})

##' @rdname sc-violin-methods
##' @aliases sc_violin,SingleCellExperiment
##' @exportMethod sc_violin
setMethod('sc_violin', 'SingleCellExperiment', 
          function(
             object, features, cells = NULL, slot = 'data', 
             .fun = NULL, mapping = NULL, ncol = 3, ...){

    d <- .extract_sce_data(object, dims = NULL, features = features)
    
    d <- tidyr::pivot_longer(d, seq(ncol(d) - length(features) + 1, ncol(d)), names_to = "features")

    if (is.numeric(features)){
        features <- rownames(object)[features]
    }

    d$features <- factor(d$features, levels = features)
    if (!is.null(.fun)) {
        d <- .fun(d)
    }
    default_mapping <- aes_string(fill="label")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
    p <- ggplot(d, aes_string("label", "value")) +
        geom_violin(mapping, ...) #+
        #ggforce::geom_sina(size=.1)

    if (length(features) > 1) {
        p <- p + facet_wrap(~features, ncol=ncol)
    }
    return(p)    

})


