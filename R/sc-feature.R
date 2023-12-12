##' @title sc_feature
##' @rdname sc-feature-methods
##' @param object Seurat object
##' @param features selected features (i.e., genes)
##' @param dims selected dimensions (must be a two-length vector) 
##' that are used in visualization
##' @param reduction reduction method, default is NULL and will 
##' use the default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping
##' @param ncol number of facet columns if 'length(features) > 1'
##' @param density whether plot the 2D weighted kernel density, default is FALSE.
##' @param grid.n number of grid points in the two directions to estimate 2D 
##' weighted kernel density, default is 100. 
##' @param joint whether joint the multiple features with \code{joint.fun}, 
##' default is FALSE.
##' @param joint.fun how to joint the multiple features if \code{joint=TRUE},
##' default is prod.
##' @param common.legend whether to use \code{facet_wrap} to display the multiple 
##' \code{features}, default is TRUE.
##' @param ... additional parameters pass to 'scattermore::geom_scattermore()'
##' @return dimension reduction plot colored by selected features
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 theme_bw
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_blank
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 rel
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
##' sce <- runTSNE(sce, assay.type = 'logcounts')
##' set.seed(123)
##' genes <- rownames(sce) |> sample(6)
##' p1 <- sc_feature(sce, genes[1], slot='logcounts', reduction = 'TSNE')
##' p2 <- sc_feature(sce, genes, slot='logcounts', reduction = 'TSNE')
##' f1 <- sc_dim(sce, slot='logcounts', reduction = 'TSNE') +
##'       sc_dim_geom_feature(sce, genes[1], color='black')
##' f2 <- sc_dim(sce, alpha=.3, slot='logcounts', reduction = 'TSNE') +
##'     ggnewscale::new_scale_color() +
##'     sc_dim_geom_feature(sce, genes, mapping=aes(color=features)) +
##'     scale_color_viridis_d()
##' p1 + p2 + f1 + f2
setGeneric('sc_feature', function(object, 
                                  features, 
                                  dims = c(1, 2), 
                                  reduction = NULL, 
                                  cells = NULL, 
                                  slot = 'data', 
                                  mapping = NULL, 
                                  ncol = 3,
                                  density = FALSE,
                                  grid.n = 100, 
                                  joint = FALSE,
                                  joint.fun = prod,
                                  common.legend = TRUE,
                                  ...)
    standardGeneric('sc_feature')
)

##' @rdname sc-feature-methods
##' @aliases sc_feature,Seurat
##' @exportMethod sc_feature
setMethod('sc_feature', 'Seurat', function(object, features, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", mapping=NULL, 
                    ncol=3, density = FALSE, grid.n = 100, joint = FALSE,
                    joint.fun = prod, common.legend = TRUE, ...) {
    d <- get_dim_data(object = object, features = features,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot, density = density, 
                    grid.n = grid.n, joint = joint, joint.fun = joint.fun)

    nm.f <- length(features)
    if (density){
        valnm <- 'density'
        if(joint){
           valnm <- "joint density"
           nm.f <- nm.f + 1
        }
    }else{
        valnm <- slot
    }

    indx.f <- seq(ncol(d)-nm.f + 1, ncol(d))

    features <- colnames(d)[indx.f]

    d2 <- tidyr::pivot_longer(d, indx.f, names_to = "features", values_to = valnm)
    d2$features <- factor(d2$features, features)

    default_mapping <- aes_string(color=valnm)
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }

    p <- sc_dim_internal(d2, mapping, ...) +
        scale_color_gradient(low='grey', high='blue')           
        #scale_color_gradient2(low='blue', mid='grey', high='red') + 

    p <- p + .feature_setting(features=features, ncol=ncol)
    
    if (!common.legend && length(features) > 1){
        p <- .split.by.feature(p, ncol)
    }
    return(p)
})

##' @rdname sc-feature-methods
##' @aliases sc_feature,SingleCellExperiment
##' @exportMethod sc_feature
setMethod("sc_feature", "SingleCellExperiment", 
          function(object, features, dims = c(1, 2), reduction = NULL, 
                   cells = NULL, slot = 'data', mapping = NULL, ncol = 3, 
                   density = FALSE, grid.n = 100, joint = FALSE, 
                   joint.fun = prod, common.legend = TRUE, ...){
    if (slot == 'data'){
        if ('logcounts' %in% assayNames(object)){
            slot <- 'logcounts'
        }else{
            slot <- 1
        }
    }
              
    d <- .extract_sce_data(object = object, features = features, dims = dims, 
                           reduction = reduction, cells = cells, slot = slot,
                           density = density, grid.n = grid.n, joint = joint,
                           joint.fun = joint.fun
         )

    nm.f <- length(features)

    if (density){
       valnm <- 'density'
       if (joint) {
           valnm <- "joint"
           nm.f <- nm.f + 1
       }
    }else{
       if (is.numeric(slot)){
           slot <- assayNames(object)[slot]
       }
       valnm <- slot
    }

    indx.f <- seq(ncol(d) - nm.f + 1, ncol(d))
    
    d2 <- tidyr::pivot_longer(
            d, 
            indx.f, 
            names_to = 'features',
            values_to = valnm
          )

    features <- colnames(d)[indx.f]

    d2$features <- factor(d2$features, features)
    
    default_mapping <- aes_string(color=valnm)
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }

    p <- sc_dim_internal(d2, mapping, ...) +
        scale_color_gradient(low='grey', high='blue')
        #scale_color_gradient2(low='blue', mid='grey', high='red') +

    p <- p + .feature_setting(features=features, ncol=ncol)

    if (!common.legend && length(features) > 1){
        p <- .split.by.feature(p, ncol)
    }
    return(p)    
})


.feature_setting <- function(features, ncol, plot.pie=FALSE) {
    if (length(features) == 1) {
        res <- list(ggtitle(features),
            theme(plot.title=element_text(size=rel(1.5), face='bold'))
        ) 
    }else if(missing(features) || is.null(features)){
        res <- theme_bw2()
    }else if(!plot.pie) {
        res <- list(facet_wrap(~features, ncol=ncol),
            theme_bw2()
        ) 
    }else if(plot.pie){
        res <- theme_bw2(legend.title=element_text(size = 16,face="bold"),
                         legend.text=element_text(size = 15),
                         legend.key.size = grid::unit(0.45, 'cm'))
    }
    return(res)
}

##' @importFrom ggplot2 %+replace%
theme_bw2 <- function(...) {
    theme_bw() %+replace% 
    theme(axis.text = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        panel.grid = element_blank()
      ) %+replace% 
    theme(...)
}

