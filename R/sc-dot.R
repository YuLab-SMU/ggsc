##' @title sc_dot
##' @rdname sc-dot-methods
##' @param object Seurat object
##' @param features selected features
##' @param group.by grouping factor
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param .fun user defined function that will be applied to selected features (default is NULL and there is no data operation)
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'ggplot2::geom_dot()'
##' @return dot plot to visualize feature expression distribution
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
setGeneric('sc_dot', function(object, features, group.by, split.by = NULL,
	                             cols = c("lightgrey", "blue"),
	                             col.min = -2.5, col.max = 2.5,
	                             dot.min = 0, dot.scale = 6,
                                 slot = "data", .fun = NULL, mapping = NULL,
                                 scale = TRUE, scale.by = 'radius',
                                 scale.min = NA, scale.max = NA,
                                 ...)
    standardGeneric('sc_dot')
)

##' @rdname sc-dot-methods
##' @aliases sc_dot,Seurat
##' @exportMethod sc_dot
setMethod("sc_dot", 'Seurat', function(object, features, 
                    group.by, split.by = NULL, cols = c("lightgrey", "blue"),
                    col.min = -2.5, col.max = 2.5,
                    dot.min = 0, dot.scale = 6,
                    slot = "data", .fun = NULL,
                    mapping = NULL,
                    scale = TRUE, scale.by = 'radius',
                    scale.min = NA, scale.max = NA,
                    ...) {
    d <- get_dim_data(object, dims=NULL, features=features)
    d <- tidyr::pivot_longer(d, 2:ncol(d), names_to = "features")
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

##' @rdname sc-dot-methods
##' @aliases sc_dot,SingleCellExperiment
##' @exportMethod sc_dot
setMethod('sc_dot', 'SingleCellExperiment', 
          function(
             object, features, group.by, split.by = NULL,
             cols = c("lightgrey", "blue"),
             col.min=-2.5, col.max=2.5, dot.min=0, dot.scale=6,
             slot = 'data', .fun = NULL, mapping = NULL, ncol = 3,
             scale = TRUE, scale.by = 'radius',
             scale.min = NA, scale.max = NA,
             ...){
    #From Seurat::DotPlot
    split.colors <- !is.null(split.by) && !any(cols %in% rownames(RColorBrewer::brewer.pal.info))
	scale.func <- switch(
	    EXPR = scale.by,
	    'size' = scale_size,
	    'radius' = scale_radius,
	    stop("'scale.by' must be either 'size' or 'radius'")
	)
    d <- .extract_sce_data(object, dims = NULL, features = features)
    d <- tidyr::pivot_longer(d, seq(ncol(d) - length(features) + 1, ncol(d)), names_to = "features")
    if (is.numeric(features)){
        features <- rownames(object)[features]
    }
    d$features <- factor(d$features, levels = features)
    if (!is.null(.fun)) {
        d <- .fun(d)
    }
    id.levels <- levels(d[[group.by]])
    if (!is.null(split.by)) {
    	splits <- d[[split.by]]
        if (split.colors) {
            if (length(unique(splits)) > length(cols)) {
                stop(paste0("Need to specify at least ", length(unique(splits)), " colors using the cols parameter"))
            }
	        cols <- cols[1:length(unique(splits))]
	        names(cols) <- unique(splits)
        }
        d[[group.by]] <- paste(d[[group.by]], splits, sep = '_')
        unique.splits <- unique(splits)
        id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)),
        	"_", rep(x = unique(x = splits), times = length(x = id.levels)))
    }
    avg.exp <- d %>% 
        dplyr::group_by(.data[[group.by]], features) %>%
        dplyr::summarise(avg.exp=mean(value),
        	      pct.exp=.PercentAbove(value, 0))

    default_mapping <- aes_string(color="avg.exp", size="pct.exp")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
	p <- ggplot(avg.exp, aes(x=features, y=.data[[group.by]])) +
    	geom_point(mapping)+
	    scale_color_gradient(low=cols[1],high=cols[2])+
	    theme_minimal()
    return(p)    

})


.PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
}