##' @title sc_feature
##' @rdname sc-feature
##' @param object Seurat object
##' @param features selected features (i.e., genes)
##' @param dims selected dimensions (must be a two-length vector) that are used in visualization
##' @param reduction reduction method, default is NULL and will use the default setting store in the object
##' @param cells selected cells to plot (default is all cells)
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param mapping aesthetic mapping
##' @param ncol number of facet columns if 'length(features) > 1'
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
sc_feature <- function(object, features, 
                    dims=c(1,2), reduction=NULL, 
                    cells=NULL, slot = "data", mapping, ncol=3, ...) {
    d <- get_dim_data(object = object, features = features,
                    dims = dims, reduction = reduction, 
                    cells = cells, slot = slot)

    d2 <- tidyr::pivot_longer(d, 4:ncol(d), names_to = "features")
    d2$features <- factor(d2$features, features)

    default_mapping <- aes_string(color="value")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }               

    p <- sc_dim_internal(d2, mapping, ...) +
        scale_color_gradient(low='grey', high='blue')           
        #scale_color_gradient2(low='blue', mid='grey', high='red') + 

    if (length(features) == 1) {
        p <- p + ggtitle(features) +
            theme(plot.title=element_text(size=rel(1.5), face='bold')) 
    } else {
        p <- p + facet_wrap(~features, ncol=ncol) +
            theme_bw() 
            theme(axis.text = element_blank(), 
                axis.title = element_blank(), 
                axis.ticks = element_blank(), 
                panel.grid = element_blank()
            ) 
    }

    return(p)
}

## features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A")
## sc_feature(pbmc, features)


