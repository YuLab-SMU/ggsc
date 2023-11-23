##' @title sc_dot
##' @rdname sc-dot-methods
##' @param object Seurat or SingleCellExperiment object
##' @param features selected features
##' @param group.by grouping factor
##' @param split.by additional split factor
##' @param cols colors of the points
##' @param col.min minimum scaled averaged expression threshold
##' @param col.max maximum scaled averaged expression threshold
##' @param dot.min the threshold of fraction of for the the smallest dot
##' @param dot.scale Scaling factor for size of points
##' @param scale whether to scale the expression value (default to TRUE)
##' @param scale.by scale the size of the points by `size` or `radius`
##' @param scale.min lower limit of scaling
##' @param scale.max upper limit of scaling
##' @param slot slot to pull expression data from (e.g., 'count' or 'data')
##' @param .fun user defined function that will be applied to selected features (default is NULL and there is no data operation)
##' @param mapping aesthetic mapping
##' @param ... additional parameters pass to 'ggplot2::geom_point()'
##' @seealso
##'  [DotPlot][Seurat::DotPlot]; 
##' @return dot plot to visualize feature expression distribution
##' @importFrom utils modifyList
##' @importFrom ggplot2 aes_string
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_point
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
##' set.seed(123)
##' genes <- rownames(sce) |> sample(6) 
##' sc_dot(sce, genes[1:5], 'Treatment', slot = 'logcounts')
##' 
setGeneric('sc_dot', function(object, features, group.by=NULL, split.by = NULL,
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
                    group.by=NULL, split.by = NULL, cols = c("lightgrey", "blue"),
                    col.min = -2.5, col.max = 2.5,
                    dot.min = 0, dot.scale = 6,
                    slot = "data", .fun = NULL,
                    mapping = NULL,
                    scale = TRUE, scale.by = 'radius',
                    scale.min = NA, scale.max = NA,
                    ...) {
    d <- get_dim_data(object, dims=NULL, features=features, slot=slot)
    d <- tidyr::pivot_longer(d, 2:ncol(d), names_to = "features")
    d$features <- factor(d$features, levels = features)
    if (!is.null(.fun)) {
        d <- .fun(d)
    }
    if (is.null(group.by)) {
    	group.by <- "ident"
    }
    return(.ReturnDotPlot(d, features, group.by, split.by, cols,
	col.min, col.max, dot.min, dot.scale, mapping, scale, scale.by,
	scale.min, scale.max, ...))
})

##' @rdname sc-dot-methods
##' @aliases sc_dot,SingleCellExperiment
##' @exportMethod sc_dot
setMethod('sc_dot', 'SingleCellExperiment', 
          function(
             object, features, group.by=NULL, split.by = NULL,
             cols = c("lightgrey", "blue"),
             col.min=-2.5, col.max=2.5, dot.min=0, dot.scale=6,
             slot = 'data', .fun = NULL, mapping = NULL,
             scale = TRUE, scale.by = 'radius',
             scale.min = NA, scale.max = NA,
             ...){
    d <- .extract_sce_data(object, dims = NULL, features = features, slot = slot)
    d <- tidyr::pivot_longer(d, seq(ncol(d) - length(features) + 1, ncol(d)), names_to = "features")
    if (is.numeric(features)){
        features <- rownames(object)[features]
    }
    d$features <- factor(d$features, levels = features)
    if (!is.null(.fun)) {
        d <- .fun(d)
    }
    if (is.null(group.by)) {
    	group.by <- "label"
    }
    return(.ReturnDotPlot(d, features, group.by, split.by, cols,
	col.min, col.max, dot.min, dot.scale, mapping, scale, scale.by,
	scale.min, scale.max, ...))
})


.ReturnDotPlot <- function(d, features, group.by, split.by, cols,
	col.min, col.max, dot.min, dot.scale, mapping, scale, scale.by,
	scale.min, scale.max, ...) {
    #Some parts in the function is adapted from Seurat::DotPlot
    #Currently, feature.groups and cluster.idents are not implemented
    feature.groups <- NULL
    split.colors <- !is.null(split.by) && !any(cols %in% rownames(RColorBrewer::brewer.pal.info))
	scale.func <- switch(
	    EXPR = scale.by,
	    'size' = scale_size,
	    'radius' = scale_radius,
	    stop("'scale.by' must be either 'size' or 'radius'")
	)
	
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
        dplyr::summarise(avg.exp=mean(expm1(value)),
        	      pct.exp=.PercentAbove(value, 0))
    ngroup <- length(id.levels)
    if (ngroup == 1) {
        scale <- FALSE
        warning(
            "Only one identity present, the expression values will be not scaled",
            call. = FALSE,
            immediate. = TRUE
        )
    } else if (ngroup < 5 & scale) {
        warning(
            "Scaling data with a low number of groups may produce misleading results",
            call. = FALSE,
            immediate. = TRUE
        )
    }
    
    .scale.fun <- function(x) {
    	if (scale) {
    		scaled <- scale(log1p(x))
    		scaled <- .MinMax(scaled, min=col.min, max=col.max)
    	} else {
    		scaled <- log1p(x)
    	}
    	return(scaled)
    }
    
    avg.exp <- avg.exp %>% 
        dplyr::group_by(features) %>%
        dplyr::mutate(avg.exp.scaled=.scale.fun(avg.exp))

    if (split.colors) {
        avg.exp <- avg.exp %>% 
            dplyr::mutate(avg.exp.scaled=as.numeric(cut(avg.exp.scaled, breaks = 20)))
    }
    avg.exp$pct.exp[avg.exp$pct.exp < dot.min] <- NA
    avg.exp$pct.exp <- avg.exp$pct.exp * 100

    if (split.colors) {
        splits.use <- unlist(x = lapply(
            X = avg.exp[[group.by]],
            FUN = function(x)
                sub(
                    paste0(".*_(",
                        paste(sort(unique(x = splits), decreasing = TRUE),
                            collapse = '|'
                        ),")$"),
                    "\\1", x
                )
            )
        )
        avg.exp$colors <- mapply(
            FUN = function(color, value) {
                return(colorRampPalette(colors = c('grey', color))(20)[value])
            },
            color = cols[splits.use],
            value = avg.exp$avg.exp.scaled
        )
    }
    color.by <- ifelse(test = split.colors, yes = 'colors', no = 'avg.exp.scaled')
    if (!is.na(x = scale.min)) {
        avg.exp[avg.exp$pct.exp < scale.min, 'pct.exp'] <- scale.min
    }
    if (!is.na(x = scale.max)) {
        avg.exp[avg.exp$pct.exp > scale.max, 'pct.exp'] <- scale.max
    }
    
    default_mapping <- aes_string(color=color.by, size="pct.exp")
    if (is.null(mapping)) {
        mapping <- default_mapping
    } else {
        mapping <- modifyList(default_mapping, mapping)
    }
	p <- ggplot(avg.exp, aes(x=features, y=.data[[group.by]])) +
    	geom_point(mapping, ...)+
    	scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max))+
        theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
        guides(size = guide_legend(title = 'Percent Expressed')) +
        labs(
            x = 'Features',
            y = ifelse(test = is.null(x = split.by), yes = 'Identity', no = 'Split Identity')
        )+
	    theme_minimal()
    if (split.colors) {
        p <- p + scale_color_identity()
    } else if (length(x = cols) == 1) {
        p <- p + scale_color_distiller(palette = cols)
    } else {
        p <- p + scale_color_gradient(low = cols[1], high = cols[2])
    }
	if (!split.colors) {
        p <- p + guides(color = guide_colorbar(title = 'Average Expression'))
    }
    return(p)
}

.PercentAbove <- function(x, threshold) {
    return(length(x = x[x > threshold]) / length(x = x))
}

.MinMax <- function(data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}