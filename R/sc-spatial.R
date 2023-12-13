##' @title sc_spatial
##' @rdname sc-spatial-methods
##' @param object Seurat object
##' @param features selected features to be visualized
##' @param sample.id the index name of sample id, which only 
##' work with SingleCellExperiment or SpatialExperiment.
##' @param image.id the index name of image id, which only work 
##' with SingleCellExperiment or SpatialExperiment.
##' @param slot if plotting a feature, which data will be used 
##' (e.g., 'data', 'counts'), the assay name if object
##' is SingleCellExperiment or SpatialExperiment.
##' @param plot.pie logical whether plot the features with pie, default is \code{FALSE}.
##' @param pie.radius.scale numeric scale to the radius of pie only work with \code{plot.pie=TRUE},
##' default is 0.3.
##' @param image.plot whether to display the issue image as background.
##' @param image.first.operation character which the first operation to 
##' image, 'rotate' or 'mirror', default is 'rotate'.
##' @param image.rotate.degree integer the degree to ratate image, default is NULL.
##' @param image.mirror.axis character the direction to mirror the image, default is 'h'.
##' @param remove.point whether to remove the spot points, it is nice 
##' if your just view the issue image, default is FALSE.
##' @param mapping aesthetic mapping, default is NULL.
##' @param ncol integer number of facet columns if 'length(features) > 1', default is 6.
##' @param density whether plot the 2D weighted kernel density, default is FALSE.
##' @param grid.n number of grid points in the two directions to estimate 2D
##' weighted kernel density, default is 100.
##' @param joint whether joint the multiple features with \code{joint.fun},
##' default is FALSE.
##' @param joint.fun how to joint the multiple features if \code{joint = TRUE},
##' default is prod.
##' @param common.legend whether to use \code{facet_wrap} to display the multiple
##' \code{features}, default is TRUE.
##' @param point.size the size of point, default is 5.
##' @param ... additional parameters.
##' @return ggplot object
##' @importFrom grid rasterGrob unit
##' @importFrom ggplot2 facet_grid annotation_custom
##' @importFrom ggplot2 xlab ylab geom_blank coord_fixed 
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom Seurat DefaultAssay
##' @export
##' @examples
##' \dontrun{
##' library(STexampleData)
##' # create ExperimentHub instance
##' eh <- ExperimentHub()
##' # query STexampleData datasets
##' myfiles <- query(eh, "STexampleData")
##' spe <- myfiles[["EH7538"]]
##' spe <- spe[, colData(spe)$in_tissue == 1]
##' set.seed(123)
##' genes <- rownames(spe) |> sample(6) 
##' p <- sc_spatial(spe, features = genes, 
##'                 image.rotate.degree = -90, 
##'                 image.mirror.axis = NULL, 
##'                 ncol = 3)
##' }
setGeneric('sc_spatial', function(object, 
                                  features = NULL, 
                                  sample.id = NULL, 
                                  image.id = NULL, 
                                  slot = "data", 
                                  plot.pie = FALSE, 
                                  pie.radius.scale = 0.3,
                                  image.plot = TRUE, 
                                  image.first.operation = 'rotate',
                                  image.rotate.degree = NULL,
                                  image.mirror.axis = NULL,
                                  remove.point = FALSE,
                                  mapping = NULL,
                                  ncol = 6,
                                  density = FALSE,
                                  grid.n = 100,
                                  joint = FALSE,
                                  joint.fun = prod,
                                  common.legend = TRUE,
                                  point.size = 5,
                                  ...) 
           standardGeneric('sc_spatial')
)

##' @importFrom grDevices as.raster
##' @rdname sc-spatial-methods
##' @aliases sc_spatial,Seurat
##' @exportMethod sc_spatial
setMethod("sc_spatial", 'Seurat', 
          function(object, features = NULL, slot = "data", 
                   plot.pie = FALSE, pie.radius.scale = .3,
                   image.plot = TRUE, image.first.operation = 'rotate', 
                   image.rotate.degree = NULL, image.mirror.axis = 'v', 
                   remove.point = FALSE, mapping = NULL, ncol = 6, 
                   density=FALSE, grid.n = 100, joint = FALSE, 
                   joint.fun = prod, common.legend = TRUE, point.size = 5, ...) {
    images <- SeuratObject::Images(object = object, 
                    assay = Seurat::DefaultAssay(object = object)
                )
    img <- object@images[[images]]@image 
    if (!is.null(img)) img <- as.raster(img)
    
    coords.da <- SeuratObject::GetTissueCoordinates(object = object[[images]])
    
    d <- get_dim_data(object = object, features = features, dims = NULL, 
                      density = density, grid.n = grid.n, joint = joint,
                      joint.fun = joint.fun, sp.coords=coords.da)

    nm.f <- length(features)

    if (density){
       valnm <- 'density'
       if (joint){
           #valnm <- "joint_density"
           nm.f <- nm.f + 1
       }
    }else{
       valnm <- slot
    }
    d <- cbind(coords.da, d)

    default_mapping <- aes_string(x = colnames(coords.da)[2], y = colnames(coords.da)[1])

    if (!is.null(features)){

        indx.f <- seq(ncol(d) - nm.f + 1, ncol(d))
        features <- colnames(d)[indx.f]
        if (plot.pie){
            d <- d[rowSums(d[,features,drop=FALSE]) != 0,,drop=FALSE]
        }
        d <- tidyr::pivot_longer(
               d, 
               indx.f, 
               names_to = 'features',
               values_to = valnm
             )
        d$features <- factor(d$features, levels=features)
        if (!plot.pie){
           default_mapping <- modifyList(default_mapping, aes_string(color = valnm))
        }else{
           colnames(d)[colnames(d) == valnm] <- 'value'
        }
    }

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    ratio <- .cal_ratio(d, mapping)

    if (!plot.pie){
        p <- ggplot(d, mapping)
    }else{
        p <- ggplot()
    }    

    if (image.plot && !is.null(img)){
        img.annot <- .build_img_annot_layer(img, 
                                            image.first.operation, 
                                            image.rotate.degree, 
                                            image.mirror.axis)
        p <- p + img.annot
    }

    if ((!remove.point && (!is.null(features) || (any(names(mapping) %in% c('color', 'colour')) && is.null(features))) && !plot.pie)){
        p <- p + sc_geom_point(pointsize = point.size, ...)
    }else if (!remove.point && plot.pie){
        rlang::check_installed('scatterpie', 'is required when `plot.pie=TRUE`')
        p <- p + scatterpie::geom_scatterpie(data=d, mapping=mapping, cols='features', long_format=TRUE, pie_scale=pie.radius.scale, ...)
    }else{
        p <- p + geom_blank()
    }

    p <- p +
         .feature_setting(features, ncol, plot.pie) +
         ylab(NULL) +
         xlab(NULL) +
         coord_fixed(ratio=ratio) +
         theme_bw2() 

    color.aes <- .check_aes_exits(p$mapping, c('color', 'colour'))
    if (!is.null(color.aes)) {
        type.color.value <- p$data |> dplyr::pull(!!color.aes)
        if (inherits(type.color.value, 'numeric')) {
            p <- p + scale_color_gradientn(colours = SpatialColors(n=100))
        }
    }

    if (plot.pie){
        Type.cols <- .set_default_cols(length(features))
        p <- p + scale_fill_manual(values=Type.cols, name='Type')
    }    

    if (!common.legend && length(features) > 1 && !plot.pie){
        ncol <- min(length(features), ncol)
        p <- .split.by.feature(p, ncol, joint)
    }    
    return(p)    
})

#' @importFrom SingleCellExperiment int_metadata
#' @importFrom ggplot2 scale_fill_manual
#' @rdname sc-spatial-methods
#' @aliases sc_spatial,SingleCellExperiment
#' @exportMethod sc_spatial
setMethod('sc_spatial', 'SingleCellExperiment', function(object, 
                                                         features = NULL, 
                                                         sample.id = NULL, 
                                                         image.id = NULL, 
                                                         slot = 1,
                                                         plot.pie = FALSE,
                                                         pie.radius.scale = .3,
                                                         image.plot = TRUE,
                                                         image.first.operation = 'rotate',
                                                         image.rotate.degree = NULL,
                                                         image.mirror.axis = 'v',
                                                         remove.point = FALSE,
                                                         mapping = NULL,
                                                         ncol = 6,
                                                         density = FALSE,
                                                         grid.n = 100,
                                                         joint = FALSE,
                                                         joint.fun = prod,
                                                         common.legend = TRUE,
                                                         ...
                                                        ){
    if (!"imgData" %in% names(int_metadata(object))){
        cli::cli_abort(c("The {.cls {class(object)}} didn't have the image data."))
    }
    
    img.da <- .extract_img_data(object, sample.id = sample.id, image.id = image.id)
    
    coords.da <- .extract_coords(object, img.da)

    if (is.numeric(features)){
        features <- rownames(object)[features]
    }

    features.da <- .extract_sce_data(object, features = features, dims = NULL, 
                                     cells = NULL, slot = slot, plot.pie = plot.pie, density=density, 
                                     grid.n = grid.n, joint = joint, joint.fun = joint.fun, 
                                     sp.coords = coords.da)

    d <- merge(coords.da, features.da, by = 0)
    rownames(d) <- d$Row.names
    d$Row.names <- NULL

    default_mapping <- aes_string(x = colnames(coords.da)[2], y = colnames(coords.da)[1])
    if (!is.null(features)){
        if (plot.pie){
            d <- d[rowSums(d[,features,drop=FALSE]) != 0,,drop=FALSE]
        }
        nm.f <- length(features)
        if (density){
           valnm <- 'density'
           if (joint){
               #valnm <- "joint_density"
               nm.f <- nm.f + 1
           }
        }else{
           if (is.numeric(slot)){
               slot <- assayNames(object)[slot]
           }
           valnm <- slot
        }

        indx.f <- seq(ncol(d)- nm.f + 1, ncol(d))
        features <- colnames(d)[indx.f]
        
        d <- tidyr::pivot_longer(d, indx.f, 
                                 names_to = 'features', values_to = valnm)
        d$features <- factor(d$features, levels=features)
        if (!plot.pie){
           default_mapping <- modifyList(default_mapping, aes_string(color = valnm))
        }else{
           colnames(d)[colnames(d) == valnm] <- 'value' 
        }
    }

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    ratio <- .cal_ratio(d, mapping)

    if (!plot.pie){
        p <- ggplot(d, mapping)
    }else{
        p <- ggplot()
    }

    if (image.plot && !is.null(img.da)){
        img.annot <- .build_img_annot_layer(img.da, 
                                            image.first.operation, 
                                            image.rotate.degree, 
                                            image.mirror.axis)
        p <- p + img.annot
    }

    if ((!remove.point && (!is.null(features) || (any(names(mapping) %in% c('color', 'colour')) && is.null(features))) && !plot.pie)){
        p <- p + sc_geom_point(pointsize = point.size, ...) 
    }else if (!remove.point && plot.pie){
        rlang::check_installed('scatterpie', 'is required when `plot.pie=TRUE`')
        p <- p + scatterpie::geom_scatterpie(data=d, mapping=mapping, cols='features', long_format=TRUE, pie_scale = pie.radius.scale, ...)
    }else{
        p <- p + geom_blank()
    }
    
    p <- p +
         .feature_setting(features, ncol, plot.pie) +
         ylab(NULL) +
         xlab(NULL) +
         coord_fixed(ratio = ratio) +
         theme_bw2() 

    color.aes <- .check_aes_exits(p$mapping, c('color', 'colour'))
    if (!is.null(color.aes)) {
        type.color.value <- p$data |> dplyr::pull(!!color.aes)
        if (inherits(type.color.value, 'numeric')) {
            p <- p + scale_color_gradientn(colours = SpatialColors(n=100))
        }
    }
    if (plot.pie){
        Type.cols <- .set_default_cols(length(features)) 
        p <- p + scale_fill_manual(values=Type.cols, name='Type')
    }

    if (!common.legend && length(features) > 1 && !plot.pie){
        ncol <- min(length(features), ncol)
        p <- .split.by.feature(p, ncol, joint)
    }
    return(p)
})

#' @importFrom SingleCellExperiment int_metadata int_colData
.extract_img_data <- function(x, sample.id = NULL, image.id = NULL){
    img.da <- int_metadata(x)[['imgData']]
    if (nrow(img.da)==0){
        return(NULL)
    }    
    if (is.null(sample.id)){
        sample.id <- unique(img.da$sample_id)[1]
    }
    img.da <- img.da[img.da$sample_id == sample.id, ]

    if (is.null(image.id)){
        img.da <- img.da[1,]
    }else{
        img.da <- img.da[img.da$image_id == image.id, ]
    }
    return(img.da)
}

.extract_coords <- function(x, image.da){
    if (is.null(image.da)){
        scaleFactor <- 1.0
    }else{
        scaleFactor <- image.da$scaleFactor
    }
    x <- int_colData(x)
    x <- x[['spatialCoords']] * scaleFactor
    return(x)
}

.build_img_annot_layer <- function(image.da, 
                                   image.first.operation = NULL, 
                                   image.rotate.degree = NULL, 
                                   image.mirror.axis = NULL){
    if (!is.null(image.first.operation)){
        image.first.operation <- match.arg(image.first.operation, c('rotate', 'mirror'))
    }else{
        image.first.operation <- 'rotate'
    }

    if (inherits(image.da, 'raster')){
        img <- image.da
    }else{
        img <- image.da[['data']][[1]] |> as.raster()
    }

    if (image.first.operation == 'rotate'){
        if (!is.null(image.rotate.degree)){
            img <- .rotate.image(img, image.rotate.degree)
        } 
        if (!is.null(image.mirror.axis)){
            img <- .mirror.image(img, image.mirror.axis)
        } 
    }else{
        if (!is.null(image.mirror.axis)){
            img <- .mirror.image(img, image.mirror.axis)
        }
        if (!is.null(image.rotate.degree)){
            img <- .rotate.image(img, image.rotate.degree)
        }
    }


    annotation_custom(grob = grid::rasterGrob(img),
                      xmin = 1,
                      ymin = 1,
                      xmax = dim(img)[2],
                      ymax = dim(img)[1]
    )
}

# reference SpatialExperiment
.rotate.image <- function(img.raster, degrees){
    stopifnot(length(degrees) == 1, is.numeric(degrees), degrees%%90 == 0)
    s <- sign(degrees)
    f <- ifelse(s == 1, function(x) t(apply(x, 2, rev)), function(x) apply(x, 1, rev))
    n <- abs(degrees/90)
    for (i in seq_len(n)) {
        x <- f(img.raster)
    }
    as.raster(x)
}

.mirror.image <- function(img.raster, mirror.axis){
    mirror.axis <- match.arg(mirror.axis, c('h', 'v'))
    x <- switch(mirror.axis, 
                v = apply(img.raster, 2, rev),
                h = t(apply(img.raster, 1, rev)))
    as.raster(x)
}

.check_aes_exits <- function(mapping, aesthetic){
    x <- match(aesthetic, names(mapping))
    x <- x[!is.na(x)]
    if (length(x)==0){
        return(NULL)
    }else{
        x <- names(mapping)[x]
        x <- mapping[[x]]
        return(x)
    }
}

##' @importFrom yulab.utils get_fun_from_pkg
SpatialColors <- yulab.utils::get_fun_from_pkg('Seurat', 'SpatialColors')
