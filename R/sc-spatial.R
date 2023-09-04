##' @title sc_spatial
##' @rdname sc-spatial-methods
##' @param object Seurat object
##' @param features selected features to be visualized
##' @param sample.id the index name of sample id, which only work with SingleCellExperiment or SpatialExperiment.
##' @param image.id the index name of image id, which only work with SingleCellExperiment or SpatialExperiment.
##' @param slot if plotting a feature, which data will be used (e.g., 'data', 'counts'), the assay name if object
##' is SingleCellExperiment or SpatialExperiment.
##' @param image.plot whether to display the issue image as background.
##' @param image.rotate.degree integer the degree to ratate image, default is NULL.
##' @param image.mirror.axis character the direction to mirror the image, default is 'h'.
##' @param remove.point whether to remove the spot points, it is nice if your just view the issue 
##' image, default is FALSE.
##' @param mapping aesthetic mapping, default is NULL.
##' @param ncol integer number of facet columns if 'length(features) > 1', default is 6.
##' @param ... additional parameters.
##' @return ggplot object
##' @importFrom grid rasterGrob unit
##' @importFrom ggplot2 facet_grid annotation_custom
##' @importFrom ggplot2 geom_point xlab ylab geom_blank coord_fixed 
##' @importFrom ggplot2 scale_color_gradientn
##' @importFrom Seurat DefaultAssay
##' @export
setGeneric('sc_spatial', function(object, features = NULL, 
                                  sample.id = NULL, image.id = NULL, 
                                  slot = "data", image.plot = TRUE, 
                                  image.rotate.degree = NULL,
                                  image.mirror.axis = NULL,
                                  remove.point = FALSE,
                                  mapping = NULL,
                                  ncol = 6,
                                  ...) 
           standardGeneric('sc_spatial')
)

##' @importFrom grDevices as.raster
##' @rdname sc-spatial-methods
##' @aliases sc_spatial,Seurat
##' @exportMethod sc_spatial
setMethod("sc_spatial", 'Seurat', function(object, features = NULL, slot = "data", image.plot = TRUE, 
                                           image.rotate.degree = NULL, image.mirror.axis = 'h', 
                                           remove.point = FALSE, mapping = NULL, ncol = 6, ...) {
    images <- SeuratObject::Images(object = object, 
                    assay = Seurat::DefaultAssay(object = object)
                )
    img <- object@images[[images]]@image |> as.raster()
    
    coord <- SeuratObject::GetTissueCoordinates(object = object[[images]])
    
    d <- get_dim_data(object = object, features = features, dims = NULL)

    d <- cbind(coord, d)

    default_mapping <- aes_string(x = colnames(coord)[2], y = colnames(coord)[1])

    if (!is.null(features)){
        d <- tidyr::pivot_longer(d, seq(ncol(d) - length(features) + 1, ncol(d)), names_to = 'features')
        default_mapping <- modifyList(default_mapping, aes_string(color = "value"))
    }

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    p <- ggplot(d, mapping)

    if (image.plot){
        img.annot <- .build_img_annot_layer(img, image.rotate.degree, image.mirror.axis)
        p <- p + img.annot
    }

    if (!remove.point || !(is.null(features) || missing(features))){
        p <- p + geom_point()
    }else{
        p <- p + geom_blank()
    }

    p <- p +
         .feature_setting(features, ncol) +
         ylab(NULL) +
         xlab(NULL) +
         coord_fixed() +
         scale_color_gradientn(colours = SpatialColors(n=100)) +
         theme_bw2() 

    return(p)    


})

#' @importFrom SingleCellExperiment int_metadata
#' @rdname sc-spatial-methods
#' @aliases sc_spatial,SingleCellExperiment
#' @exportMethod sc_spatial
setMethod('sc_spatial', 'SingleCellExperiment', function(object, 
                                                         features = NULL, 
                                                         sample.id = NULL, 
                                                         image.id = NULL, 
                                                         slot = 1,
                                                         image.plot = TRUE,
                                                         image.rotate.degree = NULL,
                                                         image.mirror.axis = 'h',
                                                         remove.point = FALSE,
                                                         mapping = NULL,
                                                         ncol = 6,
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

    features.da <- .extract_sce_data(object, features = features, dims = NULL, cells = NULL, slot = slot)

    d <- merge(coords.da, features.da, by = 0)
    rownames(d) <- d$Row.names
    d$Row.names <- NULL
    
    default_mapping <- aes_string(x = colnames(coords.da)[2], y = colnames(coords.da)[1])
    if (!is.null(features)){
        d <- tidyr::pivot_longer(d, seq(ncol(d) - length(features) + 1, ncol(d)), names_to = 'features')
        default_mapping <- modifyList(default_mapping, aes_string(color = "value"))
    }

    if (!is.null(mapping)){
        mapping <- modifyList(default_mapping, mapping)
    }else{
        mapping <- default_mapping
    }

    p <- ggplot(d, mapping) 

    if (image.plot){
        img.annot <- .build_img_annot_layer(img.da, image.rotate.degree, image.mirror.axis)
        p <- p + img.annot
    }

    if (!remove.point || !(is.null(features) || missing(features))){
        p <- p + geom_point() 
    }else{
        p <- p + geom_blank()
    }
    
    p <- p +
         .feature_setting(features, ncol) +
         ylab(NULL) +
         xlab(NULL) +
         coord_fixed() +
         scale_color_gradientn(colours = SpatialColors(n=100))

    return(p)
})

#' @importFrom SingleCellExperiment int_metadata int_colData
.extract_img_data <- function(x, sample.id = NULL, image.id = NULL){
    img.da <- int_metadata(x)[['imgData']]
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
    x <- int_colData(x)
    x <- x[['spatialCoords']] * image.da$scaleFactor
    return(x)
}

.build_img_annot_layer <- function(image.da, image.rotate.degree = NULL, image.mirror.axis = NULL){
    if (inherits(image.da, 'raster')){
        img <- image.da
    }else{
        img <- image.da[['data']][[1]] |> as.raster()
    }

    if (!is.null(image.rotate.degree)){
        img <- .rotate.image(img, image.rotate.degree)
    }
    if (!is.null(image.mirror.axis)){
        img <- .mirror.image(img, image.mirror.axis)
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
                h = apply(img.raster, 2, rev),
                v = t(apply(img.raster, 1, rev)))
    as.raster(x)
}

##' @importFrom yulab.utils get_fun_from_pkg
SpatialColors <- yulab.utils::get_fun_from_pkg('Seurat', 'SpatialColors')
