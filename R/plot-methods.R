##' @title plot_lisa_feature
##' @rdname plot-lisa-feature
##' @param spe SpatialExperiment or SingleCellExperiment object.
##' @param lisa.res the result returned by \code{SVP::runLISA()}.
##' @param features selected features to be visualized, default is NULL.
##' @param assay.type the assay name where data will be used from
##' (e.g., 'data', 'counts'), default is \code{'logcounts'}.
##' @param geom the function of geometric layer, default is \code{geom_bgpoint}, 
##' other option is \code{sc_geom_point}.
##' @param pointsize numeric the size of point, default is \code{2}.
##' @param hlpointsize numeric the size of point which contains corresbonding
##' spatially variable gene(i.e., SVG), default is \code{1.8}.
##' @param clustertype cell type which is from the result of \code{lisa.res},
##' default is \code{'High'}.
##' @param hlcolor the color of circular line which enfolds the point 
##' that contains SVG, default is \code{'black'}.
##' @param gap_line_width numeric the line width of gap between the background and
##' top point point layer, default is \code{.1}.
##' @param bg_line_width numeric the line width of background point layer,
##' default is \code{0.3}.
##' @param facet_name the name of facet used in \code{facet_wrap()},
##' default is \code{NULL}.
##' @param reduction reduction method, default is \code{NULL} and will 
##' use the default setting store in the object
##' @param image.plot logical whether display the image of spatial experiment, default
##' is FALSE.
##' @param label_wrap_width numeric maximum number of characters before wrapping the strip.
##' default is \code{30}.
##' @param ncol numeric the column number, default is 6.
##' @param ... additional parameters pass to \code{scattermore::geom_scattermore()}
##' \itemize{
##'     \item \code{bg_colour} the colour of background point, default is \code{NA}.
##'      this character also can be set in \code{mappint}.
##'     \item \code{alpha} the transparency of colour, default is 1.
##'  }
##' @return ggplot object
##' @importFrom ggplot2 theme element_rect label_wrap_gen
##' @importFrom stats as.formula
##' @export
##' @examples
##' \dontrun{
##' library(ggplot2)
##' library(SingleCellExperiment) |> suppressPackageStartupMessages()
##' library(SpatialExperiment) |> suppressPackageStartupMessages()
##' library(STexampleData)
##' # create ExperimentHub instance
##' eh <- ExperimentHub()
##' # query STexampleData datasets
##' myfiles <- query(eh, "STexampleData")
##' ah_id <- myfiles$ah_id[myfiles$title == 'Visium_humanDLPFC']
##' spe <- myfiles[[ah_id]]
##' spe <- spe[, colData(spe)$in_tissue == 1]
##' spe <-scater::logNormCounts(spe)
##' genes <- c('MOBP', 'PCP4', 'SNAP25', 'HBB', 'IGKC', 'NPY')
##' target.features <- rownames(spe)[match(genes, rowData(spe)$gene_name)]
##' library(SVP)
##' lisa.res1 <- runLISA(spe,
##'                      assay.type='logcounts',
##'                      features=target.features[seq(2)],
##'                      weight.method='knn', 
##'                      k=50)
##' plot_lisa_feature(spe, lisa.res=lisa.res1, features=target.features[seq(2)],
##'                   pointsize=2, hlpointsize=2, gap_line_width=.1)
##' }
plot_lisa_feature <- function(spe,
                         lisa.res,
                         features = NULL,
                         assay.type = 'logcounts',
                         geom = geom_bgpoint,
                         pointsize = 2,
                         hlpointsize = 1.8,
                         clustertype = 'High',
                         hlcolor = c('black'),
                         gap_line_width = .1,
                         bg_line_width = .3,
                         facet_name = NULL,
                         reduction = NULL,
                         image.plot = FALSE,
                         label_wrap_width = 30,
                         ncol = 6,
                         ...
                         ){
    if (missing(lisa.res) || is.null(lisa.res)){
        if (is.null(features)){
            cli::cli_abort("The {.var features} should not be `NULL`, when {.var lisa.res} is missing or NULL.")
        }else{
            features.nm <- features
        }
    }else if(inherits(lisa.res, 'SimpleList') || inherits(lisa.res, "list")){
        if (!is.null(features)){
            lisa.res <- lisa.res[features]
        }
        names(lisa.res) <- gsub("_", " ", names(lisa.res))
        features.nm <- names(lisa.res)
        lisa.res <- lisa.res |>
                    lapply(function(x)x|>tibble::rownames_to_column(var='.BarcodeID')) |>
                    dplyr::bind_rows(.id='features') |>
                    dplyr::mutate(features = factor(.data$features, levels=features.nm)) 
    }
    rownames(spe) <- gsub("_", " ", rownames(spe))
    if (is.null(reduction)){
        cnm <- SpatialExperiment::spatialCoordsNames(spe)
        p <- sc_spatial(spe,
                        features.nm,
                        mapping = aes(x=!!rlang::sym(cnm[1]), y = !!rlang::sym(cnm[2])),
                        pointsize = pointsize,
                        slot = assay.type,
                        gap_colour = NA,
                        image.plot = image.plot,
                        geom = geom,
                        ...
        )
    }else{
        p <- sc_feature(
               spe,
               features = features.nm,
               reduction = reduction,
               geom = geom,
               pointsize = pointsize,
               slot = assay.type,
               ...
        )
    }
    if (missing(lisa.res) || is.null(lisa.res)){
        message("The lisa result is not provided")
        return(p)
    }

    if (inherits(p, 'patchwork')){
        `%add+%` <- `&`
    }else{
        `%add+%` <- `+`
    }
    p1 <- p %add+% sc_geom_annot(
                    data = lisa.res,
                    mapping = aes(bg_colour = !!rlang::sym("cluster.test"), subset = !!rlang::sym("cluster.test") %in% clustertype),
                    pointsize = hlpointsize,
                    gap_line_width = gap_line_width,
                    bg_line_width = bg_line_width
          ) %add+%
          scale_bg_colour_manual(
            values = hlcolor,
            guide = guide_legend(
              theme = theme(
                legend.title = element_text(size = 8),
                legend.text = element_text(size = 6),
                legend.key.width  = grid::unit(.3, "cm"),
                legend.key.height  = grid::unit(.3, "cm")
              ),
              order = 1
            )
          ) %add+%
          guides(colour = guide_colorbar(
            theme = theme(
              legend.title = element_text(size=8),
              legend.text = element_text(size=6),
              legend.key.width  = grid::unit(.4, "cm"),
              legend.key.height = grid::unit(1.5, "cm")
              )
            )
          ) %add+%
          theme(
            strip.background.x=element_rect(color="white")
          )
    spe$sample_id |> unique() |> length() -> len
    if (len > 1 && !is.null(facet_name)){
        if (length(facet_name) > 1){
           tmpf <- paste0(facet_name, collapse="~") |> as.formula()
        }else{
           tmpf <- as.formula("~sample_id")
        }
    }else{
        tmpf <- as.formula("~features")
    }
    p1 <- p1 %add+%
              facet_wrap(tmpf, labeller = label_wrap_gen(label_wrap_width), ncol=ncol) %add+%
              theme(strip.background.x=element_rect(color="white"))
    return(p1)
}
    	
