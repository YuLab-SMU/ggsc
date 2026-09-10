## Convert a 'DataFrame' or 'matrix' into a 'data.frame' without touching
## the column names (i.e. what `check.names = FALSE` used to do).
##
## `as.data.frame()` no longer accepts `check.names` for 'DataFrame' objects
## (recent S4Vectors uses `make.names` instead and raises an error for
## `check.names`), so the original column names are restored afterwards,
## which works with all S4Vectors versions.
.as_data_frame <- function(x){
    nms <- colnames(x)
    df <- as.data.frame(x)
    if (!is.null(nms) && length(nms) == ncol(df)){
        colnames(df) <- nms
    }
    df
}

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

    xx <- colData(object) |> .as_data_frame() |> suppressWarnings()
    reduced.dat <- NULL
    if (!is.null(dims)){
        if (length(reducedDimNames(object)) == 0){
            cli::cli_abort(c("The {.cls {class(object)}} didn't contain the results of reduction."))
        }
        if (is.null(reduction)){
            reduction <- 1
        }
        reduced.dat <- reducedDims(object)[[reduction]][,dims] |>
            .as_data_frame()
        xx <- cbind(reduced.dat, xx)
    }

    if (!is.null(features)){
        if (slot == 'data'){
            if ('logcounts' %in% assayNames(object)){
                slot <- 'logcounts'
            }else{
                slot <- 1
            }
        }

        tmp <- .FetchDataFromSCE(object, features, assay.type=slot)

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
                 .as_data_frame()
        }
        xx <- cbind(xx[,!colnames(xx) %in% colnames(tmp),drop=FALSE], tmp)
    }
    xx <- cbind(data.frame(.BarcodeID=rownames(xx)), xx)
    return(xx)
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
        if (is.numeric(features)){
            features <- features[features <= nrow(object)]
            features <- rownames(object)[features]
        }
        if ('layer' %in% names(formals(utils::getS3method('FetchData', 'Seurat', envir = asNamespace('SeuratObject'))))) {
            tmp <- SeuratObject::FetchData(object, vars = features, cells = cells, layer = slot)
        } else {
            tmp <- SeuratObject::FetchData(object, vars = features, cells = cells, slot = slot)
        }
        xx <- xx[, !colnames(xx) %in% colnames(tmp),drop=FALSE]
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
    xx <- cbind(data.frame(.BarcodeID=rownames(xx)), xx)
    return(xx)
}


.FetchDataFromSCE <- function(x, features, assay.type = 1){
    if (is.numeric(features)){
        features <- features[features <= nrow(x)]
        features <- rownames(x)[features]
    }
    f1 <- intersect(features, rownames(x))
    y <- NULL
    if (length(f1) > 0){
        y <- assay(x, assay.type)[f1, ,drop=FALSE]
    }

    features <- setdiff(features, rownames(x))
    if (length(features)==0){
        if (is.null(y)){
        cli::cli_abort("The {.var features} is/are not in the {.cls class(x)}.")
        }else{
            return(y)
        }
    }

    meta.data <- colData(x) |>
                 .as_data_frame() |>
                 suppressWarnings()

    nm2 <- lapply(seq(ncol(meta.data)), function(x)is.numeric(meta.data[,x])) |>
        unlist()  
    nm2 <- colnames(meta.data)[nm2]
    
    f2 <- intersect(features, nm2)
    if (length(f2) > 0){
        y2 <- meta.data[,f2,drop=FALSE] |> t()
        y <- rbind(y, y2)
    }

    features <- setdiff(features, nm2)

    if (length(features) == 0){
        return(y)
    }

    rds <- reducedDims(x)
    rdsnm <- lapply(rds, colnames) 
    f3 <- intersect(features, unlist(rdsnm))
    
    if (length(f3) > 0){
        ind1 <- lapply(rdsnm, function(x) f3 %in% x)
        ind2 <- lapply(ind1, any) |> unlist()
        y3 <- lapply(rds[ind2], function(x)x[,colnames(x) %in% f3,drop=FALSE])
        y3 <- do.call('cbind', y3) |> t()
        y <- rbind(y, y3)
    }
    return(y)
}
