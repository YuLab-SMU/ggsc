.buildWkde <- function(w, coords, n = 400, joint = FALSE, joint.fun = prod){
   rlang::check_installed(c('ks', 'Matrix'), 'for the 2D weighted kernel density estimation.')
   if (inherits(w, 'matrix')){
     w <- Matrix::Matrix(w, sparse = TRUE)
   }
   lims <- c(range(coords[,1]), range(coords[,2]))
   h <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
   res <- CalWkdeCpp(x=as.matrix(coords[,seq(2),drop=FALSE]), w=w, l=lims, h = h, n = n)
   colnames(res) <- rownames(w)
   rownames(res) <- colnames(w)
   if (joint && !is.null(joint.fun)){
     oldcnm <- colnames(res)
     clnm <- paste(colnames(res), collapse="+")
     joint.res <- apply(res, 1, joint.fun) 
     res <- cbind(res, joint.res)
     colnames(res) <- c(oldcnm, clnm)
   }
   return(res)
}

.split.by.feature <- function(p, ncol, joint = FALSE){
   rlang::check_installed('aplot', 'for split ggplot object by features.')
   p <- p$data |> dplyr::group_split(.data$features) |>
           lapply(function(i){
              p$data <- i
              p <- .add_class(p, "ggsc")
              return(p)
            })

   indx <- length(p)
   
   if (joint){
      p[[indx]] <- p[[indx]] + ggplot2::labs(colour = 'joint_density')
   }
   
   p <- aplot::plot_list(gglist = p, ncol = ncol)
   return(p)
}

#' @importFrom ggfun get_aes_var
.cal_pie_radius <- function(data, mapping){
   x <- ggfun::get_aes_var(mapping, 'x') 
   y <- ggfun::get_aes_var(mapping, 'y')
   r = (max(data[[x]], na.rm=TRUE) - min(data[[x]], na.rm=TRUE)) * (max(data[[y]], na.rm=TRUE) - min(data[[y]], na.rm=TRUE))
   r = sqrt(r / nrow(data) / pi) * .85
   return(r)
}

.cal_ratio <- function(data, mapping){
   x <- ggfun::get_aes_var(mapping, 'x')
   y <- ggfun::get_aes_var(mapping, 'y')
   1*max(data[[x]], na.rm=TRUE)/max(data[[y]], na.rm=TRUE)
}

.set_default_cols <- function(n){
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
              "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
              "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
    grDevices::colorRampPalette(col2)(n)
}

# This was refering to the stat_ellipse of ggplot2
.calculate_ellipse <- function(data, vars, type = 't', level= NULL, segments=50){
  dfn <- 2
  dfd <- nrow(data) - 1
  if (is.null(level)){
     level <- .9
  }

  if (!type %in% c("t", "norm", "euclid")) {
    cli::cli_inform("Unrecognized ellipse type")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else if (dfd < 3) {
    cli::cli_inform("Too few points to calculate an ellipse")
    ellipse <- matrix(NA_real_, ncol = 2)
  } else {
    if (type == "t") {
      v <- MASS::cov.trob(data[,vars])
    } else if (type == "norm") {
      v <- stats::cov.wt(data[,vars])
    } else if (type == "euclid") {
      v <- stats::cov.wt(data[,vars])
      v$cov <- diag(rep(min(diag(v$cov)), 2))
    }
    shape <- v$cov
    center <- v$center
    chol_decomp <- chol(shape)
    if (type == "euclid") {
      radius <- level/max(chol_decomp)
    } else {
      radius <- sqrt(dfn * stats::qf(level, dfn, dfd))
    }
    angles <- (0:segments) * 2 * pi/segments
    unit.circle <- cbind(cos(angles), sin(angles))
    ellipse <- t(center + radius * t(unit.circle %*% chol_decomp))
  }

  colnames(ellipse) <- vars
  res <- stats::cov.wt(ellipse)
  res <- res$center |> as.matrix() |> t() |> data.frame()
  return(res)
}


#' @importFrom ggplot2 ggplot_build
.check_colour <- function(x, y){
  gb <- ggplot_build(x)
  flag1 <- gb$plot$scales$has_scale('colour')
  if (flag1){
     flag1 <- inherits(gb$plot$scales$get_scales("colour"), "ScaleContinuous")
  }
  flag2 <- any(c("color", "colour") %in% names(y$mapping)) || any(c("color", "colour") %in% names(y))
  flag1 && !flag2
}


