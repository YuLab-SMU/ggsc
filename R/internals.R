.buildWkde <- function(w, coords, n = 400, joint = FALSE, joint.fun = prod){
   rlang::check_installed(c('ks', 'Matrix'), 'for the 2D weighted kernel density estimation.')
   if (inherits(w, 'matrix')){
     w <- Matrix::Matrix(w, sparse = TRUE)
   }
   lims <- c(range(coords[,1]), range(coords[,2]))
   h <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
   res <- CalWkdeCpp(x=as.matrix(coords), w=w, l=lims, h = h, n = n)
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
