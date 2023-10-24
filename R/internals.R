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

.split.by.feature <- function(p, ncol){
   p <- p$data |> dplyr::group_split(.data$features) |>
           lapply(function(i){
              p$data <-i
              return(p)
            })
   p <- aplot::plot_list(gglist = p, ncol = ncol)
   return(p)
}
