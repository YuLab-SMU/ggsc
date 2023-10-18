.buildWkde <- function(w, coords, n = 400){
   rlang::check_installed(c('ks', 'Matrix'), 'for the 2D weighted kernel density estimation.')
   if (inherits(w, 'matrix')){
       w <- Matrix::Matrix(w, sparse=TRUE)
   }
   lims <- c(range(coords[,1]), range(coords[,2]))
   h <- c(ks::hpi(coords[,1]), ks::hpi(coords[,2]))
   res <- CalWkdeCpp(x=as.matrix(coords), w=w, l=lims, h = h, n = n)
   colnames(res) <- rownames(w)
   rownames(res) <- colnames(w)
   return(res)
}
