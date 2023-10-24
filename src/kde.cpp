#include <RcppArmadillo.h>
#include <RcppParallel.h>
using namespace RcppParallel;
using namespace Rcpp;
using namespace arma;
using namespace std;

arma::vec Quantile(arma::vec x, arma::vec probs) {
  const size_t n=x.n_elem, np=probs.n_elem;
  if (n==0) return x;
  if (np==0) return probs;
  arma::vec index = (n-1.0)*probs, y=sort(x), x_hi(np), qs(np);
  arma::vec lo = arma::floor(index), hi = arma::ceil(index);

  for (size_t i=0; i<np; ++i) {
    qs[i] = y[lo[i]];
    x_hi[i] = y[hi[i]];
    if ((index[i]>lo[i]) && (x_hi[i] != qs[i])) {
      double h;
      h = index[i]-lo[i];
      qs[i] = (1.-h)*qs[i] + h*x_hi[i];
    }
  }
  return qs;
}

uvec findIntervalCpp(arma::vec x, arma::vec breaks) {
  uvec out(x.size());

  vec::iterator it, pos;
  uvec::iterator out_it;

  for(it = x.begin(), out_it = out.begin(); it != x.end(); 
      ++it, ++out_it) {
    pos = std::upper_bound(breaks.begin(), breaks.end(), *it);
    *out_it = std::distance(breaks.begin(), pos);
  }
  return (out - 1);
}

double BandwidthNrdCpp(arma::vec x){
    arma::vec p = {0.25, 0.75};
    arma::vec r = Quantile(x, p);
    double h = (r[1] - r[0])/1.349;
    double w = pow(x.n_elem, -0.2);
    double s = sqrt(var(x));
    double v = 4 * 1.06 * std::min(s, h) * w;
    return (v);
}

arma::vec Kde2dWeightedCpp(arma::mat x,
                   arma::rowvec w,
                   arma::mat ax,
                   arma::mat ay,
                   arma::vec h,
                   arma::uvec indx,
                   arma::uvec indy
                   ){
    int n = ax.n_rows;

    w = w/sum(w) * w.n_elem;

    ax = ax / h[0];
    ay = ay / h[1];

    arma::mat v = repelem(w, n, 1);
    arma::mat u = arma::normpdf(ax) % v;
    arma::mat day = arma::normpdf(ay) % v;
    arma::mat daym = day.t();

    arma::mat z = (u * daym)/(accu(v) * h[0] * h[1]);
    
    arma::mat sz = z.submat(indx, indy);
    arma::vec res = sz.diag();
    return (res);
}

arma::mat outergrid(arma::vec grid, arma::vec x){
    arma::mat gxm = repelem(grid, 1, x.n_elem);
    arma::mat xm = repelem(x, 1, grid.n_elem);

    arma::mat ax = gxm - xm.t();
    
    return(ax);
}


struct CalWkde : public Worker{
  const arma::mat& x;
  const arma::mat& w;
  const arma::mat& ax;
  const arma::mat& ay;
  const arma::vec& H;
  const arma::uvec& indx;
  const arma::uvec& indy;

  mat& result;

  CalWkde(const arma::mat& x, const arma::mat& w, const arma::mat& ax,
         const arma::mat ay, const arma::vec& H, const arma::uvec& indx, 
         const arma::uvec& indy, mat& result)
  : x(x), w(w), ax(ax), ay(ay), H(H), indx(indx), indy(indy), result(result) { }

  void operator()(std::size_t begin, std::size_t end){
    for (uword i = begin; i < end; i++){
        result.col(i) = Kde2dWeightedCpp(x, w.row(i), ax, ay, H, indx, indy);
    }
  }
};


//' Two-Dimensional Weighted Kernel Density Estimation And Mapping the Result To Original Dimension
//' @param x The 2-D coordinate matrix
//' @param w The weighted sparse matrix, the number columns the same than the number rows than x.
//' @param l The limits of the rectangle covered by the grid as c(xl, xu, yl, yu)
//' @param h The vector of bandwidths for x and y directions, defaults to normal reference bandwidth
//' (see bandwidth.nrd), A scalar value will be taken to apply to both directions (see ks::hpi).
//' @param adjust numeric value to adjust to bandwidth, default is 1.
//' @param n number of grid points in the two directions, default is 400.
// [[Rcpp::export]]
arma::mat CalWkdeCpp(arma::mat& x, arma::sp_mat& w, arma::vec& l, Nullable<NumericVector> h, 
        double adjust = 1.0, int n = 400) {

  arma::mat wv = conv_to<arma::mat>::from(w);

  arma::mat result(x.n_rows, w.n_rows);

  arma::vec gx = arma::linspace(l[0], l[1], n);
  arma::vec gy = arma::linspace(l[2], l[3], n);

  arma::vec H(x.n_cols);
  if (h.isNull()){
    for (uword j=0; j < x.n_cols;j ++){
       H[j] = BandwidthNrdCpp(x.col(j)) / 4 * adjust;
    }
  }else{
    H = as<arma::vec>(h);
  }

  //mapping to original coords
  arma::uvec indx = findIntervalCpp(x.col(0), gx);
  arma::uvec indy = findIntervalCpp(x.col(1), gy);
  
  arma::mat ax = outergrid(gx, x.col(0));
  arma::mat ay = outergrid(gy, x.col(1));

  uword num = wv.n_rows;
  CalWkde calWkde(x, wv, ax, ay, H, indx, indy, result);
  parallelFor(0, num, calWkde);

  return (result);
}
