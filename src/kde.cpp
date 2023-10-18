#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

Rcpp::NumericVector Quantile(Rcpp::NumericVector x, Rcpp::NumericVector probs) {
  const size_t n=x.size(), np=probs.size();
  if (n==0) return x;
  if (np==0) return probs;
  Rcpp::NumericVector index = (n-1.)*probs, y=x.sort(), x_hi(np), qs(np);
  Rcpp::NumericVector lo = Rcpp::floor(index), hi = Rcpp::ceiling(index);

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

arma::vec extractDensity(arma::mat x, arma::vec gx, arma::vec gy, arma::mat z){
  arma::uvec newx = findIntervalCpp(x.col(0), gx);
  arma::uvec newy = findIntervalCpp(x.col(1), gy);
  arma::mat sz = z.submat(newx, newy);
  arma::vec res = sz.diag();
  return (res);
}

double BandwidthNrdCpp(arma::vec x){
    NumericVector p = {0.25, 0.75};
    NumericVector r = Quantile(as<NumericVector>(wrap(x)), p);
    double h = (r[1] - r[0])/1.349;
    double w = pow(x.n_elem, -0.2);
    double s = sqrt(var(x));
    double v = 4 * 1.06 * std::min(s, h) * w;
    return (v);
}

arma::vec Kde2dWeightedCpp(arma::mat x,
                   arma::rowvec w,
                   arma::vec gx,
                   arma::vec gy,
                   arma::vec h
                   ){

    int nx = x.n_rows;
    int n = gx.n_elem;

    w = w/sum(w) * w.n_elem;

    NumericMatrix ax = outer(as<NumericVector>(wrap(gx)), as<NumericVector>(wrap(x.col(0))), std::minus<double>());
    NumericMatrix ay = outer(as<NumericVector>(wrap(gy)), as<NumericVector>(wrap(x.col(1))), std::minus<double>());

    ax = ax / h[0];
    ay = ay / h[1];

    NumericVector v = rep_each(as<NumericVector>(wrap(w)), n);
    NumericVector dax = Rcpp::dnorm(as<NumericVector>(ax));

    NumericVector day = Rcpp::dnorm(as<NumericVector>(ay));
    day.attr("dim") = Dimension(n, nx);
    NumericMatrix daym = as<NumericMatrix>(day);
    daym = transpose(daym);


    NumericVector vx = dax * v;
    vx.attr("dim") = Dimension(n, nx);
    NumericMatrix u = as<NumericMatrix>(vx);

    arma::mat z = (as<arma::mat>(u) * as<arma::mat>(daym))/(sum(w) * h[0] * h[1]);
    
    arma::vec res = extractDensity(x, gx, gy, z);
    return (res);
}

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

  for (uword i = 0; i < w.n_rows; i++){
    result.col(i) = Kde2dWeightedCpp(x, wv.row(i), gx, gy, H);
  }

  return (result);
}
