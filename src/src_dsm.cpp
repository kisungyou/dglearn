/* ROUTINES FOR DSM-RELATED FUNCTIONS
 * 
 * (01) src_dsm_dist2diff  : scaled distance between two diffusion operators
 * (02) src_dsm_spheremean : compute a frechet mean on the sphere 
 * 
 * 
 * ---------------------------------------------------------------------------
 * special routines for sphere-valued row vectors.
 * 
 * sphere_proj : projection onto the tangent space. 
 * sphere_dist : geodesic distance
 * sphere_exp  : exponential map
 * sphere_log  : logarithmic map
 * 
 */


#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace std;


// =============================================================================
arma::rowvec sphere_proj(arma::rowvec x, arma::rowvec u){
  return(u-x*(arma::dot(x,u)));
}
double sphere_dist(arma::rowvec x, arma::rowvec y){
  arma::rowvec vecxy = x-y;
  double dotxy = arma::dot(x, y);
  
  if (arma::norm(vecxy, 2) < arma::datum::eps){
    return(0.0);
  } else if (std::sqrt(dotxy*dotxy) >= (1.0-arma::datum::eps)){
    return(arma::datum::pi);
  } else {
    return(std::acos(arma::dot(x, y)));  
  } 
}
arma::rowvec sphere_exp(arma::rowvec x, arma::rowvec d, double t){
  double nrm_td = arma::norm(t*d, 2); // theta
  arma::rowvec out;
  if (nrm_td < 1e-15){ // very close
    out = x;
  } else {
    out = std::cos(nrm_td)*x + ((std::sin(nrm_td))/nrm_td)*t*d;
    out /= arma::norm(out, 2);
  }
  return(out);
}
arma::rowvec sphere_log(arma::rowvec x, arma::rowvec y){
  arma::rowvec v = sphere_proj(x,y-x);
  double di = sphere_dist(x,y);
  if (di > 1e-6){
    double nv = arma::norm(v, 2);
    v = v*(di/nv);
  }
  return(v);
}
arma::rowvec sphere_initialize(arma::mat X, arma::vec weight){
  int N = X.n_rows;
  int P = X.n_cols;
  
  arma::rowvec mvec(P,fill::zeros);
  for (int n=0; n<N; n++){
    mvec += X.row(n)*weight(n);
  }
  arma::rowvec output = mvec/arma::norm(mvec, 2);
  return(output);
}


// =============================================================================
// (01) src_dsm_dist2diff 
// [[Rcpp::export]]
double src_dsm_dist2diff(arma::mat &Px, arma::mat &Py){
  int m = Px.n_rows;
  double mm = static_cast<double>(m);
  double pi = arma::datum::pi;
  double epsthr = 100*arma::datum::eps;
  
  double tmpval = 0.0;
  double tmpsum = 0.0;
  arma::rowvec sqrtPx(m,fill::zeros);
  arma::rowvec sqrtPy(m,fill::zeros);
  
  arma::rowvec normPx(m,fill::zeros);
  arma::rowvec normPy(m,fill::zeros);
  
  for (int i=0; i<m; i++){
    // normalize the rows in case
    normPx = Px.row(i)/arma::accu(Px.row(i));
    normPy = Py.row(i)/arma::accu(Py.row(i));
    
    // take square roots
    sqrtPx = arma::sqrt(normPx);
    sqrtPy = arma::sqrt(normPy);
    
    // contribute to the summation
    if (arma::norm(sqrtPx-sqrtPy,2) > epsthr){
      tmpval = std::acos(arma::dot(sqrtPx, sqrtPy));
      tmpsum += tmpval*tmpval; 
    }
  }
  
  double output = std::sqrt(4.0*tmpsum/(mm*pi*pi));
  return(output);
}




// =============================================================================
// (02) src_dsm_spheremean
// [[Rcpp::export]]
arma::rowvec src_dsm_spheremean(arma::mat &data, arma::vec &weight, int maxiter, double abstol){
  // PREPARE
  int N = data.n_rows;
  int P = data.n_cols;
  
  arma::rowvec Stmp(P,fill::zeros);
  arma::rowvec Sold = sphere_initialize(data, weight);
  arma::rowvec Snew(P,fill::zeros);
  
  // MAIN
  double Sinc = 0.0;
  for (int it=0; it<maxiter; it++){
    // MAIN-1. compute the gradient
    Stmp.fill(0.0);
    for (int n=0; n<N; n++){
      Stmp += 2.0*weight(n)*sphere_log(Sold, data.row(n));
    }
    
    // MAIN-2. compute an update
    Snew = sphere_exp(Sold, Stmp, 1.0);
    
    // MAIN-3. decision rule
    Sinc = sphere_dist(Sold, Snew);
    
    // MAIN-4. what should we do?
    Sold = Snew;
    if (Sinc < abstol){
      break;
    }
  }
  
  // RETURN
  return(Sold);
}