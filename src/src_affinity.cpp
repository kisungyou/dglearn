/* ROUTINES FOR AFFINITY MATRIX COMPUTATION
 * 
 * (01) src_affgaussian : affinity construction using Gaussian kernel
 * 
 */

#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;
using namespace std;


// (01) src_affgaussian ========================================================
// [[Rcpp::export]]
arma::mat src_affgaussian(arma::mat &X, bool use_auto, double bandwidth){
  // parameters
  int n = X.n_rows;
  
  // compute pairwise distances
  arma::vec upper_dist(n*(n-1)/2, fill::zeros);
  arma::mat mat_dist(n,n,fill::zeros);
  
  int counter = 0;
  for (int i=0; i<(n-1); i++){
    for (int j=(i+1); j<n; j++){
      mat_dist(i,j) = arma::norm(X.row(i)-X.row(j), 2);
      mat_dist(j,i) = mat_dist(i,j);
      
      upper_dist(counter) = mat_dist(i,j);
      counter += 1;
    }
  }
  
  // bandwidth for exp(-|x-y|^2/(2*(sig^2)))
  double par_band = 0.0;
  if (use_auto==false){
    par_band = bandwidth;
  } else {
    par_band = arma::median(upper_dist);
  }
  
  // compute
  return(arma::exp(-(arma::pow(mat_dist, 2)/(2.0*par_band*par_band))));
}