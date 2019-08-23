#include <RcppArmadillo.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::depends(RcppArmadillo)]]

typedef std::vector<double> stdvec;

// [[Rcpp::export]]
List sgccak_cpp(List A, 
  NumericMatrix C, 
  NumericVector c1,
  String scheme, 
  LogicalVector scale, 
  double long tol,
  String init, 
  LogicalVector bias, 
  LogicalVector verbose) {
  
  // Initialize variables
  int J = A.size();
  NumericVector pjs;
  for(int i = 0; i < J; ++i) {
    NumericMatrix x = A[i];
    pjs.push_back(x.ncol());
  }
  
  // initialize loadings
  R_xlen_t n = J;
  List a(n);
  if (init == "svd") {
  for(int i = 0; i < J; ++i){
    const arma::mat& X = A[i];
    arma::mat U, V;
    arma::vec S;
    arma::svd(U, S, V, X, "standard");
    arma::vec LOADINGS = V.col(1);
    a[i] = arma::conv_to<stdvec>::from(LOADINGS);
  }
  } else if (init == "random") {
    // a <- lapply(pjs, rnorm)
  for(int i = 0; i < J; ++i){
    a[i] = rnorm(pjs[i]);
  }
  } else {
    stop("init should be either random or svd.");
  }
  
  return List::create(Named("J", J), 
    Named("pjs", pjs),
    Named("a", a));
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
A <- list(mrna = matrix(rnorm(50*100), nrow = 50, ncol = 100),
  mirna = matrix(rnorm(50*100), nrow = 50, ncol = 100))
C <- matrix(c(0, 1, 1, 0), nrow = 2, ncol = 2)
c1 <- c(0.5, 0.5)
scheme = "centroid"
scale = TRUE
init = "random"
tol = 0.0001
bias = TRUE
verbose = TRUE
result <- sgccak_cpp(A, C, c1, scheme, scale, tol, init, bias, verbose)
result



*/
