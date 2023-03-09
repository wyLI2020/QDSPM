#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fc_Rrho(arma::mat Lambda1, arma::vec lambda1, arma::mat aTa, arma::vec kappa) {
  // R_{NT}(\rho) = (\Lambda_{NT} \kappa - \lambda_{NT})^{\prime} a^{\prime}a (\Lambda_{NT} \kappa - \lambda_{NT})
  arma::mat R1 = (Lambda1 * kappa - lambda1).t() * aTa * (Lambda1 * kappa - lambda1);
  return R1;
}

