#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//////////////////////////////////////// Two-stage Hybrid Bootstrap ////////////////////////////////////////

// [[Rcpp::export]]
arma::vec fc_estbth_BASErhohat(int Me, int Mpi, arma::vec tauseq, int N, int T, arma::vec y0, arma::mat Z, arma::mat Xbreve, arma::mat W1, arma::mat W2, arma::mat W3, arma::mat secA1, arma::mat secA2, arma::mat Bhatfo_upd, arma::vec phihatf, arma::vec thetahat, arma::vec muhat_T, arma::vec Vhat, double rho_ini) {
  // fr_estbth_BASErhohat in DSPM_BTS_fr.R, estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  Environment OPT = Environment::namespace_env("QDSPM");
  Function f = OPT["fr_estbth_BASErhohat"];
  NumericVector estB_mid = f(Named("Me")=Me, Named("Mpi")=Mpi, Named("tauseq")=tauseq, Named("N")=N, Named("T")=T, Named("y0")=y0, Named("Z")=Z, Named("Xbreve")=Xbreve, Named("W1")=W1, Named("W2")=W2, Named("W3")=W3, Named("secA1")=secA1, Named("secA2")=secA2, Named("Bhatfo_upd")=Bhatfo_upd, Named("phihatf")=phihatf, Named("thetahat")=thetahat, Named("muhat_T")=muhat_T, Named("Vhat")=Vhat, Named("rho_ini")=rho_ini);
  arma::vec estB(estB_mid.begin(), estB_mid.size(), false);
  return estB;
}

// [[Rcpp::export]]
arma::mat fc_estBnum_BASErhohat(int Bnum, int dimALL, int Me, int Mpi, arma::vec tauseq, int N, int T, arma::vec y0, arma::mat Z, arma::mat Xbreve, arma::mat W1, arma::mat W2, arma::mat W3, arma::mat secA1, arma::mat secA2, arma::mat Bhatfo_upd, arma::vec phihatf, arma::vec thetahat, arma::vec muhat_T, arma::vec Vhat, double rho_ini) {
  // estimates in Bnum bootstrap, the two-stage hybrid bootstrapping procedure
  arma::mat estBnum(dimALL, Bnum); estBnum.fill(0.0);
  for(int b = 0; b < Bnum; b++) {
    estBnum.col(b) = fc_estbth_BASErhohat(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini);
  }
  return estBnum;
}

// [[Rcpp::export]]
arma::vec fc_estbth_BASErhocheck(int Me, int Mpi, arma::vec tauseq, int N, int T, arma::vec y0, arma::mat Z, arma::mat Xbreve, arma::mat W1, arma::mat W2, arma::mat W3, arma::mat secA1, arma::mat secA2, arma::mat Bcheck_upd, arma::vec phicheckf, arma::vec thetacheck, arma::vec mucheck_T, arma::vec Vcheck, double rho_ini) {
  // fr_estbth_BASErhocheck in DSPM_BTS_fr.R, estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  Environment OPT = Environment::namespace_env("QDSPM");
  Function f = OPT["fr_estbth_BASErhocheck"];
  NumericVector estB_mid = f(Named("Me")=Me, Named("Mpi")=Mpi, Named("tauseq")=tauseq, Named("N")=N, Named("T")=T, Named("y0")=y0, Named("Z")=Z, Named("Xbreve")=Xbreve, Named("W1")=W1, Named("W2")=W2, Named("W3")=W3, Named("secA1")=secA1, Named("secA2")=secA2, Named("Bcheck_upd")=Bcheck_upd, Named("phicheckf")=phicheckf, Named("thetacheck")=thetacheck, Named("mucheck_T")=mucheck_T, Named("Vcheck")=Vcheck, Named("rho_ini")=rho_ini);
  arma::vec estB(estB_mid.begin(), estB_mid.size(), false);
  return estB;
}

// [[Rcpp::export]]
arma::mat fc_estBnum_BASErhocheck(int Bnum, int dimALL, int Me, int Mpi, arma::vec tauseq, int N, int T, arma::vec y0, arma::mat Z, arma::mat Xbreve, arma::mat W1, arma::mat W2, arma::mat W3, arma::mat secA1, arma::mat secA2, arma::mat Bcheck_upd, arma::vec phicheckf, arma::vec thetacheck, arma::vec mucheck_T, arma::vec Vcheck, double rho_ini) {
  // estimates in Bnum bootstrap, the two-stage hybrid bootstrapping procedure
  arma::mat estBnum(dimALL, Bnum); estBnum.fill(0.0);
  for(int b = 0; b < Bnum; b++) {
    estBnum.col(b) = fc_estbth_BASErhocheck(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini);
  }
  return estBnum;
}

