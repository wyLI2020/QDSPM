#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat fc_asmat(arma::vec vec1, int nrow, int ncol) {
  // Fill matrix with elements of vector
  arma::mat vec1_mat(nrow*ncol, 1); vec1_mat.col(0) = vec1;
  vec1_mat.reshape(nrow, ncol);
  return vec1_mat;
}

// [[Rcpp::export]]
arma::vec fc_asvec(arma::mat mat1) {
  // Matrix straighten
  int nrow = mat1.n_rows;
  int ncol = mat1.n_cols;
  mat1.reshape(nrow*ncol, 1);
  return mat1.col(0);
}

// [[Rcpp::export]]
arma::mat fc_iota(int n) {
  // \iota_{n}, n*1 matrix of ones
  arma::mat iota_n(n, 1); iota_n.fill(1.0);
  return iota_n;
}

// [[Rcpp::export]]
arma::mat fc_I(int n) {
  // I_{n}, n*n identity matrix
  arma::mat I_n(n,n); I_n.eye(n,n);
  return I_n;
}

arma::mat fc_J(int n) {
  // J_{n} = I_{n} - 1/n \iota_{n} \iota_{n}^{\prime}
  arma::mat I_n = fc_I(n);
  arma::mat iota_n = fc_iota(n);
  arma::mat J_n = I_n - 1.0/n * iota_n * iota_n.t();
  return J_n;
}

// [[Rcpp::export]]
arma::mat fc_B(double rho, int N, arma::mat W3) {
  // B(\rho) = I_{N} - \rho W_{3}
  arma::mat I_N = fc_I(N);
  arma::mat B = I_N - rho * W3;
  return B;
}

// [[Rcpp::export]]
arma::mat fc_S(double lambda, int N, arma::mat W1) {
  // S(\lambda) = I_{N} - \lambda W_{1}
  arma::mat I_N = fc_I(N);
  arma::mat S = I_N - lambda * W1;
  return S;
}

// [[Rcpp::export]]
arma::mat fc_MatrixMultip(arma::mat mat1, arma::mat mat2) {
  // matrix multiplication
  arma::mat mat1_multip_mat2 = mat1 * mat2;
  return mat1_multip_mat2;
}

// [[Rcpp::export]]
arma::mat fc_Zbar(arma::mat Z, int N, int T, int q) {
  // \bar{Z} = (\bar{z}_{1}^{\prime}, ..., \bar{z}_{N}^{\prime})^{\prime} with \bar{z}_{i} = 1/T \sum_{t=1}^{T} z_{it}
  arma::mat Zbar(N, q); Zbar.fill(0.0);
  for(int l = 0; l < q; l++) {
    arma::vec Z_l = Z.col(l);
    arma::mat Z_l_mat = fc_asmat(Z_l, N, T);
    arma::vec Zbar_l = mean(Z_l_mat, 1);
    Zbar.col(l) = Zbar_l;
  }
  return Zbar;
}

// [[Rcpp::export]]
arma::mat fc_ITWZ(int N, int q, arma::mat Z, arma::mat W) {
  // (I_{T} \otimes W) Z or (I_{T} \otimes W) Y_{-1} or (I_{T-2} \otimes W) \mathbb{Z} or (I_{T-2} \otimes W) \mathbb{Y}_{-2}
  int T = Z.n_rows / N;
  arma::mat ITWZ(N*T, q); ITWZ.fill(0.0); // (I_{T} \otimes W) Z
  for(int l = 0; l < q; l++) {
    arma::vec Z_l = Z.col(l);
    arma::mat Z_l_mat = fc_asmat(Z_l, N, T);
    arma::mat ITWZ_l_mat = W * Z_l_mat;
    arma::vec ITWZ_l = fc_asvec(ITWZ_l_mat);
    ITWZ.col(l) = ITWZ_l;
  }
  return ITWZ;
}

////////////////////////////// Estimations //////////////////////////////

// [[Rcpp::export]]
arma::vec fc_phitilde(int N, int T, arma::vec mathbbY, arma::mat mathbbZbreve, arma::mat Q) {
  // calculate an initial estimator of \phi = (\lambda, \alpha, \beta, \gamma^{\prime})^{\prime}, \widetilde{\phi} = (\breve{\mathbb{Z}}^{\prime} \mathcal{P}_{JQ,N(T-2)} \breve{\mathbb{Z}})^{-1} \breve{\mathbb{Z}}^{\prime} \mathcal{P}_{JQ,N(T-2)} \mathbb{Y}
  arma::mat J_N = fc_J(N);
  int dimQ = Q.n_cols;
  arma::mat mathcalJQ(N*(T-2),dimQ); mathcalJQ.fill(0.0); // \mathcal{J}_{N(T-2)} Q_{N(T-2)} = (I_{T-2} \otimes J_{N}) Q_{N(T-2)}
  for(int l = 0; l < dimQ; l++) {
    arma::vec Q_l = Q.col(l);
    arma::mat Q_l_mat = fc_asmat(Q_l, N, T-2);
    arma::mat mathcalJQ_l_mat = J_N * Q_l_mat;
    arma::vec mathcalJQ_l = fc_asvec(mathcalJQ_l_mat);
    mathcalJQ.col(l) = mathcalJQ_l;
  }
  arma::vec QTmathcalJmathbbY = mathcalJQ.t() * mathbbY; // Q_{N(T-2)}^{\prime} \mathcal{J}_{N(T-2)} \mathbb{Y}
  arma::mat QTmathcalJmathbbZbreve = mathcalJQ.t() * mathbbZbreve; // Q_{N(T-2)}^{\prime} \mathcal{J}_{N(T-2)} \breve{\mathbb{Z}}
  arma::mat QTmathcalJQinverse = (Q.t() * mathcalJQ).i(); // (Q_{N(T-2)}^{\prime} \mathcal{J}_{N(T-2)} Q_{N(T-2)})^{-1}
  arma::vec phitilde = (QTmathcalJmathbbZbreve.t() * QTmathcalJQinverse * QTmathcalJmathbbZbreve).i() * QTmathcalJmathbbZbreve.t() * QTmathcalJQinverse * QTmathcalJmathbbY;
  return phitilde;
}

// [[Rcpp::export]]
arma::vec fc_etaest(arma::vec phiest, arma::vec Y, arma::mat Zbreve) {
  // calculate the approximate of \eta, \widetilde{\eta} = Y - \breve{Z} \widetilde{\phi} or \widehat{\eta} = Y - \breve{Z} \widehat{\phi}_{f}
  arma::vec etaest = Y - Zbreve * phiest;
  return etaest;
}

// [[Rcpp::export]]
arma::vec fc_thetaest(arma::vec etaest, int N, int T) {
  // calculate the approximate of \theta = (\theta_{1}, ..., \theta_{N})^{\prime}, \widetilde{\theta} = (\widetilde{\theta}_{1}, ..., \widetilde{\theta}_{N})^{\prime} with \widetilde{\theta}_{i} = 1/T \sum_{t=1}^{T} \widetilde{\eta}_{it} or \widehat{\theta}
  arma::mat etaest_mat = fc_asmat(etaest, N, T);
  arma::vec thetaest = mean(etaest_mat, 1);
  return thetaest;
}

// [[Rcpp::export]]
arma::vec fc_etaestMINUSiotaTthetaest(arma::vec etaest, arma::vec thetaest, int T) {
  // calculate \widetilde{\eta} - \iota_{T} \otimes \widetilde{\theta} or \widehat{\eta} - \iota_{T} \otimes \widehat{\theta}
  arma::mat iota_T = fc_iota(T);
  arma::vec etaest_minus_iotaTthetaest = etaest - kron(iota_T, thetaest);
  return etaest_minus_iotaTthetaest;
}

// [[Rcpp::export]]
arma::vec fc_muTest(arma::vec etaest_minus_iotaTthetaest, int N, int T) {
  // calculate the approximate of (\mu_{10}, ..., \mu_{T0})^{\prime}, (\widetilde{\mu}_{1}, ..., \widetilde{\mu}_{T})^{\prime}  with \widetilde{\mu}_{t} = 1/N \sum_{i=1}^{N} (\widetilde{\eta}_{it} - \widetilde{\theta}_{i}) or (\widehat{\mu}_{1}, ..., \widehat{\mu}_{T})^{\prime}
  arma::mat etaest_minus_iotaTthetaest_mat = fc_asmat(etaest_minus_iotaTthetaest, N, T);
  arma::vec muest_T = mean(etaest_minus_iotaTthetaest_mat.t(), 1);
  return muest_T;
}

// [[Rcpp::export]]
arma::mat fc_secA1(arma::mat W3, int N) {
  // In Stage 1, A_{1} = I_{T} \otimes (W_{3}^{\prime} W_{3} - diagW3) \triangleq (I_{T} \otimes secA1) with diagW3 \triangleq diag\{w_{3,.1}^{\prime}w_{3,.1}, ..., w_{3,.N}^{\prime}w_{3,.N}\}
  arma::mat diagW3(N,N); diagW3.fill(0.0);
  for(int i = 0; i < N; i++) {
    arma::vec ithcol = W3.col(i);
    diagW3(i,i) = dot(ithcol, ithcol);
  }
  arma::mat secA1 = W3.t() * W3 - diagW3; // secA1 \triangleq (W_{3}^{\prime} W_{3} - diagW3)
  return secA1;
}

// [[Rcpp::export]]
arma::mat fc_Lambda1lambda1(arma::vec Vtilde, int N, int T, arma::mat W3, arma::mat secA1, arma::mat secA2) {
  // In Stage 1, c(\Lambda_{1}, \lambda_{1})
  arma::mat Vtilde_mat = fc_asmat(Vtilde, N, T);
  arma::mat M_Vtilde_mat = W3 * Vtilde_mat; // M \widetilde{V} = (I_{T} \otimes W_{3}) \widetilde{V}
  arma::mat A1_Vtilde_mat = secA1 * Vtilde_mat; // A_{1} \widetilde{V} = (I_{T} \otimes secA1) \widetilde{V}
  arma::mat A2_Vtilde_mat = secA2 * Vtilde_mat; // A_{1} \widetilde{V} = (I_{T} \otimes secA2) \widetilde{V}
  arma::mat A1_M_Vtilde_mat = secA1 * M_Vtilde_mat; // A_{1} M \widetilde{V}
  arma::mat A2_M_Vtilde_mat = secA2 * M_Vtilde_mat; // A_{2} M \widetilde{V}
  arma::mat A2T_M_Vtilde_mat = secA2.t() * M_Vtilde_mat; // A_{2}^{\prime} M \widetilde{V}
  arma::vec M_Vtilde = fc_asvec(M_Vtilde_mat); 
  arma::vec A1_Vtilde = fc_asvec(A1_Vtilde_mat); 
  arma::vec A2_Vtilde = fc_asvec(A2_Vtilde_mat); 
  arma::vec A1_M_Vtilde = fc_asvec(A1_M_Vtilde_mat); 
  arma::vec A2_M_Vtilde = fc_asvec(A2_M_Vtilde_mat); 
  arma::vec A2T_M_Vtilde = fc_asvec(A2T_M_Vtilde_mat); 
  arma::mat Lambda1lambda1(2,3); Lambda1lambda1.fill(0.0);
  Lambda1lambda1(0,0) = dot(Vtilde, A1_M_Vtilde) * 2.0 / (N*T*1.0);
  Lambda1lambda1(1,0) = dot(Vtilde, A2_M_Vtilde+A2T_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(0,1) = - dot(M_Vtilde, A1_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(1,1) = - dot(M_Vtilde, A2_M_Vtilde) / (N*T*1.0);
  Lambda1lambda1(0,2) = dot(Vtilde, A1_Vtilde) / (N*T*1.0);
  Lambda1lambda1(1,2) = dot(Vtilde, A2_Vtilde) / (N*T*1.0);
  return Lambda1lambda1;
}

// [[Rcpp::export]]
arma::vec fc_varepsilonest(arma::mat Best, arma::vec Vest, int N, int T) {
  // calculate the approximate of \varepsilon, \widetilde{\varepsilon} = (I_{T} \otimes B(\widehat{\rho}_{I})) \widetilde{V} or \widehat{\varepsilon} = (I_{T} \otimes B(\widehat{\rho})) \widehat{V}
  arma::mat Vest_mat = fc_asmat(Vest, N, T);
  arma::mat varepsilonest_mat = Best * Vest_mat; 
  arma::vec varepsilonest = fc_asvec(varepsilonest_mat);
  return varepsilonest;
}

// [[Rcpp::export]]
arma::vec fc_sigmavarepsilonsquareest(arma::vec varepsilonest, int N, int T) {
  // calculate the estimate of (\sigma_{\varepsilon 10}^{2}, ..., \sigma_{\varepsilon N0}^{2})^{\prime}, (\widetilde{\sigma}_{\varepsilon 1}^{2}, ..., \widetilde{\sigma}_{\varepsilon N}^{2})^{\prime} with \widetilde{\sigma}_{\varepsilon i}^{2} = 1/T \sum_{t=1}^{T} \widetilde{\varepsilon}_{it}^{2}
  arma::vec varepsilonest_square = varepsilonest % varepsilonest; // (\widetilde{\varepsilon}_{11}^{2}, ..., \widetilde{\varepsilon}_{NT}^{2})^{\prime}
  arma::mat varepsilonest_square_mat = fc_asmat(varepsilonest_square, N, T);
  arma::vec sigmavarepsilonsquareest = mean(varepsilonest_square_mat, 1); // \widetilde{\sigma}_{\varepsilon i}^{2} = 1/T \sum_{t=1}^{T} \widetilde{\varepsilon}_{it}^{2}
  return sigmavarepsilonsquareest;
}

// [[Rcpp::export]]
arma::mat fc_Sigmagest(arma::vec sigmavarepsilonsquareest, arma::mat secA1, arma::mat secA2, int N) {
  // In Stage 1, calculate the estimate of \Sigma_{g}
  arma::mat Sigmagest(2,2); Sigmagest.fill(0.0);
  double sum00 = 0.0; double sum01 = 0.0; double sum11 = 0.0; 
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      sum00 = sum00 + (secA1(i,j) + secA1(j,i)) * (secA1(i,j) + secA1(j,i)) * sigmavarepsilonsquareest(i) * sigmavarepsilonsquareest(j);
      sum01 = sum01 + (secA1(i,j) + secA1(j,i)) * (secA2(i,j) + secA2(j,i)) * sigmavarepsilonsquareest(i) * sigmavarepsilonsquareest(j);
      sum11 = sum11 + (secA2(i,j) + secA2(j,i)) * (secA2(i,j) + secA2(j,i)) * sigmavarepsilonsquareest(i) * sigmavarepsilonsquareest(j);
    }
  }
  Sigmagest(0,0) = sum00 / (2.0*N);
  Sigmagest(0,1) = sum01 / (2.0*N);
  Sigmagest(1,0) = sum01 / (2.0*N);
  Sigmagest(1,1) = sum11 / (2.0*N);
  return Sigmagest;
}

// [[Rcpp::export]]
arma::vec fc_kappacheck(arma::mat Lambda1, arma::vec lambda1) {
  // In Stage 1, calculate the linear least square estimator of \kappa = (\rho, \rho^{2})^{\prime}, \widecheck{\kappa}, using \widetilde{V}
  arma::vec kappacheck = Lambda1.i() * lambda1;
  return kappacheck;
}

// [[Rcpp::export]]
arma::vec fc_phihatf(int N, int T, int q, arma::vec Y, arma::mat Zbreve, arma::vec thetatilde, arma::vec mutilde_NT, arma::mat Bhat, arma::mat H) {
  // In Stage 1, calculate the FGSP2SLS estimator of \phi, \widehat{\phi}_{f}
  arma::mat iota_T = fc_iota(T);
  arma::vec Y_minus_iotaTthetatilde_minus_mutildeNT = Y - kron(iota_T, thetatilde) - mutilde_NT; // Y - \iota_{T} \otimes \widetilde{\theta} - \widetilde{\mu}
  arma::mat Y_minus_iotaTthetatilde_minus_mutildeNT_mat = fc_asmat(Y_minus_iotaTthetatilde_minus_mutildeNT, N, T);
  arma::mat Ystarhat_mat = Bhat * Y_minus_iotaTthetatilde_minus_mutildeNT_mat;
  arma::vec Ystarhat = fc_asvec(Ystarhat_mat); // Y_{\star}(\widehat{\rho}, \widetilde{\theta}, \widetilde{\mu}) = (I_{T} \otimes B(\widehat{\rho})) (Y - \iota_{T} \otimes \widetilde{\theta} - \widetilde{\mu})
  arma::mat Zbrevestarhat(N*T, q+3); Zbrevestarhat.fill(0.0); // \breve{Z}_{\star}(\widehat{\rho}) = (I_{T} \otimes B(\widehat{\rho})) \breve{Z}
  for(int l = 0; l < (q+3); l++) {
    arma::vec Zbreve_l = Zbreve.col(l);
    arma::mat Zbreve_l_mat = fc_asmat(Zbreve_l, N, T);
    arma::mat Zbrevestarhat_l_mat = Bhat * Zbreve_l_mat;
    arma::vec Zbrevestarhat_l = fc_asvec(Zbrevestarhat_l_mat);
    Zbrevestarhat.col(l) = Zbrevestarhat_l;
  }
  arma::vec HTYstarhat = H.t() * Ystarhat; // H^{\prime} Y_{\star}(\widehat{\rho}, \widetilde{\theta}, \widetilde{\mu})
  arma::mat HTZbrevestarhat = H.t() * Zbrevestarhat; // H^{\prime} \breve{Z}_{\star}(\widehat{\rho})
  arma::mat HTHinverse = (H.t() * H).i(); // (H^{\prime} H)^{-1}
  arma::vec phihatf = (HTZbrevestarhat.t() * HTHinverse * HTZbrevestarhat).i() * HTZbrevestarhat.t() * HTHinverse * HTYstarhat;
  return phihatf;
}  

////////////////////////////// ASD by Two-stage Hybrid Bootstrap //////////////////////////////

// [[Rcpp::export]]
arma::vec fc_varepsilonast(arma::vec east, arma::vec varepsilonhat) {
  // \varepsilon^{*} = (\varepsilon^{*}_{11}, ..., \varepsilon^{*}_{NT})^{\prime} with \varepsilon^{*}_{it} = \widehat{\varepsilon}_{it} e^{*}_{it}
  arma::vec varepsilonast = varepsilonhat % east;
  return varepsilonast;
}

// [[Rcpp::export]]
arma::mat f_Yastmat_Generation(arma::vec varepsilonast, int N, int T, int q, arma::vec y0, arma::vec phihatf, arma::vec thetahat, arma::vec muhat_T, arma::cube Zarray, arma::mat W1, arma::mat W2, arma::mat Bhat) {
  // the matrix of (y^{*}_{0}, y^{*}_{1}, ..., y^{*}_{T})
  double lambdahat = phihatf(0);
  double alphahat = phihatf(1);
  double betahat = phihatf(2);
  arma::vec gammahat = phihatf.subvec(3, q+2);
  arma::mat Shat = fc_S(lambdahat, N, W1);
  arma::mat varepsilonast_mat = fc_asmat(varepsilonast, N, T);
  // Y^{*}
  arma::mat Yast_mat(N, T+1); Yast_mat.fill(0.0);
  Yast_mat.col(0) = y0;
  for(int t = 1; t < (T+1); t++) {
    double muthat = muhat_T(t-1); // \widehat{\mu}_{t}
    arma::vec yasttminus1 = Yast_mat.col(t-1); // y^{*}_{t-1}
    arma::mat Zt = Zarray.slice(t-1); // Z_{t}
    arma::vec varepsilonastt = varepsilonast_mat.col(t-1); // \varepsilon^{*}_{t}
    Yast_mat.col(t) = Shat.i() * (thetahat + muthat + alphahat * yasttminus1 + betahat * W2 * yasttminus1 + Zt * gammahat + Bhat.i() * varepsilonastt); // y^{*}_{t}
  }
  return Yast_mat;
}

