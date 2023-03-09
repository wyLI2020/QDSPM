\name{QDSPM-Functions}
\alias{LossFUN_rho}
\alias{f_Yastmat_Generation}
\alias{fc_B}
\alias{fc_I}
\alias{fc_ITWZ}
\alias{fc_Lambda1lambda1}
\alias{fc_MatrixMultip}
\alias{fc_Rrho}
\alias{fc_S}
\alias{fc_Sigmagest}
\alias{fc_Zbar}
\alias{fc_asmat}
\alias{fc_asvec}
\alias{fc_estBnum_BASErhocheck}
\alias{fc_estBnum_BASErhohat}
\alias{fc_estbth_BASErhocheck}
\alias{fc_estbth_BASErhohat}
\alias{fc_etaest}
\alias{fc_etaestMINUSiotaTthetaest}
\alias{fc_iota}
\alias{fc_kappacheck}
\alias{fc_muTest}
\alias{fc_phihatf}
\alias{fc_phitilde}
\alias{fc_secA1}
\alias{fc_sigmavarepsilonsquareest}
\alias{fc_thetaest}
\alias{fc_varepsilonast}
\alias{fc_varepsilonest}
\alias{fr_estDSPM_ASDBTS4M_BASErhocheck}
\alias{fr_estDSPM_ASDBTS4M_BASErhohat}
\alias{fr_estbth_BASErhocheck}
\alias{fr_estbth_BASErhohat}
\alias{fr_phihatf}
\alias{fr_phitilde}
\alias{fr_rhocheck}
\alias{fr_rhohatfo}
\alias{fr_varphihat}
\alias{fr_varphihatB}
\title{Estimation and bootstrapping procedures}
\description{
  Estimation and bootstrapping procedures
}
\usage{
  fr_estDSPM_ASDBTS4M_BASErhohat(tauseq, N, T, q, p, Y, Ylag1, Z, X, W1, W2, W3, rho_ini, Bnum);
}
\arguments{
  \item{tauseq}{vector, the quantile(s) to be estimated}
  \item{N}{integer, the cross-sectional dimension}
  \item{T}{integer, the time dimension}
  \item{q}{integer, the number of time-varying explanatory variables}
  \item{p}{integer, the number of time-invariant explanatory variables}
  \item{Y}{\emph{NT}-dimensional vector, response}
  \item{Ylag1}{\emph{NT}-dimensional vector, the first lag of response}
  \item{Z}{(\emph{NT}, q) matrix, time-varying regressors}
  \item{X}{(\emph{N}, p) matrix, time-invariant regressors}
  \item{W1}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{W2}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{W3}{(\emph{N}, \emph{N}) matrix, the spatial weights matrix}
  \item{rho_ini}{scalar, the initial value for parameter}
  \item{Bnum}{integer, the number of bootstrap samples}
}
\value{
  \code{fr_estDSPM_ASDBTS4M_BASErhohat} returns a list.
}
\details{
  With the preproccess in following "Examples", we can get the estimators of the DSP model.
}
\examples{
  library(QDSPM)
  data("Dataset_AirQuality")
  W <- AirQualityData$W
  Y <- AirQualityData$Y
  Ylag1 <- AirQualityData$Ylag1
  Z <- AirQualityData$Z
  X <- AirQualityData$X
  N <- nrow(X); T <- length(Y) / N
  q <- ncol(Z); p <- ncol(X)
  set.seed(6)
  result_list_BASErhohat <- fr_estDSPM_ASDBTS4M_BASErhohat(tauseq=c(0.1,0.25,0.75,0.9), N, T, q, p, Y, Ylag1, Z, X, W, W, W, rho_ini=0.5, Bnum=1000)
  result_mat_BASErhohat <- result_list_BASErhohat$EST
  muhat_T <- result_list_BASErhohat$muhat
  result_BASErhohat <- matrix(NA, nrow = nrow(result_mat_BASErhohat), ncol = 4)
  result_BASErhohat[,c(1,2)] <- result_mat_BASErhohat[,c(1,2)]
  result_BASErhohat[,3] <- result_BASErhohat[,1] / result_BASErhohat[,2]
  result_BASErhohat[,4] <- 2*(1-pnorm(abs(result_BASErhohat[,1] / result_BASErhohat[,2])))
  rownames(result_BASErhohat) <- rownames(result_mat_BASErhohat); colnames(result_BASErhohat) <- c("est", "ASD", "z", "p-value")
  round(result_BASErhohat, 3)
}
% \references{
  % See the documentation for Armadillo, and RcppArmadillo, for more details.
% }
\author{Wenyu Li}

