# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

fc_estbth_BASErhohat <- function(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini) {
    .Call(`_QDSPM_fc_estbth_BASErhohat`, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
}

fc_estBnum_BASErhohat <- function(Bnum, dimALL, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini) {
    .Call(`_QDSPM_fc_estBnum_BASErhohat`, Bnum, dimALL, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
}

fc_estbth_BASErhocheck <- function(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini) {
    .Call(`_QDSPM_fc_estbth_BASErhocheck`, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
}

fc_estBnum_BASErhocheck <- function(Bnum, dimALL, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini) {
    .Call(`_QDSPM_fc_estBnum_BASErhocheck`, Bnum, dimALL, Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
}

fc_asmat <- function(vec1, nrow, ncol) {
    .Call(`_QDSPM_fc_asmat`, vec1, nrow, ncol)
}

fc_asvec <- function(mat1) {
    .Call(`_QDSPM_fc_asvec`, mat1)
}

fc_iota <- function(n) {
    .Call(`_QDSPM_fc_iota`, n)
}

fc_I <- function(n) {
    .Call(`_QDSPM_fc_I`, n)
}

fc_B <- function(rho, N, W3) {
    .Call(`_QDSPM_fc_B`, rho, N, W3)
}

fc_S <- function(lambda, N, W1) {
    .Call(`_QDSPM_fc_S`, lambda, N, W1)
}

fc_MatrixMultip <- function(mat1, mat2) {
    .Call(`_QDSPM_fc_MatrixMultip`, mat1, mat2)
}

fc_Zbar <- function(Z, N, T, q) {
    .Call(`_QDSPM_fc_Zbar`, Z, N, T, q)
}

fc_ITWZ <- function(N, q, Z, W) {
    .Call(`_QDSPM_fc_ITWZ`, N, q, Z, W)
}

fc_phitilde <- function(N, T, mathbbY, mathbbZbreve, Q) {
    .Call(`_QDSPM_fc_phitilde`, N, T, mathbbY, mathbbZbreve, Q)
}

fc_etaest <- function(phiest, Y, Zbreve) {
    .Call(`_QDSPM_fc_etaest`, phiest, Y, Zbreve)
}

fc_thetaest <- function(etaest, N, T) {
    .Call(`_QDSPM_fc_thetaest`, etaest, N, T)
}

fc_etaestMINUSiotaTthetaest <- function(etaest, thetaest, T) {
    .Call(`_QDSPM_fc_etaestMINUSiotaTthetaest`, etaest, thetaest, T)
}

fc_muTest <- function(etaest_minus_iotaTthetaest, N, T) {
    .Call(`_QDSPM_fc_muTest`, etaest_minus_iotaTthetaest, N, T)
}

fc_secA1 <- function(W3, N) {
    .Call(`_QDSPM_fc_secA1`, W3, N)
}

fc_Lambda1lambda1 <- function(Vtilde, N, T, W3, secA1, secA2) {
    .Call(`_QDSPM_fc_Lambda1lambda1`, Vtilde, N, T, W3, secA1, secA2)
}

fc_varepsilonest <- function(Best, Vest, N, T) {
    .Call(`_QDSPM_fc_varepsilonest`, Best, Vest, N, T)
}

fc_sigmavarepsilonsquareest <- function(varepsilonest, N, T) {
    .Call(`_QDSPM_fc_sigmavarepsilonsquareest`, varepsilonest, N, T)
}

fc_Sigmagest <- function(sigmavarepsilonsquareest, secA1, secA2, N) {
    .Call(`_QDSPM_fc_Sigmagest`, sigmavarepsilonsquareest, secA1, secA2, N)
}

fc_kappacheck <- function(Lambda1, lambda1) {
    .Call(`_QDSPM_fc_kappacheck`, Lambda1, lambda1)
}

fc_phihatf <- function(N, T, q, Y, Zbreve, thetatilde, mutilde_NT, Bhat, H) {
    .Call(`_QDSPM_fc_phihatf`, N, T, q, Y, Zbreve, thetatilde, mutilde_NT, Bhat, H)
}

fc_varepsilonast <- function(east, varepsilonhat) {
    .Call(`_QDSPM_fc_varepsilonast`, east, varepsilonhat)
}

f_Yastmat_Generation <- function(varepsilonast, N, T, q, y0, phihatf, thetahat, muhat_T, Zarray, W1, W2, Bhat) {
    .Call(`_QDSPM_f_Yastmat_Generation`, varepsilonast, N, T, q, y0, phihatf, thetahat, muhat_T, Zarray, W1, W2, Bhat)
}

fc_Rrho <- function(Lambda1, lambda1, aTa, kappa) {
    .Call(`_QDSPM_fc_Rrho`, Lambda1, lambda1, aTa, kappa)
}

