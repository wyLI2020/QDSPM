########################################### Estimation ###########################################
fr_estDSPM_ASDBTS4M_BASErhohat <- function(tauseq, N, T, q, p, Y, Ylag1, Z, X, W1, W2, W3, rho_ini, Bnum) {
  ## Estimation of \rho, \phi and \varphi(\tau), that is \widehat{\rho}, \check{\rho}, \widetilde{\phi}, \widetilde{\phi}_{f}, and \widehat{\varphi}(\tau)
  # N <- nrow(X)
  # T <- nrow(Z) / N
  # q <- ncol(Z)
  # p <- ncol(X)
  Zbreve <- cbind(as.vector(fc_MatrixMultip(W1, matrix(Y, N, T))), Ylag1, as.vector(fc_MatrixMultip(W2, matrix(Ylag1, N, T))), Z) # \breve{Z} = ((I_{T} \otimes W_{1}) Y, Y_{-1}, (I_{T} \otimes W_{2}) Y_{-1}, Z)
  Xbreve <- cbind(X, fc_Zbar(Z, N, T, q))
  ############################## Estimation Procedure ##############################
  ### Pre-estimation
  phitilde <- fr_phitilde(N, T, q, Y, Ylag1, Z, W1, W2) # An initial estimator \widetilde{\phi} of \phi = (\lambda, \alpha, \beta, \gamma^{\prime})^{\prime} 
  etatilde <- fc_etaest(phitilde, Y, Zbreve) # An approximate of \eta, \widetilde{\eta}, using \widetilde{\phi}
  thetatilde <- fc_thetaest(etatilde, N, T) # An approximate of \theta, \widetilde{\theta}
  etatilde_minus_iotaTthetatilde <- fc_etaestMINUSiotaTthetaest(etatilde, thetatilde, T) # \widetilde{\eta} - \iota_{T} \otimes \widetilde{\theta}
  mutilde_T <- fc_muTest(etatilde_minus_iotaTthetatilde, N, T) # (\widetilde{\mu}_{1}, ..., \widetilde{\mu}_{T})^{\prime}
  mutilde_NT <- rep(mutilde_T, each=N) # An approximate of \mu = (\mu_{10} \iota_{N}^{\prime}, ..., \mu_{T0} \iota_{N}^{\prime})^{\prime}, (\widetilde{\mu}_{1} \iota_{N}^{\prime}, \ldots, \widetilde{\mu}_{T} \iota_{N}^{\prime})^{\prime}
  Vtilde <- etatilde_minus_iotaTthetatilde - mutilde_NT # An approximate of V, \widetilde{V}
  ### Stage 1
  #### \widehat{\rho}_{fo} based on \widetilde{V}
  secA1 <- fc_secA1(W3, N) # A_{1} \triangleq (I_{T} \otimes secA1)
  secA2 <- W3 # A_{2} = M = I_{T} \otime W_{3}
  Lambda1lambda1 <- fc_Lambda1lambda1(Vtilde, N, T, W3, secA1, secA2)
  Lambda1 <- Lambda1lambda1[1:2,1:2] # \Lambda_{NT}
  lambda1 <- Lambda1lambda1[,3] # \lambda_{NT}
  rhohatfo <- fr_rhohatfo(N, T, W3, secA1, secA2, Lambda1, lambda1, Vtilde, rho_ini) # The feasible optimal GMM (FOGMM) estimator of \rho, \widehat{\rho}_{fo}, using \widetilde{V} and weighting matrix \widetilde{\Sigma}_{g}^{-1}
  #### \widehat{\phi}_{f} based on \widehat{\rho}_{fo}, \widetilde{\theta} and \widetilde{\mu}
  Bhatfo <- fc_B(rhohatfo, N, W3) # B(\widehat{\rho}_{fo})
  phihatf <- fr_phihatf(N, T, q, Y, Ylag1, Zbreve, thetatilde, mutilde_NT, Bhatfo, Z, W3) # The FGSP2SLS estimator of \phi, \widehat{\phi}_{f}, using \widehat{\rho}_{fo}
  #### \widehat{\theta}, \widehat{\mu} and \widehat{V} based on \widehat{\phi}_{f}
  etahat <- fc_etaest(phihatf, Y, Zbreve) # An approximate of \eta, \widehat{\eta}, using \widehat{\phi}_{f}
  thetahat <- fc_thetaest(etahat, N, T) # An approximate of \theta, \widehat{\theta}
  etahat_minus_iotaTthetahat <- fc_etaestMINUSiotaTthetaest(etahat, thetahat, T) # \widehat{\eta} - \iota_{T} \otimes \widehat{\theta}
  muhat_T <- fc_muTest(etahat_minus_iotaTthetahat, N, T) # (\widehat{\mu}_{1}, ..., \widehat{\mu}_{T})^{\prime}
  muhat_NT <- rep(muhat_T, each=N) # An approximate of \mu = (\mu_{10} \iota_{N}^{\prime}, ..., \mu_{T0} \iota_{N}^{\prime})^{\prime}, (\widehat{\mu}_{1} \iota_{N}^{\prime}, \ldots, \widehat{\mu}_{T} \iota_{N}^{\prime})^{\prime}
  Vhat <- etahat_minus_iotaTthetahat - muhat_NT # An approximate of V, \widehat{V}
  #### update \widehat{\rho}_{fo} based on \widehat{V}, denoted by \widehat{\rho}_{fo, upd}
  Lambda1lambda1_upd <- fc_Lambda1lambda1(Vhat, N, T, W3, secA1, secA2)
  Lambda1_upd <- Lambda1lambda1_upd[1:2,1:2] # updated \Lambda_{NT}
  lambda1_upd <- Lambda1lambda1_upd[,3] # updated \lambda_{NT}
  rhohatfo_upd <- fr_rhohatfo(N, T, W3, secA1, secA2, Lambda1_upd, lambda1_upd, Vhat, rho_ini) # updated \widehat{\rho}_{fo} based on \widehat{\phi}_{f}
  Bhatfo_upd <- fc_B(rhohatfo_upd, N, W3) # B(\widehat{\rho}_{fo, upd})
  ### Stage 2
  #### \widehat{\varphi}(\tau) based on \widehat{\theta}
  varphihat <- fr_varphihat(tauseq, thetahat, Xbreve, p, q) # \widehat{\varphi}(\tau) for \tau in tauseq
  ### estimates
  estALL <- c(phitilde, rhohatfo, phihatf, rhohatfo_upd, as.vector(varphihat))
  dimALL <- length(estALL)
  ############################## Two-stage Hybrid Bootstrap ##############################
  ### ASD by the two-stage hybrid bootstrapping procedure, given \widehat{\phi}_{f}, \widehat{\rho}_{fo, upd}, \widehat{\theta}, \widehat{\mu} and \widehat{V}
  y0 <- Ylag1[1:N]
  estALL_BTSM11_Bnum <- fc_estBnum_BASErhohat(Bnum, dimALL, Me=1, Mpi=1, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
  estALL_BTSM12_Bnum <- fc_estBnum_BASErhohat(Bnum, dimALL, Me=1, Mpi=2, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
  estALL_BTSM21_Bnum <- fc_estBnum_BASErhohat(Bnum, dimALL, Me=2, Mpi=1, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
  estALL_BTSM22_Bnum <- fc_estBnum_BASErhohat(Bnum, dimALL, Me=2, Mpi=2, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini)
  asdALL_BTSM11 <- apply(estALL_BTSM11_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM12 <- apply(estALL_BTSM12_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM21 <- apply(estALL_BTSM21_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM22 <- apply(estALL_BTSM22_Bnum, MARGIN = 1, FUN = sd)
  ############################## Results ##############################
  ### results
  results <- c(estALL, asdALL_BTSM11, asdALL_BTSM12, asdALL_BTSM21, asdALL_BTSM22)
  results_mat <- matrix(results, nrow = dimALL, ncol = 5)
  namephitilde <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "tilde")
  namephihatf <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "hatf")
  namephicheckf <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "checkf")
  namevarphihat <- paste0("varphihat", rep(0:(p+q), length(tauseq)), "_tau", rep(tauseq, each = p+q+1))
  namevarphicheck <- paste0("varphicheck", rep(0:(p+q), length(tauseq)), "_tau", rep(tauseq, each = p+q+1))
  rownames(results_mat) <- c(namephitilde, "rhohatfo", namephihatf, "rhohatfo_upd", namevarphihat)
  colnames(results_mat) <- c("est", "ASD_BTSM11", "ASD_BTSM12", "ASD_BTSM21", "ASD_BTSM22")
  results_list <- list(EST=results_mat, muhat=muhat_T)
  return(results_list)
}

fr_estDSPM_ASDBTS4M_BASErhocheck <- function(tauseq, N, T, q, p, Y, Ylag1, Z, X, W1, W2, W3, rho_ini, Bnum) {
  ## Estimation of \rho, \phi and \varphi(\tau), that is \widehat{\rho}, \check{\rho}, \widetilde{\phi}, \widetilde{\phi}_{f}, and \widehat{\varphi}(\tau)
  # N <- nrow(X)
  # T <- nrow(Z) / N
  # q <- ncol(Z)
  # p <- ncol(X)
  Zbreve <- cbind(as.vector(fc_MatrixMultip(W1, matrix(Y, N, T))), Ylag1, as.vector(fc_MatrixMultip(W2, matrix(Ylag1, N, T))), Z) # \breve{Z} = ((I_{T} \otimes W_{1}) Y, Y_{-1}, (I_{T} \otimes W_{2}) Y_{-1}, Z)
  Xbreve <- cbind(X, fc_Zbar(Z, N, T, q))
  ############################## Estimation Procedure ##############################
  ### Pre-estimation
  phitilde <- fr_phitilde(N, T, q, Y, Ylag1, Z, W1, W2) # An initial estimator \widetilde{\phi} of \phi = (\lambda, \alpha, \beta, \gamma^{\prime})^{\prime} 
  etatilde <- fc_etaest(phitilde, Y, Zbreve) # An approximate of \eta, \widetilde{\eta}, using \widetilde{\phi}
  thetatilde <- fc_thetaest(etatilde, N, T) # An approximate of \theta, \widetilde{\theta}
  etatilde_minus_iotaTthetatilde <- fc_etaestMINUSiotaTthetaest(etatilde, thetatilde, T) # \widetilde{\eta} - \iota_{T} \otimes \widetilde{\theta}
  mutilde_T <- fc_muTest(etatilde_minus_iotaTthetatilde, N, T) # (\widetilde{\mu}_{1}, ..., \widetilde{\mu}_{T})^{\prime}
  mutilde_NT <- rep(mutilde_T, each=N) # An approximate of \mu = (\mu_{10} \iota_{N}^{\prime}, ..., \mu_{T0} \iota_{N}^{\prime})^{\prime}, (\widetilde{\mu}_{1} \iota_{N}^{\prime}, \ldots, \widetilde{\mu}_{T} \iota_{N}^{\prime})^{\prime}
  Vtilde <- etatilde_minus_iotaTthetatilde - mutilde_NT # An approximate of V, \widetilde{V}
  ### Stage 1
  #### \check{\rho} based on \widetilde{V}
  secA1 <- fc_secA1(W3, N) # A_{1} \triangleq (I_{T} \otimes secA1)
  secA2 <- W3 # A_{2} = M = I_{T} \otime W_{3}
  Lambda1lambda1 <- fc_Lambda1lambda1(Vtilde, N, T, W3, secA1, secA2)
  Lambda1 <- Lambda1lambda1[1:2,1:2] # \Lambda_{NT}
  lambda1 <- Lambda1lambda1[,3] # \lambda_{NT}
  rhocheck <- fr_rhocheck(Lambda1, lambda1) # The linear least square estimator of \rho, \check{\rho}
  #### \check{\phi}_{f} based on \check{\rho}, \widetilde{\theta} and \widetilde{\mu}
  Bcheck <- fc_B(rhocheck, N, W3) # B(\check{\rho})
  phicheckf <- fr_phihatf(N, T, q, Y, Ylag1, Zbreve, thetatilde, mutilde_NT, Bcheck, Z, W3) # The FGSP2SLS estimator of \phi, \check{\phi}_{f}, using \check{\rho}
  #### \check{\theta}, \check{\mu} and \check{V} based on \check{\phi}_{f}
  etacheck <- fc_etaest(phicheckf, Y, Zbreve) # An approximate of \eta, \check{\eta}, using \check{\phi}_{f}
  thetacheck <- fc_thetaest(etacheck, N, T) # An approximate of \theta, \check{\theta}
  etacheck_minus_iotaTthetacheck <- fc_etaestMINUSiotaTthetaest(etacheck, thetacheck, T) # \check{\eta} - \iota_{T} \otimes \check{\theta}
  mucheck_T <- fc_muTest(etacheck_minus_iotaTthetacheck, N, T) # (\check{\mu}_{1}, ..., \check{\mu}_{T})^{\prime}
  mucheck_NT <- rep(mucheck_T, each=N) # An approximate of \mu = (\mu_{10} \iota_{N}^{\prime}, ..., \mu_{T0} \iota_{N}^{\prime})^{\prime}, (\check{\mu}_{1} \iota_{N}^{\prime}, \ldots, \check{\mu}_{T} \iota_{N}^{\prime})^{\prime}
  Vcheck <- etacheck_minus_iotaTthetacheck - mucheck_NT # An approximate of V, \check{V}
  #### update \check{\rho} based on \check{V}, denoted by \check{\rho}_{upd}
  Lambda1lambda1_upd <- fc_Lambda1lambda1(Vcheck, N, T, W3, secA1, secA2)
  Lambda1_upd <- Lambda1lambda1_upd[1:2,1:2] # updated \Lambda_{NT}
  lambda1_upd <- Lambda1lambda1_upd[,3] # updated \lambda_{NT}
  rhocheck_upd <- fr_rhocheck(Lambda1_upd, lambda1_upd) # updated \check{\rho} based on \check{\phi}_{f}
  Bcheck_upd <- fc_B(rhocheck_upd, N, W3) # B(\check{\rho}_{upd})
  ### Stage 2
  #### \check{\varphi}{\tau} based on \check{\theta}
  varphicheck <- fr_varphihat(tauseq, thetacheck, Xbreve, p, q) # \check{\varphi}(\tau) for \tau in tauseq
  ### estimates
  estALL <- c(phitilde, rhocheck, phicheckf, rhocheck_upd, as.vector(varphicheck))
  dimALL <- length(estALL)
  ############################## Two-stage Hybrid Bootstrap ##############################
  ### ASD by the two-stage hybrid bootstrapping procedure, given \check{\phi}_{f}, \check{\rho}_{upd}, \check{\theta}, \check{\mu} and \check{V}
  y0 <- Ylag1[1:N]
  estALL_BTSM11_Bnum <- fc_estBnum_BASErhocheck(Bnum, dimALL, Me=1, Mpi=1, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
  estALL_BTSM12_Bnum <- fc_estBnum_BASErhocheck(Bnum, dimALL, Me=1, Mpi=2, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
  estALL_BTSM21_Bnum <- fc_estBnum_BASErhocheck(Bnum, dimALL, Me=2, Mpi=1, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
  estALL_BTSM22_Bnum <- fc_estBnum_BASErhocheck(Bnum, dimALL, Me=2, Mpi=2, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini)
  asdALL_BTSM11 <- apply(estALL_BTSM11_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM12 <- apply(estALL_BTSM12_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM21 <- apply(estALL_BTSM21_Bnum, MARGIN = 1, FUN = sd)
  asdALL_BTSM22 <- apply(estALL_BTSM22_Bnum, MARGIN = 1, FUN = sd)
  ############################## Results ##############################
  ### results
  results <- c(estALL, asdALL_BTSM11, asdALL_BTSM12, asdALL_BTSM21, asdALL_BTSM22)
  results_mat <- matrix(results, nrow = dimALL, ncol = 5)
  namephitilde <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "tilde")
  namephihatf <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "hatf")
  namephicheckf <- paste0(c("lambda", "alpha", "beta", paste0("gamma", 1:q)), "checkf")
  namevarphihat <- paste0("varphihat", rep(0:(p+q), length(tauseq)), "_tau", rep(tauseq, each = p+q+1))
  namevarphicheck <- paste0("varphicheck", rep(0:(p+q), length(tauseq)), "_tau", rep(tauseq, each = p+q+1))
  rownames(results_mat) <- c(namephitilde, "rhocheck", namephicheckf, "rhocheck_upd", namevarphicheck)
  colnames(results_mat) <- c("est", "ASD_BTSM11", "ASD_BTSM12", "ASD_BTSM21", "ASD_BTSM22")
  results_list <- list(EST=results_mat, mucheck=mucheck_T)
  return(results_list)
}

