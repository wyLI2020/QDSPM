############################## Estimations ##############################

fr_phitilde <- function(N, T, q, Y, Ylag1, Z, W1, W2) {
  # intial estimator \widetilde{\phi} 
  wrY <- Y[(2*N+1):(N*T)] # Y^{\wr} = (y_{3}^{\prime}, ..., y_{T}^{\prime})^{\prime}
  wrYlag1 <- Y[(N+1):(N*(T-1))] # Y^{\wr}_{-1} = (y_{2}^{\prime}, ..., y_{T-1}^{\prime})^{\prime}
  wrYlag2 <- Y[1:(N*(T-2))] # Y^{\wr}_{-2} = (y_{1}^{\prime}, ..., y_{T-2}^{\prime})^{\prime}
  wrYlag3 <- Ylag1[1:(N*(T-2))] # Y^{\wr}_{-3} = (y_{0}^{\prime}, ..., y_{T-3}^{\prime})^{\prime}
  wrZ <- as.matrix(Z[(2*N+1):(N*T),]) # Z^{\wr} = (Z_{3}^{\prime}, ..., Z_{T}^{\prime})^{\prime}
  wrZlag1 <- as.matrix(Z[(N+1):(N*(T-1)),]) # Z^{\wr}_{-1} = (Z_{2}^{\prime}, ..., Z_{T-1}^{\prime})^{\prime}
  wrZlag2 <- as.matrix(Z[1:(N*(T-2)),]) # Z^{\wr}_{-2} = (Z_{1}^{\prime}, ..., Z_{T-2}^{\prime})^{\prime}
  mathbbY <- wrY - wrYlag1 # \mathbb{Y} = Y^{\wr} - Y^{\wr}_{-1}
  mathbbYlag1 <- wrYlag1 - wrYlag2 # \mathbb{Y}_{-1} = Y^{\wr}_{-1} - Y^{\wr}_{-2}
  mathbbYlag2 <- wrYlag2 - wrYlag3 # \mathbb{Y}_{-2} = Y^{\wr}_{-2} - Y^{\wr}_{-3}
  mathbbZ <- wrZ - wrZlag1 # \mathbb{Z} = Z^{\wr} - Z^{\wr}_{-1}
  mathbbZlag1 <- wrZlag1 - wrZlag2 # \mathbb{Z} = Z^{\wr}_{-1} - Z^{\wr}_{-2}
  mathbbZbreve <- cbind(as.vector(fc_MatrixMultip(W1, matrix(mathbbY, N, T-2))), mathbbYlag1, as.vector(fc_MatrixMultip(W2, matrix(mathbbYlag1, N, T-2))), mathbbZ) # \mathbb{\breve{Z}} = ((\iota_{T-2} \otimes W_{1})\mathbb{Y}, \mathbb{Y}_{-1}, (\iota_{T-2} \otimes W_{2})\mathbb{Y}_{-1}, \mathbb{Z})
  W1square <- fc_MatrixMultip(W1, W1)
  Q <- cbind(mathbbZ, fc_ITWZ(N, q, mathbbZ, W1), mathbbZlag1, fc_ITWZ(N, q, mathbbZlag1, W1), mathbbYlag2, fc_ITWZ(N, 1, matrix(mathbbYlag2, ncol=1), W1), fc_ITWZ(N, 1, matrix(mathbbYlag2, ncol=1), W1square)) # instrument variables matrix Q_{N(T-2)} = (\mathbb{Z}, (I_{T-2} \otimes W_{1}) \mathbb{Z}, \mathbb{Z}_{-1}, (I_{T-2} \otimes W_{1}) \mathbb{Z}_{-1}, \mathbb{Y}_{-2}, (I_{T-2} \otimes W_{1}) \mathbb{Y}_{-2}, , (I_{T-2} \otimes W_{1}^{2}) \mathbb{Y}_{-2}) when W_{1}=W_{2}
  phitilde <- fc_phitilde(N, T, mathbbY, mathbbZbreve, Q) # \widetilde{\phi} 
  return(phitilde)
}

fr_rhohatfo <- function(N, T, W3, secA1, secA2, Lambda1, lambda1, Vtilde, rho_ini) {
  # \widehat{\rho}_{fo}
  rhohatI <- optim(par = rho_ini, fn = LossFUN_rho, method = "L-BFGS-B", Lambda1=Lambda1, lambda1=lambda1, aTa=diag(2))$par # The GMM estimator of \rho, \widehat{\rho}_{I} = \widehat{\rho}(I_{2}), using \widetilde{V} and weighting matrix I_{2}
  BhatI <- fc_B(rhohatI, N, W3) # B(\widehat{\rho}_{I})
  varepsilontilde <- fc_varepsilonest(BhatI, Vtilde, N, T) # \widetilde{\varepsilon} = (I_{T} \otimes B(\widehat{\rho}_{I})) \widetilde{V}
  sigmavarepsilonsquaretilde <- fc_sigmavarepsilonsquareest(varepsilontilde, N, T) # (\widetilde{\sigma}_{\varepsilon 1}^{2}, ..., \widetilde{\sigma}_{\varepsilon N}^{2})^{\prime}
  Sigma_gtilde <- fc_Sigmagest(sigmavarepsilonsquaretilde, secA1, secA2, N) # \widetilde{\Sigma}_{g}
  rhohatfo <- optim(par = rho_ini, fn = LossFUN_rho, method = "L-BFGS-B", Lambda1=Lambda1, lambda1=lambda1, aTa=solve(Sigma_gtilde))$par # The feasible optimal GMM (FOGMM) estimator of \rho, \widehat{\rho}_{fo}, using \widetilde{V} and weighting matrix \widetilde{\Sigma}_{g}^{-1}
  return(rhohatfo)
}

fr_rhocheck <- function(Lambda1, lambda1) {
  # \check{\rho}
  kappacheck <- fc_kappacheck(Lambda1, lambda1) # The linear least square estimator of \kappa = (\rho, \rho^{2})^{\prime}, \check{\kappa}, using \widetilde{V}
  rhocheck <- kappacheck[1] # \check{\rho}
  return(rhocheck)
}

fr_phihatf <- function(N, T, q, Y, Ylag1, Zbreve, thetatilde, mutilde_NT, Bhat, Z, W3) {
  # \widehat{\phi}_{f}
  W3square <- fc_MatrixMultip(W3, W3)
  H <- cbind(Ylag1, fc_ITWZ(N, 1, matrix(Ylag1, ncol=1), W3), fc_ITWZ(N, 1, matrix(Ylag1, ncol=1), W3square), Z, fc_ITWZ(N, q, Z, W3)) # instrument variables matrix H = (Y_{-1}, (I_{T} \otimes W_{3}) Y_{-1}, (I_{T} \otimes W_{3}^{2}) Y_{-1}, Z, (I_{T} \otimes W_{3}) Z) when W_{1}=W_{2}=W_{3}
  phihatf <- fc_phihatf(N, T, q, Y, Zbreve, thetatilde, mutilde_NT, Bhat, H)
  return(phihatf)
}

fr_varphihat <- function(tauseq, thetahat, Xbreve, p, q) {
  # \widehat{\varphi}(\tau) = (\widehat{\varphi}_{0}(\tau), \widehat{\underline{\varphi}}^{\prime}(\tau), \widehat{\overline{\varphi}}^{\prime}(\tau))^{\prime}, for \tau in tauseq
  varphihat <- matrix(nrow = 1+p+q, ncol = length(tauseq)) # \widehat{\varphi}(\tau) for \tau in tauseq
  if (length(tauseq) == 1) {
    varphihat[,1] <- rq(thetahat ~ Xbreve, tau = tauseq)$coefficients
  } else {
    varphihat <- rq(thetahat ~ Xbreve, tau = tauseq)$coefficients
  }
  return(varphihat)
}

############################## ASD by Two-stage Hybrid Bootstrap ##############################

fr_varphihatB <- function(tauseq, piast, thetahatB, Xbreve, p, q) {
  # \widehat{\varphi}_{B}(\tau) = (\widehat{\varphi}_{0B}(\tau), \widehat{\underline{\varphi}}_{B}^{\prime}(\tau), \widehat{\overline{\varphi}}_{B}^{\prime}(\tau))^{\prime} for \tau in tauseq, by random weights bootstrap
  varphihatB <- matrix(nrow = 1+p+q, ncol = length(tauseq)) # \widehat{\varphi}_{B}(\tau) for \tau in tauseq
  if (length(tauseq) == 1) {
    varphihatB[,1] <- rq(thetahatB ~ Xbreve, tau = tauseq, weights = piast)$coefficients
  } else {
    varphihatB <- rq(thetahatB ~ Xbreve, tau = tauseq, weights = piast)$coefficients
  }
  return(varphihatB)
}


