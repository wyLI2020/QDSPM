############################## Two-stage Hybrid Bootstrap ##############################

fr_estbth_BASErhohat <- function(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bhatfo_upd, phihatf, thetahat, muhat_T, Vhat, rho_ini) {
  ## estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  q <- ncol(Z)
  p <- ncol(Xbreve) - q
  ### preparation
  #### random weights {e^{*}_{it}, i=1,...,N, t=1,...,T} in Stage 1, i.i.d. with zero mean, one variance and finite fourth moment
  if (Me == 1) { # If Me=1, we use the standard normal distribution weights
    east <- rnorm(N*T)
  } else if (Me == 2) { # If Me=2, we use the two-point distribution, takes the value of -1 or 1 with probability 0.5
    east <- rbinom(N*T, 1, 0.5)*2 - 1
  }
  #### random weights {\pi_{i}, i=1,...,N} in Stage 2, i.i.d. non-negative with mean and variance both equal to 1
  if (Mpi == 1) { # If Mpi=1, we use the standard exponential distribution weights
    piast <- rexp(N)
  } else if (Mpi == 2) { # If Mpi=2, we use the Rademacher distribution, which takes the value of 0 or 2 with probability 0.5
    piast <- rbinom(N, 1, 0.5)*2
  }
  #### transform Z to Zarray
  Zarraymid <- array(as.vector(t(Z)), dim = c(q, N, T))
  Zarray <- aperm(Zarraymid,c(2,1,3)) # dim = c(N, q, T)
  ############################## given \widehat{\phi}_{f}, \widehat{\rho}_{fo, upd}, \widehat{\theta}, \widehat{\mu} and \widehat{V} ##############################
  ### Stage 1, , the recursive-design wild bootstrap
  #### bootstrap samples
  varepsilonhat <- fc_varepsilonest(Bhatfo_upd, Vhat, N, T) # \widehat{\varepsilon} = (I_{T} \otimes B(\widehat{\rho}_{fo})) \widehat{V}
  varepsilonast_h <- fc_varepsilonast(east, varepsilonhat) # \varepsilon^{*}
  Yast_h_mat <- f_Yastmat_Generation(varepsilonast_h, N, T, q, y0, phihatf, thetahat, muhat_T, Zarray, W1, W2, Bhatfo_upd) # (y^{*}_{0}, y^{*}_{1}, ..., y^{*}_{T}) 
  Yast_h <- as.vector(Yast_h_mat[,2:(T+1)]) # Y^{*}
  Yastlag1_h <- as.vector(Yast_h_mat[,1:T]) # Y^{*}_{-1}
  Zbreveast_h <- cbind(as.vector(fc_MatrixMultip(W1, matrix(Yast_h, N, T))), Yastlag1_h, as.vector(fc_MatrixMultip(W2, matrix(Yastlag1_h, N, T))), Z) # \breve{Z}^{*} = ((I_{T} \otimes W_{1}) Y^{*}, Y^{*}_{-1}, (I_{T} \otimes W_{2}) Y^{*}_{-1}, Z)
  #### bootstrap parameter estimates
  ##### Pre-estimation
  phitildeB_h <- fr_phitilde(N, T, q, Yast_h, Yastlag1_h, Z, W1, W2) # \widetilde{\phi}_{B} 
  etatildeB_h <- fc_etaest(phitildeB_h, Yast_h, Zbreveast_h) # \widetilde{\eta}_{B}
  thetatildeB_h <- fc_thetaest(etatildeB_h, N, T) # \widetilde{\theta}_{B}
  etatildeB_minus_iotaTthetatildeB_h <- fc_etaestMINUSiotaTthetaest(etatildeB_h, thetatildeB_h, T) # \widetilde{\eta}_{B} - \iota_{T} \otimes \widetilde{\theta}_{B}
  mutildeB_T_h <- fc_muTest(etatildeB_minus_iotaTthetatildeB_h, N, T) # (\widetilde{\mu}_{1B}, ..., \widetilde{\mu}_{TB})^{\prime}
  mutildeB_NT_h <- rep(mutildeB_T_h, each=N) # (\widetilde{\mu}_{1B} \iota_{N}^{\prime}, \ldots, \widetilde{\mu}_{TB} \iota_{N}^{\prime})^{\prime}
  VtildeB_h <- etatildeB_minus_iotaTthetatildeB_h - mutildeB_NT_h # \widetilde{V}_{B}
  #### Stage 1
  ##### \widehat{\rho}_{foB} based on \widetilde{V}_{B}
  LambdaBlambdaB_h <- fc_Lambda1lambda1(VtildeB_h, N, T, W3, secA1, secA2)
  LambdaB_h <- LambdaBlambdaB_h[1:2,1:2] # \Lambda_{B}
  lambdaB_h <- LambdaBlambdaB_h[,3] # \lambda_{B}
  rhohatfoB <- fr_rhohatfo(N, T, W3, secA1, secA2, LambdaB_h, lambdaB_h, VtildeB_h, rho_ini) # \widehat{\rho}_{foB}
  ##### \widehat{\phi}_{fB} based on \widehat{\rho}_{foB}, \widetilde{\theta}_{B} and \widetilde{\mu}_{B}
  BhatfoB <- fc_B(rhohatfoB, N, W3) # B(\widehat{\rho}_{foB})
  phihatfB <- fr_phihatf(N, T, q, Yast_h, Yastlag1_h, Zbreveast_h, thetatildeB_h, mutildeB_NT_h, BhatfoB, Z, W3) # \widehat{\phi}_{fB}
  ##### \widehat{\theta}_{B}, \widehat{\mu}_{B} and \widehat{V}_{B} based on \widehat{\phi}_{fB}
  etahatB <- fc_etaest(phihatfB, Yast_h, Zbreveast_h) # \widehat{\eta}_{B}
  thetahatB <- fc_thetaest(etahatB, N, T) # \widehat{\theta}_{B}
  etahatB_minus_iotaTthetahatB <- fc_etaestMINUSiotaTthetaest(etahatB, thetahatB, T) # \widehat{\eta}_{B} - \iota_{T} \otimes \widehat{\theta}_{B}
  muhatB_T <- fc_muTest(etahatB_minus_iotaTthetahatB, N, T) # (\widehat{\mu}_{1B}, ..., \widehat{\mu}_{TB})^{\prime}
  muhatB_NT <- rep(muhatB_T, each=N) # (\widehat{\mu}_{1B} \iota_{N}^{\prime}, \ldots, \widehat{\mu}_{TB} \iota_{N}^{\prime})^{\prime}
  VhatB <- etahatB_minus_iotaTthetahatB - muhatB_NT # \widehat{V}_{B}
  ##### update \widehat{\rho}_{foB} based on \widehat{V}_{B}, denoted by \widehat{\rho}_{foB, upd}
  LambdaBlambdaB_h_upd <- fc_Lambda1lambda1(VhatB, N, T, W3, secA1, secA2)
  LambdaB_h_upd <- LambdaBlambdaB_h_upd[1:2,1:2] # \Lambda_{B, upd}
  lambdaB_h_upd <- LambdaBlambdaB_h_upd[,3] # \lambda_{B, upd}
  rhohatfoB_upd <- fr_rhohatfo(N, T, W3, secA1, secA2, LambdaB_h_upd, lambdaB_h_upd, VhatB, rho_ini) # \widehat{\rho}_{foB, upd}
  #### Stage 2, the random weight bootstrap
  varphihatB <- fr_varphihatB(tauseq, piast, thetahatB, Xbreve, p, q) # \widehat{\varphi}_{B}(\tau) for \tau in tauseq
  ### results
  estB <- c(phitildeB_h, rhohatfoB, phihatfB, rhohatfoB_upd, as.vector(varphihatB))
  return(estB)
}

fr_estbth_BASErhocheck <- function(Me, Mpi, tauseq, N, T, y0, Z, Xbreve, W1, W2, W3, secA1, secA2, Bcheck_upd, phicheckf, thetacheck, mucheck_T, Vcheck, rho_ini) {
  ## estimation in the b-th step of the two-stage hybrid bootstrapping procedure
  q <- ncol(Z)
  p <- ncol(Xbreve) - q
  ### preparation
  #### random weights {e^{*}_{it}, i=1,...,N, t=1,...,T} in Stage 1, i.i.d. with zero mean, one variance and finite fourth moment
  if (Me == 1) { # If Me=1, we use the standard normal distribution weights
    east <- rnorm(N*T)
  } else if (Me == 2) { # If Me=2, we use the two-point distribution, takes the value of -1 or 1 with probability 0.5
    east <- rbinom(N*T, 1, 0.5)*2 - 1
  }
  #### random weights {\pi_{i}, i=1,...,N} in Stage 2, i.i.d. non-negative with mean and variance both equal to 1
  if (Mpi == 1) { # If Mpi=1, we use the standard exponential distribution weights
    piast <- rexp(N)
  } else if (Mpi == 2) { # If Mpi=2, we use the Rademacher distribution, which takes the value of 0 or 2 with probability 0.5
    piast <- rbinom(N, 1, 0.5)*2
  }
  #### transform Z to Zarray
  Zarraymid <- array(as.vector(t(Z)), dim = c(q, N, T))
  Zarray <- aperm(Zarraymid,c(2,1,3)) # dim = c(N, q, T)
  ############################## given \check{\phi}_{f}, \check{\rho}_{upd}, \check{\theta}, \check{\mu} and \check{V} ##############################
  ### Stage 1, , the recursive-design wild bootstrap
  #### bootstrap samples
  varepsiloncheck <- fc_varepsilonest(Bcheck_upd, Vcheck, N, T) # \check{\varepsilon} = (I_{T} \otimes B(\check{\rho})) \check{V}
  varepsilonast_c <- fc_varepsilonast(east, varepsiloncheck) # \varepsilon^{*}
  Yast_c_mat <- f_Yastmat_Generation(varepsilonast_c, N, T, q, y0, phicheckf, thetacheck, mucheck_T, Zarray, W1, W2, Bcheck_upd) # (y^{*}_{0}, y^{*}_{1}, ..., y^{*}_{T}) 
  Yast_c <- as.vector(Yast_c_mat[,2:(T+1)]) # Y^{*}
  Yastlag1_c <- as.vector(Yast_c_mat[,1:T]) # Y^{*}_{-1}
  Zbreveast_c <- cbind(as.vector(fc_MatrixMultip(W1, matrix(Yast_c, N, T))), Yastlag1_c, as.vector(fc_MatrixMultip(W2, matrix(Yastlag1_c, N, T))), Z) # \breve{Z}^{*} = ((I_{T} \otimes W_{1}) Y^{*}, Y^{*}_{-1}, (I_{T} \otimes W_{2}) Y^{*}_{-1}, Z)
  #### bootstrap parameter estimates
  ##### Pre-estimation
  phitildeB_c <- fr_phitilde(N, T, q, Yast_c, Yastlag1_c, Z, W1, W2) # \widetilde{\phi}_{B} 
  etatildeB_c <- fc_etaest(phitildeB_c, Yast_c, Zbreveast_c) # \widetilde{\eta}_{B}
  thetatildeB_c <- fc_thetaest(etatildeB_c, N, T) # \widetilde{\theta}_{B}
  etatildeB_minus_iotaTthetatildeB_c <- fc_etaestMINUSiotaTthetaest(etatildeB_c, thetatildeB_c, T) # \widetilde{\eta}_{B} - \iota_{T} \otimes \widetilde{\theta}_{B}
  mutildeB_T_c <- fc_muTest(etatildeB_minus_iotaTthetatildeB_c, N, T) # (\widetilde{\mu}_{1B}, ..., \widetilde{\mu}_{TB})^{\prime}
  mutildeB_NT_c <- rep(mutildeB_T_c, each=N) # (\widetilde{\mu}_{1B} \iota_{N}^{\prime}, \ldots, \widetilde{\mu}_{TB} \iota_{N}^{\prime})^{\prime}
  VtildeB_c <- etatildeB_minus_iotaTthetatildeB_c - mutildeB_NT_c # \widetilde{V}_{B}
  #### Stage 1
  ##### \check{\rho}_{B} based on \widetilde{V}_{B}
  LambdaBlambdaB_c <- fc_Lambda1lambda1(VtildeB_c, N, T, W3, secA1, secA2)
  LambdaB_c <- LambdaBlambdaB_c[1:2,1:2] # \Lambda_{B}
  lambdaB_c <- LambdaBlambdaB_c[,3] # \lambda_{B}
  rhocheckB <- fr_rhocheck(LambdaB_c, lambdaB_c) # \check{\rho}_{B}
  ##### \check{\phi}_{fB} based on \check{\rho}_{B}, \widetilde{\theta}_{B} and \widetilde{\mu}_{B}
  BcheckB <- fc_B(rhocheckB, N, W3) # B(\check{\rho}_{B})
  phicheckfB <- fr_phihatf(N, T, q, Yast_c, Yastlag1_c, Zbreveast_c, thetatildeB_c, mutildeB_NT_c, BcheckB, Z, W3) # \check{\phi}_{fB}
  ##### \check{\theta}_{B}, \check{\mu}_{B} and \check{V}_{B} based on \check{\phi}_{fB}
  etacheckB <- fc_etaest(phicheckfB, Yast_c, Zbreveast_c) # \check{\eta}_{B}
  thetacheckB <- fc_thetaest(etacheckB, N, T) # \check{\theta}_{B}
  etacheckB_minus_iotaTthetacheckB <- fc_etaestMINUSiotaTthetaest(etacheckB, thetacheckB, T) # \check{\eta}_{B} - \iota_{T} \otimes \check{\theta}_{B}
  mucheckB_T <- fc_muTest(etacheckB_minus_iotaTthetacheckB, N, T) # (\check{\mu}_{1B}, ..., \check{\mu}_{TB})^{\prime}
  mucheckB_NT <- rep(mucheckB_T, each=N) # (\check{\mu}_{1B} \iota_{N}^{\prime}, \ldots, \check{\mu}_{TB} \iota_{N}^{\prime})^{\prime}
  VcheckB <- etacheckB_minus_iotaTthetacheckB - mucheckB_NT # \check{V}_{B}
  ##### update \check{\rho}_{B} based on \check{V}_{B}, denoted by \check{\rho}_{B, upd}
  LambdaBlambdaB_c_upd <- fc_Lambda1lambda1(VcheckB, N, T, W3, secA1, secA2)
  LambdaB_c_upd <- LambdaBlambdaB_c_upd[1:2,1:2] # \Lambda_{B, upd}
  lambdaB_c_upd <- LambdaBlambdaB_c_upd[,3] # \lambda_{B, upd}
  rhocheckB_upd <- fr_rhocheck(LambdaB_c_upd, lambdaB_c_upd) # \check{\rho}_{B, upd}
  #### Stage 2, the random weight bootstrap
  varphicheckB <- fr_varphihatB(tauseq, piast, thetacheckB, Xbreve, p, q) # \check{\varphi}_{B}(\tau) for \tau in tauseq
  ### results
  estB <- c(phitildeB_c, rhocheckB, phicheckfB, rhocheckB_upd, as.vector(varphicheckB))
  return(estB)
}

