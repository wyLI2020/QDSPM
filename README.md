# QDSPM

## Installation

```R
#install.packages("devtools")
library(devtools)
install_github("wyLI2020/QDSPM")
```

## Usage

```R
fr_estDSPM_ASDBTS4M_BASErhohat(tauseq, N, T, q, p, Y, Ylag1, Z, X, W1, W2, W3, rho_ini, Bnum)
```

- **tauseq**:  vector, the quantile(s) to be estimated
- **N**: integer, the cross-sectional dimension
- **T**: integer, the time dimension
- **q**: integer, the number of time-varying explanatory variables
- **p**: integer, the number of time-invariant explanatory variables
- **Y**:  $NT$​-dimensional vector, response
- **Ylag1**:  $NT$-dimensional vector, the first lag of response
- **Z**:  $(NT, q)$ matrix, time-varying regressors
- **X**:  $(N, p)$ matrix, time-invariant regressors
- **W1**:  $(N, N)$ matrix, the spatial weights matrix
- **W2**:  $(N, N)$ matrix, the spatial weights matrix
- **W3**:  $(N, N)$ matrix, the spatial weights matrix
- **rho_ini**:  scalar, the initial value for parameter $\rho$ 
- **Bnum**:  integer, the number of bootstrap samples

## Example

```R
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
```

