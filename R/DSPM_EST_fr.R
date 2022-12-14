LossFUN_rho <- function(rho, Lambda1, lambda1, aTa) {
  # (\Lambda_{NT} \kappa - \lambda_{NT})^{\prime} a^{\prime}a (\Lambda_{NT} \kappa - \lambda_{NT})
  kappa <- c(rho, rho*rho)
  # result <- t(Lambda1 %*% kappa - lambda1) %*% aTa %*% (Lambda1 %*% kappa - lambda1)
  result <- fc_Rrho(Lambda1, lambda1, aTa, kappa)
  return(result)
}
# # plot
# FUNloss <- function(x) {
#   result <- rep(0, length(x))
#   for (i in 1:length(x)) {
#     result[i] <- LossFUN_rho(x[i], Lambda1, lambda1, diag(2))
#   }
#   return(result)
# }
# x <- seq(-1, 1, by=0.01)
# y <- FUNloss(x)
# plot(x,y)

