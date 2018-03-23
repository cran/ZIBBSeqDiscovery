constrained.estimate <-
function(m, p, q, n, betastart, bvarstart, psi.start, 
                                 eta.start, gamma.start, Y, X, Y.c, ziMatrix, 
                                 gn = 3)
{
  Beta.model <- betastart
  Bvar.model <- bvarstart
  psi.model <- psi.start
  eta.model <- eta.start
  Coef <- gamma.start
  
  for (j in (1:m)) {
    temp <- try(optim(par = c(Beta.model[, j], eta.model[, j]), 
                      fn = constrained.loglikelihood, X = X, Y.col = Y[, j], 
					  Y.c = Y.c, coeff = Coef, ziMatrix = ziMatrix, 
					  hessian = TRUE), TRUE)
    if (class(temp) != "try-error") {
      Beta.model[, j] <- temp$par[1:p]
      eta.model[, j] <- temp$par[(p+1):(p+q)]
      var <- try(diag(solve(temp$hessian)), TRUE)
      if (class(var) != "try-error") {
        #Bvar.free[, j] <- var[1:p] / sum(1 - z[j,], na.rm = TRUE)
        Bvar.model[, j] <- var[1:p]
        if (class(Bvar.model[, j]) == "try-error") {
          Bvar.model[, j] <- NA
        } else if (sum(Bvar.model[, j] <= 0) > 0) {
          Bvar.model[, j] <- NA
        }
      }
    }
  }
  
  x <- apply(X %*% as.matrix(Beta.model), 2, mean)
  if (gn == 1) {
    psi.model <- Coef[1] + x * Coef[2]
  }
  if (gn == 2) {
    x2 <- x^2
    psi.model <- Coef[1] + x * Coef[2] + x2 * Coef[3]
  }
  if (gn == 3) {
    x2 <- x^2
    x3 <- x^3
    psi.model <- Coef[1] + x * Coef[2] + x2 * Coef[3] + x3 * Coef[4]
  }
  
  #out <- list(Beta.model, Bvar.model, psi.model, Coef)
  #names(out) <- c("betahat", "bvar", "psi", "gamma")
  out <- list(Beta.model, Bvar.model, psi.model, eta.model)
  names(out) <- c("betahat", "bvar", "psi", "eta")
  return(out)
}
