free.estimate <-
function(m, p, q, n, betastart, bvarstart, psi.start, eta.start, 
                          Y, X, Y.c, ziMatrix)
{
  Beta.free <- betastart
  Bvar.free <- bvarstart
  psi.free <- psi.start
  eta.free <- eta.start
  #pt.free <- rep(NA, m)
  
  for (j in (1:m)) {
    temp <- try(optim(par = c(betastart[, j], psi.start[j], eta.start[, j]), 
                      fn = free.loglikelihood, X = X, Y.col = Y[, j], 
                      Y.c = Y.c, ziMatrix = ziMatrix, hessian = TRUE,
                      control = list(reltol = 1e-200, maxit = 100)), TRUE)
    #print(class(temp))
    if (class(temp) != "try-error") {
      Beta.free[, j] <- temp$par[1:p]
      #print(Beta.free[, j])
      psi.free[j] <- temp$par[(p + 1)]
      eta.free[, j] <- temp$par[(p+2):(p+q+1)]
      var <- try(diag(solve(temp$hessian)), TRUE)
      if (class(var) != "try-error") {
        #Bvar.free[, j] <- var[1:p] / sum(1 - z[j,], na.rm = TRUE)
        Bvar.free[, j] <- var[1:p]
        if (class(Bvar.free[, j]) == "try-error") {
          Bvar.free[, j] <- NA
        } else if (sum(Bvar.free[, j] <= 0) > 0) {
          Bvar.free[, j] <- NA
        }
      }
    }
  }
  
  out <- list(Beta.free, Bvar.free, psi.free, eta.free)
  names(out) <- c("betahat", "bvar", "psi", "eta")
  return(out)
}
