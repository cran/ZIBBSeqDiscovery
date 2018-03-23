free.loglikelihood <-
function(para, X, Y.col, Y.c, ziMatrix) 
{
  #library(base)
  p <- dim(X)[2]
  n <- length(Y.col)
  q <- dim(ziMatrix)[2]
  
  B.col <- para[1:p]
  B.current <- matrix(B.col, p, 1)
  p.current <- exp(X %*% B.current) / (1 + exp(X %*% B.current))
  
  psi <- para[p + 1]
  phi <- exp(psi) / (1 + exp(psi))
  s1 <- as.vector(p.current * (1 - phi) / phi)
  s2 <- as.vector((1 - phi) / phi * (1 - p.current))
  
  eta.col <- para[(p+2):(p+q+1)]
  eta.current <- matrix(eta.col, q, 1)
  patzero.current <- exp(ziMatrix %*% eta.current) / 
    (1 + exp(ziMatrix %*% eta.current))
  
  loglike <- rep(0, n)
  zeroInd <- Y.col == 0
  
  loglike[zeroInd] <- log(patzero.current[zeroInd] + 
                            (1 - patzero.current[zeroInd]) * 
                            beta(s1[zeroInd], Y.c[zeroInd] + s2[zeroInd]) / 
                            beta(s1[zeroInd], s2[zeroInd]))
  
  loglike[!zeroInd] <- log(1 - patzero.current[!zeroInd]) + 
    lchoose(Y.c[!zeroInd], Y.col[!zeroInd]) + 
    lbeta(Y.col[!zeroInd] + s1[!zeroInd], 
          Y.c[!zeroInd] - Y.col[!zeroInd] + s2[!zeroInd]) - 
    lbeta(s1[!zeroInd], s2[!zeroInd])
  
  #loglike <- (1 - z) * 
  #  (lchoose(Y.c, Y.col) + lbeta(Y.col + s1, Y.c - Y.col + s2) - lbeta(s1, s2))
  return(-sum(loglike))
}
