regression <-
function(y, x) 
{
  n <- length(y)
  p <- dim(x)[2]
  xTx.inv <- solve(t(x) %*% x)
  betahat <- as.vector(xTx.inv %*% t(x) %*% y)
  resid <- y - x %*% betahat
  varhat <- (sum(resid^2)) * diag(xTx.inv) / (n - p)
  se <- sqrt(varhat)
  pval <- 2 * pt(-abs(betahat / se), df = n - p)
  out <- list(betahat, se, betahat / se, n - p, pval)
  names(out) <- c("betahat", "se", "t", "df", "pval")
  return(out)
}
