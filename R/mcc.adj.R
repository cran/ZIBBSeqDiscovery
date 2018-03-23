mcc.adj <-
function(out.fitZIBB, dataMatrix, X, ziMatrix, K=4) {
  out <- out.fitZIBB
  z <- ziMatrix[,2]
  cutpoint <- as.vector(quantile(z, seq(0,1,b=1/K)))
  n <- length(z)
  m <- dim(dataMatrix)[1]
  newz <- rep(0,n)
  for (k in (1:K)){
    if (k<K) which.group=which((z >= cutpoint[k] & z<cutpoint[k+1]) == TRUE)
    if (k==K) which.group=which((z >= cutpoint[k] & z<=cutpoint[k+1]) == TRUE)
    newz[which.group]=k
  }
  out.cov.con <- mcc::getAmoment(as.matrix(dataMatrix), X[,2], newz)
  out.mcc <- mcc::getbetap.A(out.cov.con, A=NULL, fix.obs=FALSE)
  out.mcc$twosidedp[is.nan(out.mcc$twosidedp)] <- NA
  idx.na <- is.na(out$p)
  out$p[idx.na] <- out.mcc$twosidedp[idx.na]
  ## warning
  if (sum(rowSums(dataMatrix == 0) == (n-1)) > 0) {
    warning("Detect OTUs for which only one nonzero count across all the samples. It is expected that the p values are NAs for those OTUs.")
  }
  return(out)
}
