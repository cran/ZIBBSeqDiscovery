fitZIBB <-
function(dataMatrix, X, ziMatrix, mode = "free", gn = 3, 
                    betastart = matrix(NA,0,0), 
                    psi.start = vector(mode="numeric", length=0), 
                    eta.start = matrix(NA,0,0)) {
  # Initialization
  if (class(dataMatrix) != "matrix") {
    dataMatrix <- as.matrix(dataMatrix)
  }
  nOTU <- nrow(dataMatrix)
  nSample <- ncol(dataMatrix)
  if (is.null(rownames(dataMatrix))) {
    rownames(dataMatrix) <- paste(rep("OTU", nOTU), as.character(c(1:nOTU)), 
                                  sep="")
  }
  if (is.null(colnames(dataMatrix))) {
    colnames(dataMatrix) <- paste(rep("Sample", nSample), 
                                  as.character(c(1:nSample)), sep="")
  }
  
  p <- dim(X)[2]
  q <- dim(ziMatrix)[2]
  n1 <- length(which(X[, 2] == 1))
  n2 <- length(which(X[, 2] == 0))
  #n <- n1 + n2
  Y <- t(dataMatrix)
  Y.c <- as.vector(apply(dataMatrix, 2, sum))
  m <- dim(dataMatrix)[1]
  
  if (mode == "free") {
    # initial value for psi
    phihat <- rep(0, m)
    model.p <- rep(0, m)
    #temp.coef <- (nSample - 1) / nSample
    mean.y = mean(Y.c)
    for (j in (1:m)) {
      #mean.y <- max(Y[, j], 2)
      phat <- Y[, j] / mean.y
      model.p[j] <- mean(phat)
      #phihat[j] <- (1 / (mean.y - 1)) * 
      #  (temp.coef * var(Y[, j]) / (mean.y * model.p[j] * (1 - model.p[j])) - 1)
      phihat[j] <- (1 / (mean.y - 1)) * 
        (var(Y[, j]) / (mean.y * model.p[j] * (1 - model.p[j])) - 1)
    }
    #phihat[phihat <= 0] = 1e-09
    #phihat[phihat >= 1] = 1 - (1e-09)
    #phihat <- rep(0.1, m)
    options(warn = -1)
    psi.out <- log(phihat / (1 - phihat))
    psi.out[is.na(psi.out)] <- mean(psi.out, na.rm = TRUE)
    options(warn = 0)
    # initial value for beta
    Beta.out <- matrix(NA, p, m)
    for (j in (1:m)) {
      phat <- Y[, j] / Y.c
      phat[phat > 1] <- Y[phat > 1, j] / rowSums(Y)[phat > 1]
      phat[is.nan(phat)] <- mean(phat[!is.nan(phat)])
      phat[phat == 0] <- 1 / (2 * Y.c[phat == 0])
      phat[phat == 1] <- 1 - 1 / (2 * Y.c[phat == 1])
      temp.reg <- regression(log(phat / (1 - phat)), X)
      Beta.out[, j] <- temp.reg$betahat
    }
    # initial value for eta
    zeroLM <- lm.fit(ziMatrix, qlogis(colSums(dataMatrix == 0) / m))
    estimate.eta <- as.vector(zeroLM$coefficients)
    eta.out <- matrix(rep(estimate.eta, m), q, m)
    # initial value for others
    Bvar.out <- matrix(NA, p, m)
    pt.out <- rep(NA, m)
  } else if (mode == "constrained") {
    Beta.out <- betastart
    psi.out <- psi.start
    eta.out <- eta.start
    Bvar.out <- matrix(NA, p, m)
    pt.out <- rep(NA, m)
    
    x <- apply(X %*% as.matrix(Beta.out), 2, mean)
    y <- psi.out
    if (gn == 1) {
      Coef <- summary(lm(y ~ x))$coef[, 1]
      psi.out <- Coef[1] + x * Coef[2]
    }
    if (gn == 2) {
      x2 <- x^2
      Coef <- summary(lm(y ~ x + x2))$coef[, 1]
      psi.out <- Coef[1] + x * Coef[2] + x2 * Coef[3]
    }
    if (gn == 3) {
      x2 <- x^2
      x3 <- x^3
      Coef <- summary(lm(y ~ x + x2 + x3))$coef[, 1]
      psi.out <- Coef[1] + x * Coef[2] + x2 * Coef[3] + x3 * Coef[4]
    }
    gamma.out <- Coef
  } else {
    print("Unsupported mode!")
  }
  
  # MLE
  if (mode == "free") {
    fit <- free.estimate(m, p, q, nSample, Beta.out, Bvar.out, psi.out, eta.out, 
                         Y, X, Y.c, ziMatrix)
    Beta.out <- fit$betahat
    Bvar.out <- fit$bvar
    psi.out <- fit$psi
    eta.out <- fit$eta
  } else {
    fit <- constrained.estimate(m, p, q, nSample, Beta.out, Bvar.out, psi.out, 
                                eta.out, gamma.out, Y, X, Y.c, ziMatrix)
    Beta.out <- fit$betahat
    Bvar.out <- fit$bvar
    psi.out <- fit$psi
    eta.out <- fit$eta
  }
  
  
  # calculate p-values
  #k <- nSample - 2 - apply(Y == 0, 2, sum)
  k <- nSample - p - q
  k[k <= 0] <- 1
  if (mode == "free") {
    k <- k - 1
    k[k <= 0] <- 1
    pt.out <- 2 * pt(-abs(as.numeric(Beta.out[2, ]) / 
                            sqrt(as.numeric(Bvar.out[2, ]))), df = k)
  } else {
    pt.out <- 2 * pnorm(-abs(as.numeric(Beta.out[2, ]) / 
                               sqrt(as.numeric(Bvar.out[2, ]))))
  }
  
  #pt.out <- getnewp(as.numeric(Beta.out[2, ]) / sqrt(as.numeric(Bvar.out[2, ])))
  
  # output
  colnames(Beta.out) <- rownames(dataMatrix)
  colnames(Bvar.out) <- rownames(dataMatrix)
  colnames(eta.out) <- rownames(dataMatrix)
  names(pt.out) <- rownames(dataMatrix)
  names(psi.out) <- rownames(dataMatrix)
  if (mode == "free") {
    out <- list(Beta.out, Bvar.out, pt.out, psi.out, eta.out)
    names(out) <- c("betahat", "bvar", "p", "psi", "zeroCoef")
    #out <- list(Beta.out, Bvar.out, pt.out, psi.out, pi.out)
    #names(out) <- c("betahat", "bvar", "p", "psi", "pi")
  } else {
    out <- list(Beta.out, Bvar.out, pt.out, psi.out, eta.out, gamma.out)
    names(out) <- c("betahat", "bvar", "p", "psi", "zeroCoef", "gamma")
    #out <- list(Beta.out, Bvar.out, pt.out, psi.out, pi.out, fit$gamma)
    #names(out) <- c("betahat", "bvar", "p", "psi", "pi", "gamma")
  }
  return(out)
}
