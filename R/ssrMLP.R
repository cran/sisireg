#
# reconstruction of 2-layer MLP with partial sum optimizer
#

#
# internal: sigmoid-function
#
sigmoidR <- function(x) {
  1 / (1 + exp(-x))
}

#
# internal: error-function for partial sums
#
err_psR <- function(nb, Y, Yp, k=10, fn=4) {
  s <- 0
  for (nbi in nb) {
    tmp <- 0
    for (i in unlist(nbi)) {
      tmp <- tmp + sign(Y[i] - Yp[i])
    }
    s <- s + max(0, abs(tmp)-fn)
  }
  return (s)
}

#
# internal: check partial sums for a given coordinate i
#
check_psR <- function(nb, Y, Yp, i, k=10, fn=4) {
  for (nbi in nb) {
    if (i %in% unlist(nbi)) {
      s <- 0
      for (m in unlist(nbi)) {
        s <- s + sign(Y[m] - Yp[m])
      }
      if (abs(s) > fn)  return (FALSE)
    }
  }
  return (TRUE)
}

#
# internal: factor calculation for partial sums
#
fac_psR <- function(nb, Y, Yp, i, k=10, fn=4) {
  if (check_psR(nb, Y, Yp, i)) yfactor <- 0
  else                         yfactor <- sign(Y[i]-Yp[i]) * (abs(Y[i]-Yp[i])**0.01)
  return (yfactor)
}

#
# internal: error function for squared rediuals (L2)
#
err_l2R <- function(nb, Y, Yp, k=10, fn=4) {
  ret <- sum((Y-Yp)**2)/length(Y)
  return (ret)
}

#
# internal: factor calculation for L2
#
fac_l2R <- function(nb, Y, Yp, i, k=10, fn=4) {
  ret <- Y[i] - Yp[i]
  return (ret)
}

#
# internal: calculate output vector
#
calcOutR <- function(X, W0, W1, W2) {
  O1 <- sigmoidR(W0 %*% t(X))
  O2 <- sigmoidR(W1 %*% O1)
  y <- t(W2 %*% O2)
  return (y)
}  

#
# calcuale prediction for a given model W
#
ssrmlp_predict <- function(X, W) {
  X <- t(apply(X, 1, function(x) (x-W$minX)/(W$maxX-W$minX)))
  X <- cbind(X, rep(1, nrow(X)))
  y <- calcOutR(X, W$W0, W$W1, W$W2)
  y = y * (W$maxY - W$minY)  + W$minY
  return (y)
}

#
# training of the network
#
ssrmlp_train <- function(X, Y, std=TRUE, opt='ps', hl = NULL, W = NULL, k=10, fn=4, eta=0.75, maxIter=1000) {
  
  il <- ncol(X) + 1
  ol <- 1
  n <- length(Y)
  
  if (is.null(hl)) {
    hln <- as.integer(-(il+ol+1)/2 + sqrt((il+ol+1)**2+n))*2
    hl <- c(hln, hln)
  }
  printR('ssrMLP: number of neurons per layer: %i\n', hl[1])
  
  # standardizing to unit interval
  if (is.null(W)) {
    if (std) {
      minX <- apply(X, 2, min)
      maxX <- apply(X, 2, max)
      X <- t(apply(X, 1, function(x) (x-minX)/(maxX-minX)))
      minY <- min(Y)
      maxY <- max(Y)
      Y = (Y - minY) / (maxY - minY)
    } else {
      minX <- 0
      maxX <- 1
      minY <- 0
      maxY <- 1
    }
    W0 <- array(runif(hl[1]*il, -0.5, 0.5), dim=c(hl[1], il))
    W1 <- array(runif(hl[2]*hl[1], -0.5, 0.5), dim=c(hl[2], hl[1]))
    W2 <- array(runif(ol*hl[2], -0.5, 0.5), dim=c(ol, hl[2]))
  } else {
    printR("ssrMLP: re-training...\n")
    W0 <- W$W0
    W1 <- W$W1
    W2 <- W$W2
    
    minX <- W$minX
    maxX <- W$maxX
    minY <- W$minY
    maxY <- W$maxY
    X <- t(apply(X, 1, function(x) (x-minX)/(maxX-minX)))
    Y = (Y - minY) / (maxY - minY)
  }

  Xnb <- X  # coordinates without bias-coordinate
  X <- cbind(X, rep(1, nrow(X)))  
  
  if (opt == 'l2') {
    printR("ssrMLP: optimizing with L2-residuals\n")
    errfct <- err_l2R
    facfct <- fac_l2R
    nb <- NULL
  } else {
    printR("ssrMLP: optimizing with partial sums\n")
    errfct <- err_psR
    facfct <- fac_psR
    d <- as.matrix(dist(Xnb))
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, (il-1)*k+1))
  }
  
  Yp <- calcOutR(X, W0, W1, W2)
  
  meanE <- errfct(nb, Y, Yp)
  printR("ssrMLP: start error: %f\n", meanE)
  minError = meanE 
  minW0 <- W0
  minW1 <- W1
  minW2 <- W2
  
  printR("ssrMLP: %d iterations: [", maxIter)
  chng = FALSE
  
  # mix input
  mixSet <- sample(nrow(X))
  # iterate
  for (c in 1:maxIter) {
    for (i in mixSet) {
      x <- matrix(X[i,], nrow = 1)
      O1 <- as.vector(sigmoidR(W0 %*% t(x)))
      O2 <- as.vector(sigmoidR(W1 %*% O1))
      temp <- W2*O2*(1-O2)
      dW2 <- O2
      dW1 <- t(t(temp)%*%O1) 
      dW0 <- t(O1*(1-O1)*(temp %*% W1))%*% x 
      Yp[i] <- calcOutR(x, W0, W1, W2)  
      yfactor = facfct(nb, Y, Yp, i)
      W0 <- W0 + eta * yfactor * dW0
      W1 <- W1 + eta * yfactor * dW1
      W2 <- W2 + eta * yfactor * dW2
    }
    Yp <- calcOutR(X, W0, W1, W2)
    
    meanE <- errfct(nb, Y, Yp)
    #print(paste0('error (', c, '): ', meanE))
    if (meanE < minError) {
      chng = TRUE
      minError <- meanE
      minW0 <- W0
      minW1 <- W1
      minW2 <- W2
      #print(paste0('minError (', c, '): ', minError))
    }
    if (c %% 50 == 0) {
      if (chng) printR("X")
      else      printR("x")
      chng = FALSE
    }
  }
  printR("]\n")
  printR('ssrMLP: final minError: %f\n', minError)
  return (list(W0=minW0, W1=minW1, W2=minW2, 
               minX=minX, maxX=maxX, minY=minY, maxY=maxY))
}
