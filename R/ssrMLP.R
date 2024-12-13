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
# internal: 1st derivative from sigmoid-function
#
d1sigmoidR <- function(x) {
  sigmoidR(x)*(1-sigmoidR(x))
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
# definition of error functions and their first derivative
#

#
# internal: error-function for partial sums
#
err_psR <- function(nb, nb_dst=NULL, dim=NULL, Y, Yp, k=10, fn=4) {
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
# internal: factor calculation for partial sums (derivative of error function regarding Yp[i])
#
fac_psR <- function(nb, nb_dst=NULL, dim=NULL, Y, Yp, i, k=10, fn=4) {
  if (check_psR(nb, Y, Yp, i)) yfactor <- 0
  else                         yfactor <- sign(Y[i]-Yp[i]) * (abs(Y[i]-Yp[i])**0.01)
  return (yfactor)
}

#
# internal: error function for squared residuals (lse)
#
err_lseR <- function(nb=NULL, nb_dst=NULL, dim=NULL, Y, Yp, k=10, fn=4) {
  ret <- sum((Y-Yp)**2)/length(Y)
  return (ret)
}
#
# internal: factor calculation for lse (derivative of error function regarding Yp[i])
#
fac_lseR <- function(nb=NULL, nb_dst=NULL, dim=NULL, Y, Yp, i, k=10, fn=4) {
  ret <- Y[i] - Yp[i]
  return (ret)
}

#
# internal: error-function for partial sums combined with least squares
#
err_ps_lseR <- function(nb, dist, dim, Y, Yp, k=10, fn=4, alpha=0.0001) {
  # partial sum
  s <- 0
  for (nbi in nb) {
    tmp <- 0
    for (i in unlist(nbi)) {
      tmp <- tmp + sign(Y[i] - Yp[i])
    }
    s <- s + max(0, abs(tmp)-fn)
  }
  # additional L1
  s <- s + alpha*sum((Yp-Y)**2)
  return (s)
}
#
# internal: factor function for combined PS-LSE error 
#
fac_ps_lseR <- function(nb, dist, dim, Y, Yp, i, k=10, fn=4, alpha=0.0001) {
  yfactor <- alpha*(Y[i]-Yp[i])
  if (!check_psR(nb, Y, Yp, i)) yfactor <- yfactor + sign(Y[i]-Yp[i]) * (abs(Y[i]-Yp[i])**0.01)
  return (yfactor)
}

#
# internal: error-function for partial sums combined with L1 curvature of regression function
#
err_ps_l1R <- function(nb, dist, dim, Y, Yp, k=10, fn=4, alpha=0.0001) {
  # partial sum
  s <- 0
  for (nbi in nb) {
    tmp <- 0
    for (i in unlist(nbi)) {
      tmp <- tmp + sign(Y[i] - Yp[i])
    }
    s <- s + max(0, abs(tmp)-fn)
  }
  # additional L1: take only one observation per axis plus center
  l1 <- 0
  # iterating through 2 lists simultaneously: for-loop / nb and dist have same dimension
  for (j in 1:length(dist)) {
    nbi <- unlist(nb[j])
    dst <- unlist(dist[j])
    i0 <- nbi[1] # neighborhoods are sorted by distance ascending, center is first index in neighborhood
    # take only first dim values: dimension of Input plus one
    for (i in 1:length(nbi)) {
      if (dst[i] != 0) {
        l1 <- l1 + ((Yp[i0]-Yp[nbi[i]])/dst[i])**2
      }
    }
  }
  # final value  
  s <- s + alpha * l1
  
  return (s)
}
#
# internal: factor function (first derivative) of combined PS-L1 error function
#
fac_ps_l1R <- function(nb, dist, dim, Y, Yp, ind, k=10, fn=4, alpha=0.0001) {
  yfactor <- 0.0
  if (!check_psR(nb, Y, Yp, ind)) {
    yfactor <- sign(Y[ind]-Yp[ind]) * (abs(Y[ind]-Yp[ind])**0.01)
  }
  # use only direct neighbors: one per axis plus center
  l1 <- 0
  for (i in 1:length(nb)) {
    nbi <- unlist(nb[i])
    dst <- unlist(dist[i])
    if (ind %in% head(nbi, dim)) {
      i0 <- nbi[1] # neighborhoods are sorted by distance ascending, center is first index in neighborhood
      for (j in 2:length(nbi)) {
        if (dst[j] != 0) {
          l1 <- l1 + ((Yp[i0]-Yp[nbi[j]])/dst[j])
        }
      }
    }
  }
  yfactor <- yfactor + alpha * l1
  #print(paste(ind, ":", l1, "/", yfactor))
  return (yfactor)
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
# factor-wise sum of the model weights
#
fii_model <- function(W) {
  X <- diag(ncol(W$W0))
  O1 <- sigmoidR(W$W0 %*% t(X))
  O2 <- sigmoidR(W$W1 %*% O1)
  y <- t(W$W2 %*% O2)
  y = y / (sum(y))
  return (c(y))
} 

#
# factorwise calculation
#
fii_prediction <- function(W, x) {
  fii_model <- matrix(, nrow = ncol(W$W0), ncol = 0)
  for (i in 1:nrow(x)) {
    fii <- c((x[i,]-W$minX)/(W$maxX-W$minX),1) # norming to [0,1]
    fii[is.na(fii)] <- 0 # to avoid NaN from div / 0
    X <- diag(fii, ncol(W$W0))
    O1 <- sigmoidR(W$W0 %*% t(X))
    O2 <- sigmoidR(W$W1 %*% O1)
    y <- t(W$W2 %*% O2)
    y = y * (W$maxY - W$minY) + W$minY
    y = y / (sum(y))
    fii_model <- cbind(fii_model, y)
  }
  fii_model <- apply(fii_model, 1, function(x) mean(na.omit(x)))
  return(fii_model)
}


#
# calculate prediction for a given model W
#
ssrmlp_predict <- function(X, W) {
  il <- ncol(X)
  X <- apply(X, 1, function(x) (x-W$minX)/(W$maxX-W$minX))
  X[is.na(X)] <- 0 # to avoid NaN from div / 0
  # unbelievable: special treating of one dimension necessary:
  #    'apply' is transposing if X has more than one column and converting to array otherwise
  if (il > 1) X <- t(X) else
    X <- matrix(X)

  X <- cbind(X, rep(1, nrow(X)))
  y <- calcOutR(X, W$W0, W$W1, W$W2)
  y = y * (W$maxY - W$minY)  + W$minY
  return (y)
}

#
# training of the network
#
ssrmlp_train <- function(X, Y, std=TRUE, opt='ps', hl=NULL, W=NULL, k=10, fn=4, eta=0.75, maxIter=1000, facfct_ex=NULL, errfct_ex=NULL, alpha=NULL) {
  
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
      X <- apply(X, 1, function(x) (x-minX)/(maxX-minX))
      X[is.na(X)] <- 0 # to avoid NaN from div / 0
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
    X <- apply(X, 1, function(x) (x-minX)/(maxX-minX))
    Y = (Y - minY) / (maxY - minY)
  }
  
  # unbelievable: special treating of one dimension
  # apply is transposing if X has more than one column
  # and converting to array otherwise
  if (il > 2) X <- t(X) else
    X <- matrix(X)
  
  Xnb <- X  # coordinates without bias-coordinate
  X <- cbind(X, rep(1, nrow(X)))  # bind additional bias neuron
  
  switch(opt,
  'lse'= {
    printR("ssrMLP: optimizing with least-squares-residuals\n")
    errfct <- err_lseR
    facfct <- fac_lseR
    nb <- NULL
    nb_cur_ind <- NULL
    nb_cur_dst <- NULL
  },
  'ps_lse'= {
    printR("ssrMLP: optimizing with combined partial sum and least-squares-residuals criterion\n")
    errfct <- err_ps_lseR
    facfct <- fac_ps_lseR
    if (is.null(alpha)) {
      alpha = 0.1
    }
    # neighborhood for partial sums: k times the number of axis plus center
    d <- as.matrix(dist(Xnb))
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, (il-1)*k+1))
    nb_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), (il-1)*k+1)) # will be used for error function derivative
  },
  'ps_l1'= {
    printR("ssrMLP: optimizing with combined partial sum and l1 curvature criterion\n")
    errfct <- err_ps_l1R
    facfct <- fac_ps_l1R
    if (is.null(alpha)) {
      alpha = 0.0001
    }
    # neighborhood for partial sums: k times the number of axis plus center
    d <- as.matrix(dist(Xnb))
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, (il-1)*k+1))
    nb_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), (il-1)*k+1)) # will be used for error function derivative
  },
  'ext'= {
    printR("ssrMLP: optimizing with external error function\n")
    errfct = errfct_ex # use external error function given by parameter
    facfct = facfct_ex # use external factor function given by parameter
    # neighborhood for partial sums: k times the number of axis plus center
    d <- as.matrix(dist(Xnb))
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, (il-1)*k+1))
    nb_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), (il-1)*k+1)) # will be used for error function derivative
  },
  { # standard: 'ps'
    printR("ssrMLP: optimizing with partial sums\n")
    errfct <- err_psR
    facfct <- fac_psR
    # neighborhood for partial sums: k times the number of axis plus center
    d <- as.matrix(dist(Xnb))
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, (il-1)*k+1))
    nb_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), (il-1)*k+1)) # will be used for error function derivative
  }
  )
  
  # set if not already set
  if (is.null(alpha)) {
    alpha = 0.0
  }
  
  Yp <- calcOutR(X, W0, W1, W2)
  
  meanE <- errfct(nb, nb_dst, il, Y, Yp)
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
      yfactor = facfct(nb, nb_dst, il, Y, Yp, i)
      W0 <- W0 + eta * yfactor * dW0
      W1 <- W1 + eta * yfactor * dW1
      W2 <- W2 + eta * yfactor * dW2
    }
    Yp <- calcOutR(X, W0, W1, W2)
    
    meanE <- errfct(nb, nb_dst, il, Y, Yp)
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
