#
# rework of ssrMLP: partial sums in a pre-calculated sign list
#

#
# internal: matrix-safe multi-dimension transpose
#
tR <- function(v) {
  return (aperm(v, perm=rev(1:length(dim(v)))))
}

#
# internal: add an additional dimension
#
addDimR <- function(v, v1, v2) {
  if (is.null(dim(v1))) d1 <- length(v1)
  else                  d1 <- dim(v1)
  if (is.null(dim(v2))) d2 <- length(v2)
  else                  d2 <- dim(v2)
  return (array(v, dim=c(d1, d2)))
}


#
# internal: continuous version of signum function
#
sign2R <- function(x, a = 100) {
  return (atan(a*x)/(pi*0.5))
}

#
# internal: first derivative of continuous version of signum function
#
sign2_dR <- function(x, a = 100) {
  return ( 2*a/(pi*((a*x)^2)+1) )
}

#
# internal: partial sum functional
#
psR <- function(s_nb, fn) {
  return (sum(pmax(0, abs(s_nb) - fn)))
}

#
# definitions of error functions and their first derivatives
#

#
# internal: loss function for least squared residuals
#
loss_lseR <- function(y, yp) {
  return (sum((y-yp)**2))
}
#
# internal: first derivative of lse loss function
#
loss_lse_dR <- function(y, yp, nb, nb_dst, nb_i, i) {
  return ((y[i]-yp[i]))
}

#
# internal: loss function for L1 curvature / currently not used
#
loss_l1R <- function(Yp, nb, nb_dst) {
  l1 <- sum(unlist(
    lapply(1:length(nb), function(j) sum(unlist(
      #      lapply(2:length(unlist(nb[j])), function(i) ((Yp[unlist(nb[j])[1]]-Yp[unlist(nb[j])[i]])/unlist(nb_dst[j])[i])**2)
      #     only one difference in each corrdinate
      lapply(2:2, function(i) ((Yp[unlist(nb[j])[1]]-Yp[unlist(nb[j])[i]])/unlist(nb_dst[j])[i])**2)
    )) )
  ) )
  return (l1)
}
#
# internal: first derivative of l1 curvature loss function / currently not used
#
loss_l1_dR <- function(y, Yp, nb, nb_dst, nb_i, i) {
  # similar to l1 but consider every summand which contains index i
  l1_d <- sum(unlist(
    lapply(unlist(nb_i[i]), function(j) sum(unlist(
      #      lapply(2:length(unlist(nb[j])), function(k) (Yp[unlist(nb[j])[1]]-Yp[unlist(nb[j])[k]])/(unlist(nb_dst[j])[k]**2))
      #     only one difference in each corrdinate
      lapply(2:2, function(k) (Yp[unlist(nb[j])[1]]-Yp[unlist(nb[j])[k]])/(unlist(nb_dst[j])[k]**2))
    )) )
  ) )
  return (l1_d)
}

#
# internal: penalty function for partial sums / currently not used
#
penalty_psR <- function(y, yp, s_nb, k, fn) {
  return (sum(pmax(0, (s_nb - fn))**2 + pmin(0, (fn + s_nb))**2)/(length(y)*k))
}
#
# internal: first derivative for partial sum penalty function / currently not used
#
penalty_ps_dR <- function(y, yp, s_nb, nb_i, i, k, fn) {
  h <- s_nb[unlist(nb_i[i])]
  h_d <- sign2_dR(h)
  return (sum(2*h_d*(pmax(0, (h-fn)) + pmin(0, (h+fn))))/(length(y)*k))
}

#
# internal: function to be minimized: l1 curvature s.t. partial sums criterion / currently not used
#
opt_l1R <- function(y, yp, nb, nb_dst, s_nb, k, fn, alpha, beta) {
  return (alpha * loss_l1R(y, nb, nb_dst) + beta * penalty_psR(y, yp, s_nb, k, fn))
}
#
# internal: first derivative of optimizing function: l1 s.t. ps / currently not used
#
opt_l1_dR <- function(y, yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta) {
  return (alpha * loss_l1_dR(y, yp, nb, nb_dst, nb_i, i) + beta * penalty_ps_dR(y, yp, s_nb, nb_i, i, k, fn))
}

#
# internal: function to be minimized: squared error s.t. partial sums criterion
#
opt_lseR <- function(y, yp, nb, nb_dst, s_nb, k, fn, alpha, beta) {
  return (alpha * loss_lseR(y, yp) + beta * opt_psR(y, yp, nb, nb_dst, s_nb, k, fn, alpha, beta))
}
#
# internal: first derivative of optimizing function: lse s.t. ps
#
opt_lse_dR <- function(y, yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta) {
  return (alpha * loss_lse_dR(y, yp, nb, nb_dst, nb_i, i) + beta * opt_ps_dR(y, yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta))
}

#
# internal: function to be minimized: numerically stable and working alternative
#
opt_simpleR <- function(y, yp, nb, nb_dst, s_nb, k, fn, alpha, beta) {
  return (sum(abs(y - yp)))
}
#
# internal: first derivative of optimizing function: alternative
#
opt_simple_dR <- function (y, yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta) {
  return (sign2R(y[i] - yp[i]))
}

#
# internal: error-function for partial sums
#
opt_psR <- function(y, yp, nb, nb_dst, s_nb, k, fn, alpha, beta) {
  return (psR(s_nb, fn))
}
#
# internal: first derivative ps functional, approximating sign function with power function
#
opt_ps_dR <- function(y, yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta) {
  h <- s_nb[unlist(nb_i[i])]
  if (sum(pmax(0, h-fn) - pmin(0, h+fn)) == 0) yfactor <- 0.0
  else                                         yfactor <- sign(y[i]-yp[i]) * (abs(y[i]-yp[i])**0.01)
  return (yfactor)
}

#
# internal: check Armijo condition
#
checkArmijoR <- function(meanE, mean_eta, zeta, eta_i, yfactor, dW0, dW1, dW2) {
  return (meanE - mean_eta >= -zeta*eta_i*(
    t(yfactor*c(as.vector(dW0),as.vector(dW1),as.vector(dW2))) %*% 
      (yfactor*c(as.vector(dW0),as.vector(dW1),as.vector(dW2)))))
}


#
# internal: calculate Armijo step width
#
armijoR <- function(W0, W1, W2, optfct, yfactor, dW0, dW1, dW2, x, i, Y, Yp, Yp_old, meanE, nb, nb_i, nb_dst, s_nb, k, fn, alpha, beta) {
  zeta = 0.25
  eta_m <- 0.9
  eta_i <- 1.0 / eta_m
  # if no descending direction found: leave and return small standard value
  mean_eta <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)
  if (!checkArmijoR(meanE, mean_eta, zeta, 0.01, yfactor, dW0, dW1, dW2)) return (eta_m)
  #if (meanE < mean_eta) return (0.1) also possible but with poorer results
  
  # Armijo inequality: f(x) - f(x+eta^q*p) >= -zeta*eta^q df(x)p^T 
  # definitions: f = optfct, x = W{1,2,3}, p = yfactor*dW{1,2,3}, df(x) = yfactor*dW{1,2,3}
  repeat { 
    eta_i <- eta_i * eta_m
    # calculate W anew with current eta
    W0_eta <- W0 + eta_i * yfactor * dW0
    W1_eta <- W1 + eta_i * yfactor * dW1
    W2_eta <- W2 + eta_i * yfactor * dW2
    # calculate Yp anew with current W 
    Yp[i] <- calcOutR(x, W0_eta, W1_eta, W2_eta)  
    old_sign = sign(Y[i]-Yp_old[i])
    new_sign = sign(Y[i]-Yp[i])
    if (old_sign != new_sign) {
      s_nb[unlist(nb_i[i])] <- lapply(unlist(nb_i[i]), function(j) s_nb[j] - old_sign + new_sign)
      s_nb <- unlist(s_nb)
    }
    # calculate mean anew with Yp 
    mean_eta <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)
    # assuming that we start with inequality not true and we stop if inequality is true -> largest eta_i is found
    if (checkArmijoR(meanE, mean_eta, zeta, eta_i, yfactor, dW0, dW1, dW2) || (eta_i < 0.01)) {
      break
    } 
  }  
  return (eta_i)
}

#
# training of the network
#
ssrmlp2_train <- function(X, Y, std=TRUE, opt='ps', hl=NULL, W=NULL, k=NULL, fn=NULL, eta=0.5, accept = 10, maxIter=1000, alpha=NULL, beta=NULL) {
  
  start.time <- Sys.time()
  
  il <- ncol(X) + 1
  ol <- 1
  n <- length(Y)
  printR('ssrMLP2: Data details: n=%i, Input Layer=%i, Output Layer=%i\n', n, il-1, ol)
  
  if (is.null(k))   k = maxRunR(n)
  k_nb <- (il-1)*k+1
  if (is.null(fn)) fn = fnR(n, k_nb)*0.5
  printR('ssrMLP2: Partial Sum Parameters: k=%i / nb=%i / fn=%f\n', k, k_nb, fn)
  
  if (is.null(hl)) {
    hln <- as.integer(-(il+ol+1)/2 + sqrt((il+ol+1)**2+n))
    hl <- c(hln*2, hln)
  }
  printR('ssrMLP2: number of neurons per layer: %i / %i\n', hl[1], hl[2])
  
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
  if (il > 2) X <- t(X) 
  else        X <- matrix(X)
  
  Xnb <- X  # coordinates without bias-coordinate
  X <- cbind(X, rep(1, nrow(X)))  # bind additional bias neuron
  
  if (eta == 0) eta_i <- 0.1     # Armijo step width
  else          eta_i <- eta

  # set Lagrange parameter if not already set
  if (is.null(alpha)) {
    alpha <- 0.5
  }
  if (is.null(beta)) {
    beta <- 1-alpha
  }
  
  # options: lse + ps, simple optimizing, ps (standard)
  switch (opt,
          'lse' = {
            printR("ssrMLP2: minimizing squared errors s.t. partial sums\n")
            printR("ssrMLP2: parameters: alpha=%f, beta=%f\n", alpha, beta)
            optfct   <- opt_lseR
            optfct_d <- opt_lse_dR
          }, 
          'simple' = {
            printR("ssrMLP2: minimizing simplified loss function\n")
            optfct   <- opt_simpleR
            optfct_d <- opt_simple_dR
          },
          {
            printR("ssrMLP2: minimizing partial sums\n")
            optfct   <- opt_psR
            optfct_d <- opt_ps_dR
          }
  )
  
  # stating strictness of result
  if (accept > 0) printR("ssrMLP2: accepting %i%% deviation in PS criterion\n", accept)
  else            printR("ssrMLP2: accepting only strict PS criterion\n")
  
  # neighborhood for partial sums: k times the number of axis plus center
  d <- as.matrix(dist(Xnb))
  nb <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), k_nb))
  nb_dst <- lapply(seq_len(ncol(d)), function(i) head(as.vector(sort(d[,i])), k_nb)) # will be used for error function derivative
  # neighborhoods containing index i
  nb_i <- lapply(seq_len(ncol(d)), function(i) which(TRUE==(apply(t(matrix(unlist(nb),k_nb)), 1, function(x) i %in% x))))
  
  # init output vector
  Yp <- calcOutR(X, W0, W1, W2)
  # init vector of partial sums
  s_nb <- unlist(lapply(seq_len(ncol(d)), function(i) sum(sign(Y[unlist(nb[i])] - Yp[unlist(nb[i])]))))
  meanE <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)
  if (eta == 0) printR("ssrMLP2: step width: armijo\n")
  else          printR("ssrMLP2: step width: eta=%f\n", eta)
  
  printR("ssrMLP2: start error: %f\n", meanE)
  minError = meanE 
  minW0 <- W0
  minW1 <- W1
  minW2 <- W2
  minps <- psR(s_nb, fn)
  minYp <- Yp
  printR("ssrMLP2: start partial sum: %f\n", minps)
  
  printR("ssrMLP2: %d iterations: [", maxIter)
  chng = FALSE
  accepted = FALSE
  
  # iterate / mix input in every iteration
  for (c in 1:maxIter) {
    mixSet <- sample(nrow(X))
    Yp_old <- Yp
    s_nb_old <- s_nb
    for (i in mixSet) {
      x <- matrix(X[i,], nrow = 1)
      O1 <- as.vector(sigmoidR(W0 %*% tR(x)))
      O2 <- as.vector(sigmoidR(W1 %*% O1))
      temp <- W2*O2*(1-O2)
      temp <- addDimR(temp, 1, temp) 
      dW2 <- O2
      dW1 <- tR(temp) %*% tR(addDimR(O1, O1, 1))
      dW0 <- tR(O1*(1-O1)*(temp %*% W1)) %*% tR(addDimR(x, x, 1)) 
      # calculate i-th coordinate
      Yp[i] <- calcOutR(x, W0, W1, W2)  
      old_sign = sign(Y[i]-Yp_old[i])
      new_sign = sign(Y[i]-Yp[i])
      if (old_sign != new_sign) {
        s_nb[unlist(nb_i[i])] <- lapply(unlist(nb_i[i]), function(j) s_nb[j] - old_sign + new_sign)
        s_nb <- unlist(s_nb)
      }
      # calculate gradient  
      yfactor = optfct_d(Y, Yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta)
      # calculate step width
      if (eta == 0) {
        eta_i <- armijoR(W0, W1, W2, optfct, yfactor, dW0, dW1, dW2, x, i, Y, Yp, Yp_old, meanE, nb, nb_i, nb_dst, s_nb, k, fn, alpha, beta)
      }
      # do descent
      W0 <- W0 + eta_i * yfactor * dW0
      W1 <- W1 + eta_i * yfactor * dW1
      W2 <- W2 + eta_i * yfactor * dW2
    }
    Yp <- calcOutR(X, W0, W1, W2)
    s_nb <- unlist(lapply(seq_len(ncol(d)), function(i) sum(sign(Y[unlist(nb[i])] - Yp[unlist(nb[i])]))))
    
    meanE_old <- meanE
    meanE <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)

    if (meanE <= minError) {
      chng = TRUE
      minError <- meanE
      minW0 <- W0
      minW1 <- W1
      minW2 <- W2
      minps <- psR(s_nb, fn)
      minYp <- Yp
    }
    if (c %% 50 == 0) {
      if (chng) {
        printR("Y")
      } else {
        printR("x")
      } 
      chng = FALSE
      # check, if ps criterion already satisfied
      if (psplotnd2(Xnb, Y, minYp, resultplot=FALSE, accept=accept)) {
        accepted = TRUE
        break
      }
    }
  }
  printR("]\n")
  printR('ssrMLP2: finshed after %i iterations\n', c)
  if (accepted) printR('ssrMLP2: final result is acceptable\n')
  else          printR('ssrMLP2: final result is not acceptable!\n')
  printR('ssrMLP2: final minError: %f\n', minError)
  printR('ssrMLP2: final partial sum: %f\n', minps)
  end.time <- Sys.time()
  printR('ssrMLP2: training duration %s secs\n', difftime(end.time, start.time, units='secs'))
  return (list(W0=minW0, W1=minW1, W2=minW2, 
               minX=minX, maxX=maxX, minY=minY, maxY=maxY,
               X=Xnb, Y=Y, Yp=minYp, accepted=accepted))
}

