#
# rework of ssrMLP: partial sums in a pre-calculated sign list
#


#
# internal: sigmoid-function
#
sigmoidR <- function(x) {
  1 / (1 + exp(-x))
}

#
# internal: continuous version of signum function
#
sign2R <- function(x, a = 100) {
  return (atan(a*x)/(pi*0.5))
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
# internal: builds the k_nb-nearest-neighborhood as matrices (k_nb x n)
#           and the inverse neighborhood index nb_i (list of length n)
#
# nb[,i]     - indices of the k_nb nearest neighbors of point i (incl. i itself)
# nb_dst[,i] - corresponding distances
# nb_i[[j]]  - indices of all points whose neighborhood contains point j
#
buildNeighborhoodR <- function(Xnb, k_nb) {
  n <- nrow(Xnb)
  d <- as.matrix(dist(Xnb))

  nb     <- vapply(seq_len(n), function(i) head(sort.list(d[,i]), k_nb), integer(k_nb))
  nb_dst <- vapply(seq_len(n), function(i) head(as.vector(sort(d[,i])), k_nb), numeric(k_nb))

  # inverse index via split() (O(n log n)) instead of O(n^2) apply()/%in%
  owner <- rep(seq_len(n), each = k_nb)
  nb_i <- split(owner, as.vector(nb))
  nb_i <- nb_i[as.character(seq_len(n))]
  empty <- vapply(nb_i, is.null, logical(1))
  if (any(empty)) nb_i[empty] <- list(integer(0))

  return (list(d=d, nb=nb, nb_dst=nb_dst, nb_i=nb_i))
}

#
# internal: vectorized computation of the partial sum vector s_nb
#           for the full neighborhood matrix nb (k_nb x n)
#
computeSnbR <- function(Y, Yp, nb, k_nb) {
  return (colSums(matrix(sign(Y[nb] - Yp[nb]), nrow = k_nb)))
}

#
# internal: builds the list of weight matrices for an n-layer MLP
#
# dims - vector c(il, hl_1, ..., hl_nl, ol) of layer sizes, i.e. length nl+2
# Wl[[l]] has dimension dims[l+1] x dims[l], l = 1..nl+1
# (Wl[[nl+1]] is the linear output layer)
#
initWeightsR <- function(dims) {
  nw <- length(dims) - 1
  Wl <- vector("list", nw)
  for (l in seq_len(nw)) {
    Wl[[l]] <- array(runif(dims[l+1]*dims[l], -0.5, 0.5), dim=c(dims[l+1], dims[l]))
  }
  return (Wl)
}


#
# factorwise influence indiction (fii) - factor-wise sum of the model weights
#
fii_model <- function(W) {
  il <- ncol(W$Wl[[1]])
  X <- diag(il)
  y <- calcOutNR(X, W$Wl)
  y <- y / sum(y)
  return (c(y))
}

#
# factorwise influence indiction (fii) - factorwise calculation based on input data
#
fii_prediction <- function(W, x) {
  il <- ncol(W$Wl[[1]])
  fii_pred <- matrix(, nrow = il, ncol = 0)
  for (i in 1:nrow(x)) {
    fii <- c((x[i,]-W$minX)/(W$maxX-W$minX), 1) # norming to [0,1]
    fii[is.na(fii)] <- 0 # to avoid NaN from div / 0
    Xd <- diag(fii, il)
    y <- calcOutNR(Xd, W$Wl)
    y <- y * (W$maxY - W$minY) + W$minY
    y <- y / sum(y)
    fii_pred <- cbind(fii_pred, y)
  }
  fii_pred <- apply(fii_pred, 1, function(x) mean(na.omit(x)))
  return (fii_pred)
}


#
# internal: batch forward pass through an n-layer MLP
#
# all hidden layers use the sigmoid activation, the output layer (Wl[[nl+1]])
# is linear (no activation) - as in the original 2-layer architecture
#
calcOutNR <- function(X, Wl) {
  nl <- length(Wl) - 1
  O <- t(X)
  for (l in seq_len(nl)) {
    O <- sigmoidR(Wl[[l]] %*% O)
  }
  return (t(Wl[[nl+1]] %*% O))
}

#
# internal: forward pass for a single sample x (row vector, length dims[1])
#
# returns the list of layer activations A (A[[1]] = x, A[[l+1]] = output of
# hidden layer l) and the (scalar) network output y
#
forwardR <- function(x, Wl) {
  nl <- length(Wl) - 1
  A <- vector("list", nl + 1)
  A[[1]] <- as.vector(x)
  for (l in seq_len(nl)) {
    A[[l+1]] <- as.vector(sigmoidR(Wl[[l]] %*% A[[l]]))
  }
  y <- as.vector(Wl[[nl+1]] %*% A[[nl+1]])
  return (list(A=A, y=y))
}

#
# internal: backpropagation - gradient of the (scalar) output y w.r.t. every
#           weight matrix in Wl, for a single sample (activations A from forwardR)
#
gradR <- function(Wl, A) {
  nl <- length(Wl) - 1
  dW <- vector("list", nl + 1)
  dW[[nl+1]] <- matrix(A[[nl+1]], nrow = 1)   # dy/dW_out (linear output layer)
  delta <- as.vector(Wl[[nl+1]])              # dy/dA[[nl+1]]
  for (l in rev(seq_len(nl))) {
    deltaL <- delta * A[[l+1]] * (1 - A[[l+1]])   # dy/dnet_l
    dW[[l]] <- outer(deltaL, A[[l]])
    if (l > 1) delta <- as.vector(deltaL %*% Wl[[l]])  # dy/dA[[l]]
  }
  return (dW)
}

#
# calculate prediction for a given model W (as returned by ssrmlp_train)
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
  y <- calcOutNR(X, W$Wl)
  y = y * (W$maxY - W$minY)  + W$minY
  return (y[,1])
}

#
# training of the network
#
ssrmlp_train <- function(X, Y, std=TRUE, opt='ps', nl=2, hl=NULL, W=NULL, k=NULL, fn=NULL, eta=0.5, accept = 10, maxIter=1000, alpha=NULL, beta=NULL) {
  
  start.time <- Sys.time()
  
  il <- ncol(X) + 1
  ol <- 1
  n <- length(Y)
  printR('ssrMLP: Data details: n=%i, Input Layer=%i, Output Layer=%i\n', n, il-1, ol)
  
  if (is.null(k))   k = maxRunR(n)
  k_nb <- (il-1)*k+1
  if (is.null(fn)) fn = fnR(n, k_nb)*0.5
  printR('ssrMLP: Partial Sum Parameters: k=%i / nb=%i / fn=%f\n', k, k_nb, fn)
  
  if (nl < 1) stop("ssrMLP: nl (number of hidden layers) must be >= 1")
  if (is.null(hl)) {
    hln <- as.integer(-(il+ol+1)/2 + sqrt((il+ol+1)**2+n))
    hl <- c(rep(hln, nl - 1), hln %/% 2)  # default: all layers with equal neuron count, except last layer
  } else if (length(hl) == 1) {
    hl <- rep(hl, nl)                  # single value given: replicate to all layers
  } else if (length(hl) != nl) {
    stop("ssrMLP: length(hl) must be 1 or equal to nl")
  }
  printR('ssrMLP: number of layers: %i / neurons per layer: %s\n', nl, paste(hl, collapse=" / "))
  
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
    Wl <- initWeightsR(c(il, hl, ol))
  } else {
    printR("ssrMLP: re-training...\n")
    Wl <- W$Wl
    
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
  
  eta_i <- eta

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
            printR("ssrMLP: minimizing squared errors s.t. partial sums\n")
            printR("ssrMLP: parameters: alpha=%f, beta=%f\n", alpha, beta)
            optfct   <- opt_lseR
            optfct_d <- opt_lse_dR
          }, 
          'simple' = {
            printR("ssrMLP: minimizing simplified loss function\n")
            optfct   <- opt_simpleR
            optfct_d <- opt_simple_dR
          },
          {
            printR("ssrMLP: minimizing partial sums\n")
            optfct   <- opt_psR
            optfct_d <- opt_ps_dR
          }
  )
  
  # stating strictness of result
  if (accept > 0) printR("ssrMLP: accepting %i%% deviation in PS criterion\n", accept)
  else            printR("ssrMLP: accepting only strict PS criterion\n")
  
  # neighborhood for partial sums: k times the number of axis plus center
  # (als Matrizen statt Listen - siehe buildNeighborhoodR)
  nbdata <- buildNeighborhoodR(Xnb, k_nb)
  nb     <- nbdata$nb
  nb_dst <- nbdata$nb_dst
  nb_i   <- nbdata$nb_i
  
  # init output vector
  Yp <- calcOutNR(X, Wl)
  # init vector of partial sums (vektorisiert)
  s_nb <- computeSnbR(Y, Yp, nb, k_nb)
  meanE <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)
  printR("ssrMLP: step width: eta=%f\n", eta)
  
  printR("ssrMLP: start error: %f\n", meanE)
  minError = meanE 
  minWl <- Wl
  minps <- psR(s_nb, fn)
  minYp <- Yp
  printR("ssrMLP: start partial sum: %f\n", minps)
  
  printR("ssrMLP: %d iterations: [", maxIter)
  chng = FALSE
  accepted = FALSE
  
  # iterate / mix input in every iteration
  for (c in 1:maxIter) {
    mixSet <- sample(nrow(X))
    Yp_old <- Yp
    for (i in mixSet) {
      x <- X[i,]
      fwd <- forwardR(x, Wl)
      dW <- gradR(Wl, fwd$A)
      # calculate i-th coordinate
      Yp[i] <- fwd$y
      old_sign = sign(Y[i]-Yp_old[i])
      new_sign = sign(Y[i]-Yp[i])
      if (old_sign != new_sign) {
        idx <- nb_i[[i]]
        s_nb[idx] <- s_nb[idx] - old_sign + new_sign
      }
      # calculate gradient  
      yfactor = optfct_d(Y, Yp, nb, nb_dst, s_nb, nb_i, i, k, fn, alpha, beta)
      # do descent
      for (l in seq_along(Wl)) {
        Wl[[l]] <- Wl[[l]] + eta_i * yfactor * dW[[l]]
      }
    }
    Yp <- calcOutNR(X, Wl)
    s_nb <- computeSnbR(Y, Yp, nb, k_nb)
    
    meanE_old <- meanE
    meanE <- optfct(Y, Yp, nb, nb_dst, s_nb, k, fn, alpha, beta)

    if (meanE <= minError) {
      chng = TRUE
      minError <- meanE
      minWl <- Wl
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
      if (psplotnd(Xnb, Y, minYp, resultplot=FALSE, accept=accept)) {
        accepted = TRUE
        break
      }
    }
  }
  printR("]\n")
  printR('ssrMLP: finshed after %i iterations\n', c)
  if (accepted) printR('ssrMLP: final result is acceptable\n')
  else          printR('ssrMLP: final result is not acceptable!\n')
  printR('ssrMLP: final minError: %f\n', minError)
  printR('ssrMLP: final partial sum: %f\n', minps)
  end.time <- Sys.time()
  printR('ssrMLP: training duration %s secs\n', difftime(end.time, start.time, units='secs'))
  return (list(Wl=minWl, 
               minX=minX, maxX=maxX, minY=minY, maxY=maxY,
               X=Xnb, Y=Y, mu=minYp[,1], accepted=accepted))
}

