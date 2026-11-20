#
# runtime-optimized rework of ssrnd 
#

#
# internal: print output
#
printR <- function(...) cat(sprintf(...))


#
# internal: calculates the maximum partial sum
#            in a given neighborhood
#           based on the observations and the regression function
#
psmaxndR <- function(nb, dat, mu) {
  sum <- 0
  for (j in 1:ncol(dat)) {
    z <- dat[,j]
    m <- mu[,j]
    for (i in nb) {
      tmp <- abs(sum(sign(z[i] - m[i])))
      if (tmp > sum)  {
        sum <- tmp
      }
    }
  }
  return (sum)
}


#
# plots the partial sum statistic
#
psplotnd <- function(koord, dat, mu, text = 'Sample', maxint = NULL, resultplot = TRUE, accept = 10) {
  n <- length(koord[,1])
  if (is.null(maxint)) maxint <- as.integer(round(max(n/5,10)))
  ps <- array(data = NA, dim = maxint, dimnames = NULL)
  d <- as.matrix(dist(koord))
  # todo: beware of multiple output dimensions
  dat <- as.matrix(dat)
  mu <- as.matrix(mu)
  # calculating partial sums
  if (resultplot) printR("psplotnd: analyze %d intervals: [", maxint)
  for (k in (1:maxint)) {
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), k))
    ps[k] <- psmaxndR(nb, dat, mu)
    if (resultplot) {
      if (k %% 10 == 0) printR("%d",k/10)
      else              printR("x")
    }
  }
  if (resultplot) printR("]\n")
  fn <- fnR(n, 1:maxint)
  failed <- which(fn<ps)
  failed_fn_diff <- round(max((ps-fn)/fn)*100, digits=0)
  failed_fn_nmbr <- round(length(failed)/maxint*100, digits=0)
  if (resultplot) {
    if (all(fn >= ps)) {
      psResult = 'passed'
    } else {
      if ((failed_fn_nmbr <= accept) && (failed_fn_diff <= accept))
        psResult = paste0('accepted (#: ', failed_fn_nmbr, '%, d: ', failed_fn_diff, '%)')
      else 
        psResult = paste0('failed (#: ', failed_fn_nmbr, '%, d: ', failed_fn_diff, '%)')
    }
    # plot of the maximal partial sums
    plot(ps, xlim=c(5, maxint), ylim=c(0, max(ps,na.rm=TRUE)+5), type="h",
         main = paste0(text, ': Partial Sum Test -> ', psResult),
         xlab = 'number of neighbors', ylab = 'max. partial sum')
    # plot of the quantiles
    lines(fn, xlim=c(5, maxint), col="green", lwd = 3)
    # plot violations
    points(failed, ps[failed], col='red', pch=16)
    legend("topleft", inset=c(0.01,0.01),
           legend=c('quantiles', 'violations'), col=c("green", "red"), lty=c(1, NA), pch=c(NA, 16))
  }
  return ((failed_fn_nmbr <= accept) && (failed_fn_diff <= accept))
}


#
# internal: calculates the distance-weighted mean (linear weights)
#
wmean_ndR <- function(mu, dst) {
  if (dst[1] == 0)  return(mu[1]) # input coordinate is model coordinate
  g <- rev(dst)
  wmean <- sum(g*mu)/sum(g)
  return(wmean)
}


#
# calculates the Nearest Neighborhood Estimator
#
ssrndNN <- function(X, Y, k, weighted=FALSE) {
  il <- ncol(X)
  # normalize coordinates and target vector to [0,1]
  minX <- apply(X, 2, min)
  maxX <- apply(X, 2, max)
  X <- apply(X, 1, function(x) (x-minX)/(maxX-minX))
  minY <- min(Y)
  maxY <- max(Y)
  Y <- (Y - minY) / (maxY - minY)
  # redim input
  if (il > 1)  X <- t(X) 
  else              X <- matrix(X)
  
  # calculate distances
  d <- as.matrix(dist(X))
  g <- exp(-d)
  # build neighborhood
  nb <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), k))
  # calculate mean over neighborhood
  if (weighted)  mu <- unlist(lapply(1:length(nb), function(j) sum(Y[unlist(nb[j])]*d[unlist(nb[j])]) / sum(d[unlist(nb[j])]) ))
  else           mu <- unlist(lapply(1:length(nb), function(j) sum(Y[unlist(nb[j])]) / k ))
  
  # return complete mdl-list to be consistent with predict function
  return (list(X=X, Y=Y, mu=mu, lastmu=NULL,
               minX=minX, maxX=maxX, minY=minY, maxY=maxY, n_l1=k,
               k=k, fn=NULL, meanE=NULL, psE=NULL))
}


#
# internal: builds a fixed-size (k x n) nearest-neighbor index/distance
#           matrix from a precomputed full distance matrix d (n x n)
#
# nb[,i]     - indices of the k nearest neighbors of point i (incl. i itself)
# nb_dst[,i] - corresponding distances
#
buildNbMatR <- function(d, k) {
  n <- ncol(d)
  nb     <- vapply(seq_len(n), function(i) head(sort.list(d[,i]), k), integer(k))
  nb_dst <- vapply(seq_len(n), function(i) head(sort(d[,i]), k), numeric(k))
  return (list(nb=nb, nb_dst=nb_dst))
}

#
# internal: builds the inverse neighborhood index nb_i (list of length n)
#           from a fixed-size (k x n) neighbor index matrix nb, i.e.
#           nb_i[[j]] = indices of all points whose neighborhood contains
#           point j - via split() (O(n log n)) instead of O(n^2) apply()/%in%
#
buildInverseNbR <- function(nb, n) {
  k <- nrow(nb)
  owner <- rep(seq_len(n), each = k)
  nb_i <- split(owner, as.vector(nb))
  nb_i <- nb_i[as.character(seq_len(n))]
  empty <- vapply(nb_i, is.null, logical(1))
  if (any(empty)) nb_i[empty] <- list(integer(0))
  return (nb_i)
}

#
# internal: vectorized computation of the partial sum vector s_nb
#           for the full neighborhood matrix nb_sc (n_sc x n)
#
computeSnbNd2R <- function(Y, Yp, nb_sc, n_sc) {
  return (colSums(matrix(sign(Y[nb_sc] - Yp[nb_sc]), nrow = n_sc)))
}

#
# internal: vectorized loss function for L1 curvature (see loss_l1R in
#           ssrnd.R for the per-point reference implementation)
#
# nb (k x n) / nb_dst (k x n): row 1 is the point itself (distance 0),
# rows 2:k are its neighbors used for the curvature estimate
#
loss_l1_matR <- function(Yp, nb, nb_dst) {
  k <- nrow(nb)
  refVal   <- Yp[nb[1,]]
  refMat   <- matrix(refVal, nrow = k-1, ncol = ncol(nb), byrow = TRUE)
  otherVal <- matrix(Yp[nb[2:k, , drop=FALSE]], nrow = k-1)
  distVal  <- nb_dst[2:k, , drop=FALSE]
  return (sum(((refMat - otherVal)/distVal)**2))
}

#
# calculates the prediction for a given model (runtime-optimized: cross-
# distances between test- and model-coordinates via a single matrix
# multiplication instead of apply()/sweep() per test point)
#
ssrnd_predict <- function(xx, mdl) {
  # set model values
  X <- mdl$X
  mu <- mdl$mu
  n <- mdl$n_l1 # number of elements in l1 surface neighborhood
  il <- ncol(X) # coordinates dimension
  # normalizing input coordinates according to model
  xx <- apply(xx, 1, function(x) (x-mdl$minX)/(mdl$maxX-mdl$minX))
  if (il > 1) xx <- t(xx)
  else        xx <- matrix(xx)
  # distances of input coordinates to model coordinates:
  # ||a-b||^2 = ||a||^2 + ||b||^2 - 2*a.b, computed for all pairs at once
  dst2 <- outer(rowSums(xx**2), rowSums(X**2), "+") - 2 * xx %*% t(X)
  dst2[dst2 < 0] <- 0 # guard against numeric noise
  dst <- sqrt(dst2)
  # building neighborhoods
  nb_ind <- apply(dst, 1, function(d) head(sort.list(d), n))
  nb_dst <- apply(dst, 1, function(d) head(sort(d), n))
  # calculating weighted mean
  z <- unlist( lapply( seq_len(nrow(xx)), function(i) wmean_ndR(mu[nb_ind[,i]], nb_dst[,i]) ) )
  # redo normalizing
  z <- z * (mdl$maxY - mdl$minY)  + mdl$minY
  return (z)
}

#
# calculates the n-dimensional ssr-function for a one dimensinal output
# either mdl or (X and Y): re-training can only be done with same (X,Y)
#
# runtime-optimized (algorithmically identical) variant of ssrnd - see
# file header for details
#
ssrnd <- function(X, Y, mdl=NULL, k=NULL, fn=NULL, n_l1=NULL, iter=1000, eps=0.0001) {

  start.time <- Sys.time()
  printR("ssrnd: minimizing l1 curvature s.t. partial sums\n")

  n <- length(Y)
  il <- ncol(X)

  if (is.null(k))        k = maxRunR(length(Y))   # number of neighbors for PS criterion per dimension
  n_sc <- k*il+1                                  # neighborhood for partial sums: maxRun-length per dimension plus 1
  if (is.null(fn))      fn = fnR(n, n_sc)*0.5     # number of allowed unbalanced signs
  if (is.null(n_l1))  n_l1 = ncol(X)*2 + 1        # number of neighbors in l1 surface


  # standardizing to unit interval
  if (is.null(mdl)) {
    minX <- apply(X, 2, min)
    maxX <- apply(X, 2, max)
    X <- apply(X, 1, function(x) (x-minX)/(maxX-minX))
    minY <- min(Y)
    maxY <- max(Y)
    Y <- (Y - minY) / (maxY - minY)
    mu <- NULL
  } else {
    printR("ssrnd: re-training...\n")
    minX <- mdl$minX
    maxX <- mdl$maxX
    minY <- mdl$minY
    maxY <- mdl$maxY
    X <- mdl$X
    Y = mdl$Y
    mu <- mdl$mu
  }

  # unbelievable: special treating of one dimension
  # apply is transposing if X has more than one column
  # and converting to array otherwise
  if (il > 1)  X <- t(X)
  else         X <- matrix(X)

  # n-dimensional neighborhood
  d <- as.matrix(dist(X))

  #  one observations per axis for loss function (plus reference coordinate)
  n_l1_loss <- ncol(X) + 1
  nbLoss <- buildNbMatR(d, n_l1_loss)
  nb_l1_loss <- nbLoss$nb
  nb_l1_dst_loss <- nbLoss$nb_dst

  #  two observations per axis plus middle point for minimal surface
  nbL1 <- buildNbMatR(d, n_l1)
  nb_l1 <- nbL1$nb
  nb_l1_dst <- nbL1$nb_dst
  # neighborhoods including index i
  nb_l1_i <- buildInverseNbR(nb_l1, n)

  #  k observations per dimension for partial sum criterion (index only)
  nbSc <- buildNbMatR(d, n_sc)
  nb_sc <- nbSc$nb
  # neighborhoods including index i
  nb_sc_i <- buildInverseNbR(nb_sc, n)
  printR("ssrnd: building neighborhood: k=%i, n_l1=%s, n_sc=%i, fn=%f\n", k, n_l1, n_sc, fn)

  # init Yp with median value in each neighborhood
  if (is.null(mu)) {
    mu <- colSums(matrix(Y[nb_l1], nrow = n_l1)) / n_l1
  }
  Yp <- Y   # initial values: data itself
  s_nb <- computeSnbNd2R(Y, Yp, nb_sc, n_sc)

  meanE = rep(NA, iter+1)
  psE = rep(NA, iter+1)
  meanE[1] <- loss_l1_matR(Yp, nb_l1_loss, nb_l1_dst_loss)
  psE[1] <- psR(s_nb, fn)
  minError = meanE[1]
  minps <- psE[1]
  minmu <- Yp
  minind <- 0
  printR("ssrnd: start error: %f\n", minError)
  printR("ssrnd: start partial sum: %f\n", minps)

  printR("ssrnd: max %d iterations: [", iter)
  chng = FALSE

  for (c in 1:iter) {
    mixSet <- sample(n)
    Yp_old <- Yp
    s_nb_old <- s_nb
    for (i in mixSet) {
      x <- matrix(X[i,], nrow = 1)
      curr_s_nb <- any(abs(s_nb[nb_sc_i[[i]]]) > fn)
      # calculate new candidate
      nb_l1_0 <- nb_l1[,i][-1]
      nb_l1_dst_0 <- nb_l1_dst[,i][-1]
      Yp[i] <- sum(Yp[nb_l1_0] / nb_l1_dst_0) / sum(1/nb_l1_dst_0)
      # check partial sum
      new_sign = sign(Y[i]-Yp[i])
      old_sign = sign(Y[i]-Yp_old[i])
      # consider constraints, may also be violated before new calculation
      if (any(abs(s_nb[nb_sc_i[[i]]]-old_sign+new_sign) > fn)) {
        if ((old_sign == 0) || (curr_s_nb)) {Yp[i] <- Y[i]; new_sign <- 0}
        else if (old_sign != new_sign)      {Yp[i] <- (Y[i] + Yp_old[i])/2; new_sign <- old_sign}
        else                                {Yp[i] <- Yp_old[i]; new_sign <- old_sign}
      }
      # recalculate i-th partial sums
      if (old_sign != new_sign) {
        idx <- nb_sc_i[[i]]
        s_nb[idx] <- s_nb[idx] - old_sign + new_sign
      }
    }

    # avoid local minima
    if (c %% 20 == 0) {
      #setting back touch points
      eq <- which(Yp == Y) # loop over eq
      for (i in eq) {
        if (Yp[i] == Y[i])  {
          nb_l1_0 <- nb_l1[,i][-1]
          nb_l1_dst_0 <- nb_l1_dst[,i][-1]
          Yp[i] <- sum(Yp[nb_l1_0] / nb_l1_dst_0) / sum(1/nb_l1_dst_0)
        }
      }
      s_nb <- computeSnbNd2R(Y, Yp, nb_sc, n_sc)

    }

    # algorithm is decreasing loss function value in every iteration, PS crit could lead to momentary increase
    meanE[c+1] <- loss_l1_matR(Yp, nb_l1_loss, nb_l1_dst_loss)
    psE[c+1] <- psR(s_nb, fn)
    delta <- minError - meanE[c+1]
    if ((delta > 0) && (psE[c+1] < fn)) {
      chng = TRUE
      minError <- meanE[c+1]
      minps <- psE[c+1]
      minmu <- Yp
      minind <- c
      if (delta <= eps) break # progress too small: stop processing
    }
    if ((c-minind > 100) && (eps > 0.0)) break  # no progress for 100 iterations: stop processing

    if (c %% 50 == 0) {
      if (chng) printR("Y")
      else      printR("X")
      chng = FALSE
    }
  }
  printR("]\n")
  printR('ssrnd: finished after %i iterations\n', c-1)
  printR('ssrnd: final meanE: %f\n', meanE[c+1])
  printR('ssrnd: final partial sum: %f\n', psE[c+1])
  printR('ssrnd: final min Error(%d): %f\n', minind, minError)
  printR('ssrnd: final min partial sum: %f\n', minps)

  end.time <- Sys.time()
  printR('ssrnd: training duration %s secs\n', difftime(end.time, start.time, units='secs'))

  # returning all relevant data in a data frame for prediction
  return (list(X=X, Y=Y, mu=minmu, lastmu=Yp,
               minX=minX, maxX=maxX, minY=minY, maxY=maxY, 
               n_l1=n_l1-1, k=k, fn=fn, meanE=meanE, psE=psE))
}
