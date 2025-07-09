#
# functions for n-dimensional optimizing (rewoked)
#

#
# internal: calculates the maximum partial sum
#            in a given neighborhood
#           based on the observations and the regression function
#
psmaxnd2R <- function(nb, dat, mu) {
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
psplotnd2 <- function(koord, dat, mu, text = 'Sample', maxint = NULL, resultplot = TRUE, accept = 10) {
  n <- length(koord[,1])
  if (is.null(maxint)) maxint <- as.integer(round(max(n/5,10)))
  ps <- array(data = NA, dim = maxint, dimnames = NULL)
  d <- as.matrix(dist(koord))
  # todo: beware of multiple output dimensions
  dat <- as.matrix(dat)
  mu <- as.matrix(mu)
  # calculating partial sums
  if (resultplot) printR("psplotnd2: analyze %d intervals: [", maxint)
  for (k in (1:maxint)) {
    nb <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), k))
    ps[k] <- psmaxnd2R(nb, dat, mu)
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
wmean_nd2R <- function(mu, dst) {
  if (dst[1] == 0)  return(mu[1]) # input coordinate is model coordinate
  g <- rev(dst)
  wmean <- sum(g*mu)/sum(g)
  return(wmean)
}


#
# calculates the Nearest Neighborhood Estimator
#
ssrnd2NN <- function(X, Y, k, weighted=FALSE) {
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
# calculates the prediction for a given model
#
ssrnd2_predict <- function(mdl, xx) {
  # set model values
  X <- mdl$X
  mu <- mdl$mu
  n <- mdl$n_l1 # number of elements in l1 surface neighborhood
  il <- ncol(X) # coordinates dimension
  # normalizing input coordinates according to model
  xx <- apply(xx, 1, function(x) (x-mdl$minX)/(mdl$maxX-mdl$minX))
  if (il > 1) xx <- t(xx) 
  else        xx <- matrix(xx)
  # distances of input coordinates to model coordinates
  dst <- t(apply(xx, 1, function(x)  sqrt(rowSums((sweep(X, 2, x))^2))))
  # building neighborhoods
  nb_ind <- apply(dst, 1, function(d) head(sort.list(d), n))
  nb_dst <- apply(dst, 1, function(d) head(sort(d), n))
  # calculating weighted mean
  z <- unlist( lapply( seq_len(nrow(xx)), function(i) wmean_nd2R(mu[nb_ind[,i]], nb_dst[,i]) ) )
  # redo normalizing
  z <- z * (mdl$maxY - mdl$minY)  + mdl$minY
  return (z)
}


#
# internal: loss function for L1 curvature
#
loss_l1R <- function(Yp, nb, nb_dst) {
  l1 <- sum(unlist(
    lapply(1:length(nb), function(j) sum(unlist(
      lapply(2:length(unlist(nb[j])), function(i) ((Yp[unlist(nb[j])[1]]-Yp[unlist(nb[j])[i]])/unlist(nb_dst[j])[i])**2)
    )) )
  ) )
  return (l1)
}


#
# calculates the n-dimensional ssr-function for a one dimensinal output 
# either mdl or (X and Y): re-training can only be done with same (X,Y)
#
ssrnd2 <- function(X, Y, mdl=NULL, k=NULL, fn=NULL, n_l1=NULL, iter=1000, eps=0.0001) {
  
  start.time <- Sys.time()
  printR("ssrnd2: minimizing l1 curvature s.t. partial sums\n")
  
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
    printR("ssrnd2: re-training...\n")
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
  nb_l1_loss <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), n_l1_loss))
  nb_l1_dst_loss <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), n_l1_loss))
  
  #  two observations per axis plus middle point for minimal surface
  nb_l1 <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), n_l1))
  nb_l1_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), n_l1))
  # neighborhoods including index i
  nb_l1_i <- lapply(seq_len(ncol(d)), function(i) which(TRUE==(apply(t(matrix(unlist(nb_l1),n_l1)), 1, function(x) i %in% x))))
  
  #  k observations per dimension for partial sum criterion (index only)
  nb_sc <- lapply(seq_len(ncol(d)), function(i) head(sort.list(d[,i]), n_sc))
  # neighborhoods including index i
  nb_sc_i <- lapply(seq_len(ncol(d)), function(i) which(TRUE==(apply(t(matrix(unlist(nb_sc),n_sc)), 1, function(x) i %in% x))))
  printR("ssrnd2: building neighborhood: k=%i, n_l1=%s, n_sc=%i, fn=%f\n", k, n_l1, n_sc, fn)
  
  # init Yp with median value in each neighborhood 
  if (is.null(mu)) {
    mu <- unlist(lapply(1:length(nb_l1), function(j) sum(Y[unlist(nb_l1[j])]) / n_l1 ))
  } 
  Yp <- Y   # initial values: data itself
  s_nb <- unlist(lapply(seq_len(n), function(i) sum(sign(Y[unlist(nb_sc[i])] - Yp[unlist(nb_sc[i])]))))
  
  meanE = rep(NA, iter+1)
  psE = rep(NA, iter+1)
  meanE[1] <- loss_l1R(Yp, nb_l1_loss, nb_l1_dst_loss)
  psE[1] <- psR(s_nb, fn)
  minError = meanE[1] 
  minps <- psE[1]
  minmu <- Yp
  minind <- 0
  printR("ssrnd2: start error: %f\n", minError)
  printR("ssrnd2: start partial sum: %f\n", minps)
  
  printR("ssrnd2: max %d iterations: [", iter)
  chng = FALSE
  
  for (c in 1:iter) {
    mixSet <- sample(n)
    Yp_old <- Yp
    s_nb_old <- s_nb
    for (i in mixSet) {
      x <- matrix(X[i,], nrow = 1)
      curr_s_nb <- any(abs(s_nb[unlist(nb_sc_i[i])]) > fn)
      # calculate new candidate
      nb_l1_0 <- unlist(nb_l1[i])[-1]
      nb_l1_dst_0 <- unlist(nb_l1_dst[i])[-1]
      Yp[i] <- sum(Yp[nb_l1_0] / nb_l1_dst_0) / sum(1/nb_l1_dst_0)
      # check partial sum
      new_sign = sign(Y[i]-Yp[i])
      old_sign = sign(Y[i]-Yp_old[i])
      # consider constraints, may also be violated before new calculation
      if (any(abs(s_nb[unlist(nb_sc_i[i])]-old_sign+new_sign) > fn)) {
        if ((old_sign == 0) || (curr_s_nb)) {Yp[i] <- Y[i]; new_sign <- 0}
        else if (old_sign != new_sign)      {Yp[i] <- (Y[i] + Yp_old[i])/2; new_sign <- old_sign}
        else                                {Yp[i] <- Yp_old[i]; new_sign <- old_sign}
      }
      # recalculate i-th partial sums
      if (old_sign != new_sign) {
        s_nb[unlist(nb_sc_i[i])] <- lapply(unlist(nb_sc_i[i]), function(j) s_nb[j] - old_sign + new_sign)
        s_nb <- unlist(s_nb)
      }
    }
    
    # avoid local minima
    if (c %% 20 == 0) {
      #setting back touch points
      eq <- which(Yp == Y) # loop over eq
      for (i in eq) {
        if (Yp[i] == Y[i])  {
          nb_l1_0 <- unlist(nb_l1[i])[-1]
          nb_l1_dst_0 <- unlist(nb_l1_dst[i])[-1]
          Yp[i] <- sum(Yp[nb_l1_0] / nb_l1_dst_0) / sum(1/nb_l1_dst_0)
        }
      }
      s_nb <- unlist(lapply(seq_len(n), function(i) sum(sign(Y[unlist(nb_sc[i])] - Yp[unlist(nb_sc[i])]))))
      
    }
    
    # algorithm is decreasing loss function value in every iteration, PS crit could lead to momentary increase
    meanE[c+1] <- loss_l1R(Yp, nb_l1_loss, nb_l1_dst_loss)
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
  printR('ssrnd2: finished after %i iterations\n', c-1)
  printR('ssrnd2: final meanE: %f\n', meanE[c+1])
  printR('ssrnd2: final partial sum: %f\n', psE[c+1])
  printR('ssrnd2: final min Error(%d): %f\n', minind, minError)
  printR('ssrnd2: final min partial sum: %f\n', minps)
  
  end.time <- Sys.time()
  printR('ssrnd2: training duration %s secs\n', difftime(end.time, start.time, units='secs'))
  
  # returning all relevant data in a data frame for prediction
  return (list(X=X, Y=Y, mu=minmu, lastmu=Yp,
               minX=minX, maxX=maxX, minY=minY, maxY=maxY, n_l1=n_l1-1,
               k=k, fn=fn, meanE=meanE, psE=psE))
}

