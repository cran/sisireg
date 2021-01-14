#
# functions for n-dimensional optimizing
#

#
# calculates the maximum partial sum
#  in a given neighborhood
#  based on the observations and the regression function
#
psmaxndR <- function(nb, dat, mu, k) {
  sum <- 0
  for (j in 1:length(dat)) {
    z <- dat[,j]
    m <- unlist(mu[j])
    for (i in nb) {
      if (length(i) != k * 4 + 1) next
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
psplotnd <- function(koord, dat, mu, text = 'Sample') {
    n <- length(koord[,1])
    maxint <- max(n/20,10)
    ps <- array(data = NA, dim = maxint, dimnames = NULL)
    d <- as.matrix(dist(koord))
    # calculating partial sums
    printR("psplotnd: analyze %d intervals: [", maxint)
    for (k in (1:maxint)) {
      nb <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, 4*k+1))
      ps[k] <- psmaxndR(nb, dat, mu, k)
      if (k %% 10 == 0) printR("%d",k/10)
      else              printR("x")
    }
    printR("]\n")
    # plot of the maximal partial sums
    plot(ps, xlim=c(1, maxint), ylim=c(0, max(ps,na.rm=TRUE)+5), type="h",
         main = paste0('3D partial sums (', text, ')'),
         xlab = 'number of neighbors per quadrant', ylab = 'max. partial sum')
    # plot of the quantiles
    fn <- fnR(n, seq(5, 4*maxint+1, by = 4))
    lines(fn, xlim=c(1, maxint), col="red", lwd = 3)
    legend("topleft", inset=c(0.01,0.01),
           legend=c('quantiles'), col=c("red"), lty=1:1)

}


#
# calculates the n-dimensional ssr-function for a single output
#
ssrndSingleR <- function(koord, dat, k = NULL, fn= NULL, iter = 1000) {
  if (is.null(k))   k = maxRunR(length(dat))/2 # neighbors per quandrant
  if (is.null(fn)) fn = fnR(length(dat), k*4)

  printR(paste0("ssrnd: building neighborhood (k=",k,",fn=",fn))
  # n-dimensional neighborhood
  d <- as.matrix(dist(koord))
  #  two observations per axis for minimal surface
  n <- ncol(koord)*4 + 1
  printR(paste0(",nb_sm=",n))
  nb_ind <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, n))
  nb_dst <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i]), n))
  #  2k observations per dimension for partial sum criterion (index only)
  n <- k*4 + 1
  printR(paste0(",nb_ps=",n))
  printR(")\n")
  nb_ps_ind <- lapply(seq_len(ncol(d)), function(i) head(sort(d[,i], index.return=TRUE)$ix, n))
  printR("ssrnd: calculating...\n")
  mu <- .Call(ssrndC, koord, nb_ind, nb_dst, nb_ps_ind, dat, fn, iter)

  # returning all relevant data in a data frame for prediction
  #df <- data.frame(koord=koord, y=dat, mu=mu)
}

#
# calculates the n-dimensional ssr-function
#
ssrnd <- function(koord, dat, k = NULL, fn= NULL, iter = 1000) {
  mu <- lapply(seq_len(ncol(dat)), function(i) ssrndSingleR(koord, dat[,i], k, fn, iter))
  df <- list(koord=koord, dat=dat, mu=mu)
}


#
# calculates the exponential-distance-weighted mean
#
wmean_ndR <- function(mu, dst) {
  if (dst[1] == 0)  return(mu[1])
  g <- exp(-dst)
  wmean <- sum(g*mu)/sum(g)
  return(wmean)
}

#
# prediction for a single output
#
ssrnd_predictSingleR <- function(mu, nb_ind, nb_dst, xx) {
  lapply(seq_len(nrow(xx)), function(i) wmean_ndR(mu[nb_ind[,i]], nb_dst[,i]))
}


#
# calculates the prediction for a given model
#
ssrnd_predict <- function(df_model, xx) {
  dim <- length(df_model$koord)
  dimout <- length(df_model$dat)
  print(paste0("ssrnd_predict: ", dim, " in, ", dimout, " out"))
  n <- (dim)*4

  koord <- df_model$koord
  mu <- df_model$mu

  dst <- t(apply(xx, 1, function(x)  sqrt(rowSums((sweep(koord, 2, x))^2))))

  nb_ind <- apply(dst, 1, function(d) head(sort(d, index.return=TRUE)$ix, n))
  nb_dst <- apply(dst, 1, function(d) head(sort(d), n))

  z <- lapply(seq_len(length(mu)), function(i) ssrnd_predictSingleR(unlist(mu[i]), nb_ind, nb_dst, xx))

}

