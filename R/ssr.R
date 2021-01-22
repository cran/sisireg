#
# Functions for one-dimensional regression
#

#
# calculates the 95%-quantile for the maximum run length
#
maxRunR <- function(n) {
  k <- as.integer(3.3+1.44*log(n))
}

#
# calculates the 95%-quantile of partial sums
#
fnR <- function(n, k) {
  fn <- sqrt(1+2.33*log(n)*k)
  fn <- apply(matrix(c(fn, k), length(k), 2), 1, min)
}

#
# calculates the maximum partial sum for a given interval length
#
psmaxR <- function(dat, mu, k) {
  n <- length(dat)
  maxint <- n/5 + 5
  s <- sign(dat-mu)
  psmax <- max(abs(coredata(rollapply(zoo(s[1:n]), k, sum, align = 'center'))))
}


#
# counts the number of local extremal values in a discrete regression function
#
numberOfExtremaR <- function(mu) {

  n <- length(mu)

  anzMin <- 0
  anzMax <- 0
  slope <- sign(mu[2] - mu[1])

  for (i in (3:n)) {
    if (sign(mu[i] - mu[i-1]) != slope) {
      if (slope < 0)  anzMin <- anzMin + 1
      else            anzMax <- anzMax + 1
      slope <- sign(mu[i] - mu[i-1])
    }
  }
  return(anzMin + anzMax)

}


#
# calculates the partial sums and plots the statistic
#
psplot <- function(dat, mu, text = 'Sample') {
  n <- length(dat)
  maxint <- n/5 + 5
  ps <- array(data = NA, dim = n, dimnames = NULL)
  # claculating the partial sums
  for (k in (1:maxint)) {
    ps[k] <- psmaxR(dat, mu, k)
  }
  # plot of the maximum partial sums
  plot(ps, xlim=c(5, maxint), type="h", main = paste0('partial sums (', text, ')'),
       xlab = 'interval length', ylab = 'max. partial sum')
  # plot of the quantiles
  fn <- fnR(n, seq(1, maxint))
  lines(fn, xlim=c(5, maxint), col="red", lwd = 3)
  legend("topleft", inset=c(0.01,0.01),
         legend=c('quantiles'), col=c("red"), lty=1:1)

}


#
# check, if a given function satisfies the partial sum criterion
#
psvalid <- function(dat, mu) {
  n <- length(dat)
  maxint <- as.integer(n/5)
  # loop over the interval lengths
  for (k in (5:maxint)) {
    if (psmaxR(dat, mu, k) > fnR(n,k)) {
      print(paste0("invalid partial sum detected for interval length k=", k))
      return(FALSE)
    }
  }
  return(TRUE)
}


#
# check, if a given function satisfies the maximum run criterion
#
runvalid <- function(dat, mu, k = NULL) {
  n <- length(dat)
  if (is.null(k)) k <- maxRunR(n)

  s <- sign(dat-mu)
  runmax <- max(abs(coredata(rollapply(zoo(s[1:n]), k+1, sum, align = 'center'))))

  if (runmax > k) {
    print("invalid maximum run detected")
    return(FALSE)
  }
  return(TRUE)
}


#
# minimum statistic for equidistant regression
#
ssr_minR <- function(dat, funk, y1 = NULL, yn = NULL, ps = TRUE, simanz = 10000) {
  # start parameters
  n <- length(dat)
  if (ps) {
    k <- maxRunR(n)
    fn <- fnR(n,k)*0.66
  } else {
    fn <- maxRunR(n)
  }
  mu <- ssrR(dat, funk, y1, yn, fn, ps, simanz)
  valid <- if (ps) psvalid(dat, mu) else runvalid(dat, mu)
  anzExOpt <- numberOfExtremaR(mu)
  anzEx <- anzExOpt
  print(paste0("start: Ex=", anzEx, ", valid=", valid))
  # loop
  while (valid && (anzEx <= anzExOpt) && (fn > 0)) {
    fn <- fn - 1
    mu <- ssrR(dat, funk, y1, yn, fn, ps, simanz)
    valid <- if (ps) psvalid(dat, mu) else runvalid(dat, mu)
    anzEx <- numberOfExtremaR(mu)
    print(paste0("iteration fn=", fn, ": Ex=", anzEx, ", valid=", valid))
  }
  # calculate final function
  fn <- fn + 1
  mu <- ssrR(dat, funk, y1, yn, fn, ps, simanz)
}


#
# wrapper for the iterative solver QSOR (equidistant)
#
ssrR <- function(dat, funk, y1 = NULL, yn = NULL, fn = 0, ps = TRUE, simanz = 10000) {

  # parameters
  n <- length(dat)
  if (ps) {
    if (fn < 2) fn <- max(2, log(n) - 2);
    k <- as.integer(3.3+1.44*log(n)) # 0.95-quantile of max. run length
    k2 <- as.integer(k/2)
    k1 <- k - 2*k2
  } else {
    if (fn < 2) fn <- as.integer(3.3+1.44*log(n)) # 0.95-quantile of max. run length
    k <- fn
    k2 <- as.integer(k/2)
    k1 <- k - 2*k2
  }

  # start values s: 1:k2: median, k2+1:n-k2-1: rollapply(median), n-k2:n: median
  s <- array(dim=n)
  s[1:k2] <- rep(median(dat[1:k2]), k2)
  s[(k2+k1):(n-k2)] <- coredata(rollapply(zoo(dat[1:n]), k, median, align = 'center'))
  s[(n-k2+1):n] <- rep(median(dat[(n-k2+1):n]), k2)
  # boundary values
  if (!is.null(y1)) s[1] <- y1
  if (!is.null(yn)) s[n] <- yn

  result <- .C(ssrC,
               as.integer(funk),
               as.double(dat),
               mu = as.double(s),
               as.integer(n),
               as.integer(fn),
               as.integer(if (ps) 1 else 0),
               as.integer(simanz))
  result$mu
}


#
# minimum statistic for partial sum criterion
#
ssr_ne_minR <- function(df, simanz = 10000) {
  # sorting: quicker when performed here than in each separate iteration
  colnames(df) <- c("x", "y")
  df <- df[order(df$x),]
  # start parameters
  n <- nrow(df)
  k <- maxRunR(n)
  fn <- fnR(n,k)*0.66
  mu <- ssr_neR(df, fn, simanz)$y
  valid <- psvalid(df$y, mu)
  anzExOpt <- numberOfExtremaR(mu)
  anzEx <- anzExOpt
  print(paste0("start: Ex=", anzEx, ", valid=", valid))
  # loop
  while (valid && (anzEx <= anzExOpt) && (fn > 0)) {
    fn <- fn - 1
    mu <- ssr_neR(df, fn, simanz)$y
    valid <- psvalid(df$y, mu)
    anzEx <- numberOfExtremaR(mu)
    print(paste0("iteration fn=", fn, ": Ex=", anzEx, ", valid=", valid))
  }
  # calculate final regression function
  fn <- fn + 1
  dfl1 <- data.frame("x"=df$x, "y"=ssr_neR(df, fn, simanz)$y)
}

#
# wrapper for non-equidistant QSOR-solver, L1 only
#
ssr_neR <- function(df, fn = 0, simanz = 10000) {
  colnames(df) <- c("x", "y")
  # sorting
  df <- df[order(df$x),]
  # parameters
  n <- nrow(df)
  if (fn < 2) fn = max(2, log(n) - 2);
  k <- maxRunR(n) # 0.95-quantile of max. run length
  k2 <- as.integer(k/2)
  k1 <- k - 2*k2

  # start values s: 1:k2: median, k2+1:n-k2-1: rollapply(median), n-k2:n: median
  s <- array(dim=n)
  s[1:k2] <- rep(median(df$y[1:k2]), k2)
  s[(k2+k1):(n-k2)] <- coredata(rollapply(zoo(df$y[1:n]), k, median, align = 'center'))
  s[(n-k2+1):n] <- rep(median(df$y[(n-k2+1):n]), k2)
  # boundary values
  y1 <- s[1]
  yn <- s[n]

  result <- .C(ssr_neC,
     as.double(df$x),
     as.double(df$y),
     mu = as.double(s),
     as.integer(n),
     as.integer(fn),
     as.integer(simanz))
  data.frame("x"=df$x, "y"=result$mu)
}



#
# calculates the ssr-model
#
ssr <- function(df, y1 = NULL, yn = NULL, fn = 0, iter = 10000, minStat = FALSE, ne = TRUE, l1 = TRUE, ps = TRUE) {

  funk <- if (l1) 1 else 2
  colnames(df) <- c("x", "y")

  # non-equidistant
  if (ne) {
    # minimum statistic
    if (minStat) {
      ssr_ne_minR(df, iter)
    } else {
      ssr_neR(df, fn, iter)
    }
  } else {
    # minimum statistic
    if (minStat) {
      ssr_minR(df$y, funk, y1, yn, ps, iter)
    } else {
      ssrR(df$y, funk, y1, yn, fn, ps, iter)
    }
  }

}



#
# one-dimensional QSOR-prediction
# x = coordinates (possibly unsorted),
# mu = model regression function (sorted),
# x_ = coordinate for prediction
#
ssr_ne_pred_singleR <- function(x_, df) {
  colnames(df) <- c("x", "y")
  # sorting
  df <- df[order(df$x),]

  n <- nrow(df)
  x <- df$x
  mu <- df$y

  # ATTENTION: error handling: x[i] = x[i-1]
  if (x_ <= min(x)) {
    y_ = mu[1] - (x[1] - x_)*(mu[2] - mu[1])/(x[2] - x[1])
  } else if (x_ >= max(x)) {
    y_ = mu[n] + (x_ - x[n])*(mu[n] - mu[n-1])/(x[n] - x[n-1])
  } else {
    # searching for index
    j <- 1
    while (x_ < x[j]) j <- j + 1
    i <- j
    while (x[i] < x_ ) i <- i + 1
    # calculating point on a line-segment
    y_ = mu[i] + (x_-x[i])*(mu[j] - mu[i])/(x[j]- x[i]);
  }

  return(y_)
}


#
# prediction based on a given model
#
ssr_predict <- function(df, xx) {
  colnames(df) <- c("x", "y")
  # sorting
  df <- df[order(df$x),]

  yy <- array(as.numeric(unlist(lapply(xx, ssr_ne_pred_singleR, df))), dim = length(xx))
}
