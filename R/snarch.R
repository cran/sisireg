#
# Implementation of the S-NARCH model
#

#
# calculates the trend for period
#
sntrendR <- function(dat, period) {
  mu <- rep(NA, length(dat))
  seq = 1:period

  # start values
  medreg <- lm(dat[1:period] ~ seq, data = as.list(dat[1:period]))
  mu[seq] <- coefficients(medreg)[1] + (seq-1)*coefficients(medreg)[2]
  trend <- coefficients(medreg)[2]

  seq0 = 0:period
  for (i in (period+1):length(dat)) {
    mu[i] <- mu[i-1] + trend
    vz <- sum(sign(mu[i-seq0] - dat[i-seq0]))
    # check run
    if (abs(vz) > period) {
      s <- sign(mu[i]-dat[i])
      # back to beginning of the run
      #while (sign(mu[i]-dat[i]) == s) i <- i - 1
      j <- i - period
      # new trend ...
      m <- min(j+period-1, length(dat))
      seq <- 0:(m-j)
      medreg <- lm(dat[j:m] ~ seq, data = as.list(dat[j:m]))
      mu[j+seq] <- coefficients(medreg)[1] + seq*coefficients(medreg)[2]
      trend <- coefficients(medreg)[2]
      # ... and continue
      mu[i] <- mu[i-1] + trend
    }
  }
  return (mu)
}

#
# calculate the volatility for period
#
snvolaR <- function(dat, mu, period) {
  r <- abs(dat - mu)
  # trend of vola
  vl <- sntrendR(r, period)
  # quantile of noise distribution
  q <- quantile(abs(r-vl), 0.99)
  # vola function
  vl <- vl * q
  return (vl)
}


#
# creates the data frame with trends and vola
#
snarch <- function(dat) {

  tr20 <- sntrendR(dat, 20)
  vl20 <- snvolaR(dat, tr20, 20)
  tr10 <- sntrendR(dat, 10)
  vl10 <- snvolaR(dat, tr10, 10)
  tr5 <- sntrendR(dat, 5)
  vl5 <- snvolaR(dat, tr5, 5)
  df <- data.frame(tr20, vl20, tr10, vl10, tr5, vl5)
  return (df)
}

