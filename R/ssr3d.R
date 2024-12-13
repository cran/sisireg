#
# Functions for spatial surfaces
#

#
# help function: pointDistance from library raster not available any more
#
pointDistanceR <- function(a, b, lonlat) {
  b <- as.matrix(b)
  return (sqrt(rowSums((sweep(a, 2, b))^2)))
}

#
# calculates the k nearest neighbors in each quadrant for a reference coordinate x
#
nearneighborsQR <- function(x, koord, k = 1) {
  d <- pointDistanceR(as.matrix(koord), x, lonlat=FALSE)
  # refrence point
  nb0 <- which(d == 0)
  # upper right / north east
  xQ1 <- as.matrix(koord[which(koord$x > x[1] & koord$y >= x[2]),])
  d1 <- pointDistanceR(xQ1, x, lonlat=FALSE)
  nb1 <- which(d %in% head(sort(d1), k))
  # only coordinates lying in the quadrant
  XXQ1 <- koord[nb1,]
  XXQ1 <- XXQ1[which(XXQ1$x > x[1] & XXQ1$y >= x[2]),]
  nb1 <- as.integer(rownames(XXQ1))
  # upper left / north west
  xQ2 <- as.matrix(koord[which(koord$x <= x[1] & koord$y > x[2]),])
  d2 <- pointDistanceR(xQ2, x, lonlat=FALSE)
  nb2 <- which(d %in% head(sort(d2), k))
  # only coordinates lying in the quadrant
  XXQ2 <- koord[nb2,]
  XXQ2 <- XXQ2[which(XXQ2$x <= x[1] & XXQ2$y > x[2]),]
  nb2 <- as.integer(rownames(XXQ2))
  # lower left / south west
  xQ3 <- as.matrix(koord[which(koord$x < x[1] & koord$y <= x[2]),])
  d3 <- pointDistanceR(xQ3, x, lonlat=FALSE)
  nb3 <- which(d %in% head(sort(d3), k))
  # only coordinates lying in the quadrant
  XXQ3 <- koord[nb3,]
  XXQ3 <- XXQ3[which(XXQ3$x < x[1] & XXQ3$y <= x[2]),]
  nb3 <- as.integer(rownames(XXQ3))
  # lower right / south east
  xQ4 <- as.matrix(koord[which(koord$x >= x[1] & koord$y < x[2]),])
  d4 <- pointDistanceR(xQ4, x, lonlat=FALSE)
  nb4 <- which(d %in% head(sort(d4), k))
  # only coordinates lying in the quadrant
  XXQ4 <- koord[nb4,]
  XXQ4 <- XXQ4[which(XXQ4$x >= x[1] & XXQ4$y < x[2]),]
  nb4 <- as.integer(rownames(XXQ4))
  # aggregate and disambiguate
  if (k == 1) nb <- c(nb1, nb2, nb3, nb4)
  else        nb <- c(nb0, nb1, nb2, nb3, nb4)
  return (nb[!duplicated(nb)])
}

#
# calculates the k-quadrant-neighborhood for all coordinates in the given grid
#
nearneighboursGridQR <- function(koord, k) {
  nb <- apply(koord, 1, nearneighborsQR, koord, k)
}


#
# calculates the k-neighborhood for a given reference point
#
nearneighborsR <- function(x, koord, k = 1) {
  d <- pointDistanceR(as.matrix(koord), x, lonlat=FALSE)
  nb <- which(d %in% head(sort(d), k*4+1))
  return (nb)

}

#
# calculates the k-neighborhood for all coordinates in the given grid
#
nearneighboursGridR <- function(koord, k) {
  nb <- apply(koord, 1, nearneighborsR, koord, k)
  if (typeof(nb) != "list") {
    nb <- lapply(seq_len(ncol(nb)), function(i) nb[,i])
  }
  return (nb)
}

#
# calculates the maximum partial sum in a given neighborhood based on the observations and the regression function
#
psmax3dR <- function(nb, z, mu, k) {
  sum <- 0
  for (i in nb) {
    if (length(i) != k * 4 + 1) next
    tmp <- abs(sum(sign(z[i] - mu[i])))
    if (tmp > sum)  {
      sum <- tmp
    }

  }
  return (sum)
}

#
# print output
#
printR <- function(...) cat(sprintf(...))

#
# calculates the partial sums and plots the statistic
#
psplot3d <- function(koord, z, mu, text = 'Sample') {
  n <- length(koord[,1])
  maxint <- max(n/20,10)
  ps <- array(data = NA, dim = maxint, dimnames = NULL)
  # calculating partial sums
  printR("psplot3d: analyze %d intervals: [", maxint)
  for (k in (1:maxint)) {
     nb <- nearneighboursGridR(koord, k) # corresponds to 4k + 1 neighbors
     ps[k] <- psmax3dR(nb, z, mu, k)
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
# calculates the 3-dimensional ssr-function
#
ssr3d <- function(koord, dat, k= NULL, fn= NULL, iter = 1000) {
  if (is.null(k)) k = maxRunR(length(dat))/2 # neighbors per quandrant
  if (is.null(fn)) fn = fnR(length(dat), k*4)
  printR("ssr3d: Calculating neighborhood...\n")
  nb <- nearneighboursGridR(koord, k)  # for the partial sum
  nb1 <- nearneighboursGridQR(koord, 1) # for the surface
  printR("ssr3d: Calculating the minimal surface...\n")
  mu <- .Call(ssr3dC, koord, nb, nb1, dat, fn, iter)
  # returning all relevant data in a data frame for prediction
  df <- data.frame(koord=koord, z=dat, mu=mu)
}

#
# calculates the reciprocal-distance-weighted mean
#
wmeanR <- function(koord, mu, x) {

  d <- pointDistanceR(koord, x, lonlat=FALSE)
  if (0 %in% d) {
    return (mu[which(d == 0)])
  }
  wmean <- sum(mu/d) / sum(1/d)
}

#
# calculates the exponential-distance-weighted mean
#
wmean_expR <- function(koord, mu, x) {

  d <- pointDistanceR(koord, x, lonlat=FALSE)
  g <- exp(-d)
  wmean <- sum(g*mu)/sum(g)
}

#
# calculates the weighted mean based on the 4-point minimal surface
#
wmean_msR <- function(koord, mu, x) {

  x1 <- koord[1,]
  x2 <- koord[2,]
  x3 <- koord[3,]
  x4 <- koord[4,]

  # x-Achse
  dx2 <- pointDistanceR(x2, x, lonlat = FALSE)
  dx3 <- pointDistanceR(x3, x, lonlat = FALSE)
  d23 <- pointDistanceR(x2,x3, lonlat = FALSE)
  s23 <- (dx2 + dx3 + d23)/2
  if (s23*(s23-dx3)*(s23-dx2)*(s23-d23) < 0) {
    h23 <- 0
  } else {
    h23 <- 2/d23*sqrt(abs(s23*(s23-dx3)*(s23-dx2)*(s23-d23)))
  }

  dx1 <- pointDistanceR(x1, x, lonlat = FALSE)
  dx4 <- pointDistanceR(x4, x, lonlat = FALSE)
  d14 <- pointDistanceR(x1,x4, lonlat = FALSE)
  s14 <- (dx1 + dx4 + d14)/2
  if (s14*(s14-dx4)*(s14-dx1)*(s14-d14) < 0) {
    h14 <- 0
  } else {
    h14 <- 2/d14*sqrt(abs(s14*(s14-dx4)*(s14-dx1)*(s14-d14)))
  }

  # y-Achse
  d34 <- pointDistanceR(x3,x4, lonlat = FALSE)
  s34 <- (dx3 + dx4 + d34)/2
  if (s34*(s34-dx3)*(s34-dx4)*(s34-d34) < 0) {
    h34 <- 0
  } else {
    h34 <- 2/d34*sqrt(abs(s34*(s34-dx3)*(s34-dx4)*(s34-d34)))
  }

  d12 <- pointDistanceR(x1,x2, lonlat = FALSE)
  s12 <- (dx1 + dx2 + d12)/2
  if (s12*(s12-dx1)*(s12-dx2)*(s12-d12) < 0.001) {
    h12 <- 0
  } else {
    h12 <- 2/d12*sqrt(abs(s12*(s12-dx1)*(s12-dx2)*(s12-d12)))
  }

  # x-intercept based on x3
  rx <- h23/(h23+h14)
  # y-intercept based on x3
  ry <- h34/(h34+h12)
  if (h34+h12 - max(d14, d23) > 0.1) {
    if (abs(h34-h12) < 0.1) {
      ry <- 0.5
    } else {
      ry <- h34/(h34-h12)
    }
  }

  # line-segment west-east
  a1 <- mu[3] + rx * (mu[4] - mu[3])
  b1 <- mu[2] + rx * (mu[1] - mu[2])
  m1 <- a1+ry*(b1-a1)

  return (m1)

}


#
# calculates the prediction based on exponential weights
#
predict3dSingleR <- function(x, koord, mu) {
  nb1 <- unlist(nearneighborsQR(x, koord, k = 1))
  koord__ <- data.frame(x=koord$x[nb1], y=koord$y[nb1])
  prd <- wmean_expR(koord__, mu[nb1], x)
}

#
# calculates the prediction based on minimal surfaces
#
predict3dSingle_msR <- function(x, koord, mu) {
  nb1 <- unlist(nearneighborsQR(x, koord, k = 1))
  koord__ <- data.frame(x=koord$x[nb1], y=koord$y[nb1])
  prd <- wmean_msR(koord__, mu[nb1], x)
}


#
# calculates the prediction for a given model
#
ssr3d_predict <- function(df_model, xy, ms = FALSE) {
  koord <- data.frame(x = df_model$koord.x, y = df_model$koord.y)
  if (ms) {
    z <- apply(xy, 1, predict3dSingle_msR, koord, df_model$mu)
  } else {
    z <- apply(xy, 1, predict3dSingleR, koord, df_model$mu)
  }
}

