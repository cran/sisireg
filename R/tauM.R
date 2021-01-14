#
# calculates the trendbased correlation for a given maturity
#

#
# calculates the trend function
#
trendR <- function(x) {

  # mit ssr-L1
  df <- data.frame(x=seq(1,length(x)), y=x)
  xl1 <- ssr(df, ne = FALSE)
  trend <- diff(xl1)

  return (trend)
}


#
# calculates the trendbased correlation (Metzner's tau)
#
tauM <- function(x, y) {

  trend_x <- trendR(x)
  trend_y <- trendR(y)
  tauM <- sum(sign(trend_x)*sign(trend_y), na.rm=TRUE)/(length(trend_x))

}



