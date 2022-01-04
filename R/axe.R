#
# calculation of trends and volatility of financial time series
#


#
# internal function: calculating trend and volatility
#
trend_volaR <- function(period, df) {
  trend <- rep(NA, nrow(df))
  vola <- rep(NA, nrow(df))
  seq = 1:period

  printR(paste0("calculating ", period, "-day trend and vola (", nrow(df), "): ["))
  for (i in period:nrow(df)) {
    if (i %% 100 == 0)      printR("%d",i)
    else if (i %% 10 == 0)  printR("x")
    medreg <- lm(df$quotes[i-period+1:period] ~ seq, data = as.list(df$quotes[i-period+1:period]))
    trend[i] <- coefficients(medreg)[2]/df$quotes[i]
    vola[i]  <- quantile(abs(residuals(medreg)), 0.95)/df$quotes[i]
  }
  df[paste0("trend",period)] <- trend
  df[paste0("vola",period)] <- vola
  printR("]\n")
  return(df)
}

#
# internal function: calculates the quotes change for 'period' days
#
changesR <- function(period, df) {
  printR(paste0("calculating ", period, "-day changes...\n"))
  tmp <- c(tail(df$quotes, nrow(df)-period), rep(NA, period))
  df[paste0("chng",period)] <- (tmp - df$quotes) / df$quotes
  return(df)
}


#
# internal function: calculates the level of risk
#
risklevelR <- function(period, df) {
  printR(paste0("calculating ", period, "-day risk level...\n"))
  df[paste0("risk",period)] <- c(unlist(lapply(seq_len(nrow(df)-(period-1)), function(i) (max(df$quotes[i:(i+period-1)])-min(df$quotes[i:(i+period-1)])))), rep(NA, (period-1))) / df$quotes
  return(df)
}


#
# main function. creates the data frame with model data
#
axe <- function(quotes) {
  df1 <- data.frame(quotes = quotes)
  for (frist in list(5,10,20)) {
    df1 <- trend_volaR(frist, df1)
    df1 <- changesR(frist, df1)
    df1 <- risklevelR(frist, df1)
  }
  return (df1)
}
