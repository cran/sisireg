#
# implements the NARCH model based on the AxE model
#

#
# internal function: get current trend
#
coordR <- function(quotes) {
  n <- length(quotes)
  coord <- c(0,0,0,0,0,0)
  for (period in list(5,10,20)) {
    seq = 1:period
    medreg <- lm(quotes[n-period+1:period] ~ seq, data = as.list(quotes[n-period+1:period]))
    coord[period/5] <- (coefficients(medreg)[2]/quotes[n])
    coord[period/5+1] <- quantile(abs(residuals(medreg)), 0.95)/quotes[n]
  }
  X <- matrix(coord, ncol = 6)
  return (X)
}
  

#
# main function: creates the model from the data
#
axe_narch_model <- function(quotes, T, tgt) {
  # validation
  stopifnot(T == 5 || T == 10 || T ==20)
  stopifnot(length(quotes) > 249)
  stopifnot((tgt == 'trend') || (tgt == 'vola'))
  
  # coordinates and target values
  df1 <- axe(quotes)
  
  # Datensatz bereinigen
  ac <- df1[complete.cases(df1), ]
  # coordinates and target
  x <- data.frame(x1=ac$trend5, x2=ac$trend10, x3=ac$trend20, x4=ac$vola5, x5=ac$vola10, x6=ac$vola20)
  
  if (tgt == 'trend') y <- data.frame(y1=ac[paste0("chng",T)])
  else                y <- data.frame(y1=ac[paste0("risk",T)])
  # als Matrix
  X <- as.matrix(x)
  Y <- as.matrix(y)
  # train
  df_model <- ssrmlp_train(X, Y)
  # Retrain
  df_model <- ssrmlp_train(X, Y, W=df_model, fn=3)
  
  return (df_model)
  
}


#
# main function: calculates prognosis based on quotes and model 
#
axe_narch_predict <- function(model, quotes, tgt) {
  # validation
  stopifnot((tgt == 'trend') || (tgt == 'vola'))
  stopifnot(length(quotes) == 20)
  
  # calculate current coordinates from quotes
  X <- coordR(quotes)
  Yp <- ssrmlp_predict(X, model)
  
  if (tgt == 'trend') prediction <- (Yp + 1) * quotes[length(quotes)]
  else                prediction <- Yp * quotes[length(quotes)] / 2
  
  return (prediction)
}