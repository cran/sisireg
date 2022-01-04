% File sisireg/man/axe_narch_model.Rd
\name{axe_narch_model}
\alias{axe_narch_model}
\title{implementation of the AxE model based on the ssr-MLP}
\description{
Trains a 2-layer MLP with a given time series of quotes with price changes 
or volatility as target value. The coordinates (or independent factors) are given through
the AxE model)
}
\usage{
axe_narch_model(quotes, T, tgt)
}
\arguments{
  \item{quotes}{array with observations.}
  \item{T}{period: T = 5, 10 or 20.}
  \item{tgt}{target variable: tgt = 'trend' or 'vola'.}
}
\value{
  \item{model}{the trained model for prediction.}
}
\examples{
\donttest{
set.seed(1234)
n <- 250
s <- 13000 + cumsum(rnorm(n))
T = 20
# create model for 5-day trend
model <- axe_narch_model(s, T, 'trend')
# calculate prognosis for trend 
s_ <- s[n] + cumsum(rnorm(20))
s_T <- axe_narch_predict(model, s_, 'trend')
# plot the results
plot(seq(1:20), s_, type = "l", 
    xlim = c(0,21+T), ylim = c(min(s_, s_T)-5, max(s_, s_T)+5))
points(20+T, s_T, col='red', pch = 16)
# create model for 5-day vola
model <- axe_narch_model(s, T, 'vola')
r_T <- axe_narch_predict(model, s_, 'vola')
lines(c(20+T,20+T), c(s_T-r_T, s_T+r_T), col='orange')
}
}
\author{
  Dr. Lars Metzner
}
\encoding{UTF-8}
\references{
  Dr. Lars Metzner (2021) \emph{Adäquates Maschinelles Lernen}.
  Independently Published.
}
\keyword{model}
\keyword{regression}
\keyword{nonparametric}
