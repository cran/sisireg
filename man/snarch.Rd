% File sisireg/man/snarch.Rd
\name{snarch}
\alias{snarch}
\title{S-NARCH Model}
\description{Calculates the long-, middle- and short-term trends and vola for a financial time series.}
\usage{
snarch(dat)
}
\arguments{
  \item{dat}{financial time series.}
}
\value{data frame
  \item{tr20}{long-term trend}
  \item{vl20}{long-term vola}
  \item{tr10}{middle-term trend}
  \item{vl10}{middle-term vola}
  \item{tr5}{short-term trend}
  \item{vl5}{short-term vola}
}
\examples{
# generate test data
set.seed(1234)
x <- seq(1:250)
dat <- 13000 + cumsum(rnorm(250))
# calculate the S-NARCH model
df <- snarch(dat)
# plot the results
op <- par(mfrow=c(1,3))
plot(x,dat)
lines(x,df$tr20)
lines(x,df$tr20 - df$vl20, lty = 'dotted')
lines(x,df$tr20 + df$vl20, lty = 'dotted')
plot(x,dat)
lines(x,df$tr10)
lines(x,df$tr10 - df$vl10, lty = 'dotted')
lines(x,df$tr10 + df$vl10, lty = 'dotted')
plot(x,dat)
lines(x,df$tr5)
lines(x,df$tr5 - df$vl5, lty = 'dotted')
lines(x,df$tr5 + df$vl5, lty = 'dotted')
par(op)
}
\author{
  Dr. Lars Metzner
}
\encoding{UTF-8}
\references{
  Dr. Lars Metzner (2019) \emph{Finanzmathematische Zeitreihenanalyse}.
  Independently Published.
}
\keyword{model}
\keyword{ts}
\keyword{nonparametric}
