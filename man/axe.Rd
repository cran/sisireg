% File sisireg/man/axe.Rd
\name{axe}
\alias{axe}
\title{
Data model for the AxE-Model (Axiomatic Econometric Modeling Paradigm)
}
\description{
Calculation of the relevant data for the AxE-model from a financial time series: trend, volatiliy, change in quotes and risk level.
}
\usage{
axe(quotes)
}
\arguments{
  \item{quotes}{financial time series}
}
\value{data frame
  \item{quotes}{the given time series}
  \item{trend5}{5-day trend}
  \item{trend10}{10-day trend}
  \item{trend20}{20-day trend}
  \item{vola5}{5-day volatility}
  \item{vola10}{10-day volatility}
  \item{vola20}{20-day volatility}
  \item{chng5}{5-day price change}
  \item{chng10}{10-day price change}
  \item{chng20}{20-day price change}
  \item{risk5}{5-day risk level}
  \item{risk10}{10-day risk level}
  \item{risk20}{20-day risk level}
}
\examples{
set.seed(1234)
s <- 13000 + cumsum(rnorm(100))
df_axe <- axe(s)
op <- par(mfrow=c(3,1))
plot(s, type = "l")
plot(df_axe$trend5, type = "l")
abline(a = 0, b = 0)
plot(df_axe$vola5, type = "l")
par(op)
}
\author{
Dr. Lars Metzner
}

\references{
  Dr. Lars Metzner (2020) \emph{Trendbasierte Prognostik}.
  Independently Published.
}
\keyword{model}
\keyword{ts}
\keyword{nonparametric}
