% File sisireg/man/tauM.Rd
\name{tauM}
\alias{tauM}
\title{Trend-based Correlation}
\description{
Calculates the trend-based correlation of two time series based on the trend function (Metzner's Tau)
}
\usage{
tauM(x, y)
}
\arguments{
  \item{x}{time series.}
  \item{y}{time series.}
}
\value{trend-based correlation.}
\examples{
set.seed(1234)
s <- seq(-pi, pi, length.out = 200)
x <- s + rnorm(200)
y <- exp(s) + 5*rnorm(length(s))
op <- par(mfrow=c(1,2))
plot(x)
plot(y)
par(op)

p <- cor(x,y) # 0.5037
t <- cor(x,y, method = 'kendall') # 0.2959
tm <- tauM(x, y) # 0.0867
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
