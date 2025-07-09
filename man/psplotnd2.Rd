% File sisireg/man/psplotnd2.Rd
\name{psplotnd2}
\alias{psplotnd2}
\title{Partial Sum Plot for the multidimensional coordinates - reworked}
\description{
plots the partial sum statistic for the general n-dimensional SSR-model
}
\usage{
psplotnd2(koord, dat, mu, text = "Sample", maxint = NULL, resultplot = TRUE, accept = 10)
}
\arguments{
  \item{koord}{data frame with coordinates.}
  \item{dat}{data frame of observations.}
  \item{mu}{list of discrete regression function.}
  \item{text}{optional: title for the plot.}
  \item{maxint}{optional: upper limit for evaluated subset size.}
  \item{resultplot}{optional: switch for creating a plot, default is TRUE.}
  \item{accept}{optional: percentage of allowed deviations in number and signs in partial sums statistics to accept a result.}
}
\value{TRUE if Partial Sum criterion is passed or result is accepted, FALSE otherwise}
\examples{
\donttest{
# generate data
set.seed(42)
x_ <- y_ <- seq(-3, 3, length=20)
xy <- expand.grid(x=x_,y=y_)
x <- xy$x
y <- xy$y
z <- rnorm(400) + atan2(x, y)
# coordinates
X <- matrix(cbind(x, y), ncol = 2)
Y <- as.double(z)
# Training
W <- ssrmlp2_train(X, Y, accept=0) # passed: ps=39.04, secs 24.8, 250 Iterationen
Yp <- ssrmlp_predict(X, W)
psplotnd2(W$X, W$Y, W$Yp, 'ssrMLP2', resultplot = TRUE)
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
