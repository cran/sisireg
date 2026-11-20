% File sisireg/man/psplotnd.Rd
\name{psplotnd}
\alias{psplotnd}
\title{Partial Sum Plot for the multidimensional coordinates}
\description{
plots the partial sum statistic for the general n-dimensional SSR-model
}
\usage{
psplotnd(koord, dat, mu, text = "Sample", maxint = NULL, resultplot = TRUE, accept = 10)
}
\arguments{
  \item{koord}{data frame with coordinates.}
  \item{dat}{data frame of observations.}
  \item{mu}{list of discrete regression function.}
  \item{text}{optional: title for the plot.}
  \item{maxint}{optional: maximum of neighborhood elements}
  \item{resultplot}{optional: boolean plotting test result}, 
  \item{accept}{optional: number of accepted PS violations}
}
\value{No explicit return value: a plot is generated}
\examples{
\donttest{
# generate data
set.seed(1234)
x_ <- y_ <- seq(-3, 3, length=20)
xy <- expand.grid(x=x_,y=y_)
x <- xy$x
y <- xy$y
z <- rnorm(400) + atan2(x, y)
X <- matrix(cbind(x, y), ncol = 2)
Y <- as.double(z)
# Training
df_model <- ssrnd(X, Y, k = 4, fn = 8)
# plot partial sum statistic
psplotnd(X, Y, df_model$mu, 'ssrnd')
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
