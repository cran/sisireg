% File sisireg/man/ssrnd2NN.Rd
\name{ssrnd2NN}
\alias{ssrnd2NN}
\title{Near Neighborhood Estimator for Multi-dimensional SSR model}
\description{
Calculates the Near Neighborhood Estimator in a multi-dimensional SSR model
}
\usage{
ssrnd2NN(X, Y, k, weighted=FALSE)
}
\arguments{
  \item{X}{matrix with n-dimensional coordinates.}
  \item{Y}{array with observations.}
  \item{k}{optional: maxumum run length for the model.}
  \item{weighted}{optional: flag if a weighted mean should be used.}
}
\value{
  \item{mdl}{data list with the relevant model data.}
}
\examples{
\donttest{
# generate data
set.seed(1234)
x <- rnorm(400)
y <- rnorm(400)
z <- rnorm(400) + atan2(x, y)
X <- matrix(cbind(x, y), ncol = 2)
Y <- as.double(z)

# Training
mdl <- ssrnd2NN(X, Y, 10)
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
