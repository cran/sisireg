% File sisireg/man/ssrnd2.Rd
\name{ssrnd2}
\alias{ssrnd2}
\title{Multi-dimensional SSR model - reworked}
\description{
Calculates the multi-dimensional SSR model
}
\usage{
ssrnd2(X, Y, mdl=NULL, k=NULL, fn=NULL, n_l1=NULL, iter=1000, eps=0.0001)
}
\arguments{
  \item{X}{matrix with n-dimensional coordinates.}
  \item{Y}{array with observations.}
  \item{mdl}{optional: model from previous training.}
  \item{k}{optional: maxumum run length for the model.}
  \item{fn}{optional: quantile for partial sums.}
  \item{n_l1}{optional: subset size for L1 curvature calculation.}
  \item{iter}{optional: number of iterations for the numeric solver.}
  \item{eps}{optional: delta for ending calculation iteration.}
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
mdl <- ssrnd2(X, Y)
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
