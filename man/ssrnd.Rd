% File sisireg/man/ssrnd.Rd
\name{ssrnd}
\alias{ssrnd}
\title{Multi-dimensional SSR model}
\description{
Calculates the multi-dimensional SSR model
}
\usage{
ssrnd(X, Y, mdl = NULL, k = NULL, fn = NULL, n_l1 = NULL, iter = 1000, eps = 1e-04)
}
\arguments{
  \item{X}{data frame with n-dimensional coordinates.}
  \item{Y}{data frame with observations.}
  \item{mdl}{optional: pre-trained model.}
  \item{k}{optional: maxumum run length for the model.}
  \item{fn}{optional: quantile for partial sums.}
  \item{n_l1}{optional: number of elements in neighborhood for L1-calculation.}
  \item{iter}{optional: number of iterations for the numeric solver.}
  \item{eps}{optional: barrier parameter for optimization.}
}
\value{
  \item{df}{data frame with the relevant model data.}
}
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
