% File sisireg/man/ssrmlp_train.Rd
\name{ssrmlp_train}
\alias{ssrmlp_train}
\title{2-layer MLP with partial sum optimization}
\description{
Calculates the weights of a 2-layer MLP with respect to the partial sums 
critereon
}
\usage{
ssrmlp_train(X, Y, opt='ps', hl = NULL, W = NULL, k=10, fn=4, 
                   eta=0.75, maxIter=1000)
}
\arguments{
  \item{X}{matrix with n-dimensional coordinates.}
  \item{Y}{array with observations.}
  \item{opt}{optional: optimizing function ('l2' or 'ps'.}
  \item{hl}{optional: array tupel with number of perceptrons in each layer.}
  \item{W}{optional: previously calculates weights for refining the model.}
  \item{k}{optional: number of neighbors per quadrant.}
  \item{fn}{optional: quantile for partial sums.}
  \item{eta}{optional: constant factor of the gradient algorithm.}
  \item{maxIter}{optional: number of iterations for the numeric solver.}
}
\value{
  \item{W}{List with weight matrices.}
}
\examples{
\donttest{
# generate data
set.seed(42)
x <- rnorm(300)
y <- rnorm(300)
z <- rnorm(300) + atan2(x, y)
# coordinates
X <- matrix(cbind(x,y), ncol = 2)
Y <- as.double(z)
# scaling
min_ <- apply(X, 2, min)
max_ <- apply(X, 2, max)
X <- t(apply(X, 1, function(x) (x-min_)/(max_-min_)))
Y = (Y - min(Y)) / (max(Y) - min(Y))
# Training
W <- ssrmlp_train(X, Y)
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
