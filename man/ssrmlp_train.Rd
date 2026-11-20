% File sisireg/man/ssrmlp_train.Rd
\name{ssrmlp_train}
\alias{ssrmlp_train}
\title{2-layer MLP with partial sum optimization}
\description{
Calculates the weights of a 2-layer MLP with respect to the partial sums 
critereon
}
\usage{
ssrmlp_train(X, Y, std = TRUE, opt = "ps", nl = 2, hl = NULL, W =
                 NULL, k = NULL, fn = NULL, eta = 0.5, accept = 10,
                 maxIter = 1000, alpha = NULL, beta = NULL)
}
\arguments{
  \item{X}{matrix with n-dimensional coordinates.}
  \item{Y}{array with observations.}
  \item{std}{optional: standardizing values if TRUE.}
  \item{opt}{optional: optimizing function ('ps', 'lse', 'simple').}
  \item{nl}{optional: number of hidden layers.}
  \item{hl}{optional: array tupel with number of perceptrons in each layer.}
  \item{W}{optional: previously calculates weights for refining the model.}
  \item{k}{optional: number of neighbors per quadrant.}
  \item{fn}{optional: quantile for partial sums.}
  \item{eta}{optional: constant factor of the gradient algorithm.}
  \item{accept}{optional: number of accepted PS violations.}
  \item{maxIter}{optional: number of iterations for the numeric solver.}
  \item{alpha}{optional: weight parameter for error function.}
  \item{beta}{optional: weight parameter for error function.}
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
