% File sisireg/man/ssrmlp_predict.Rd
\name{ssrmlp_predict}
\alias{ssrmlp_predict}
\title{Prediction function for the ssrMLP}
\description{
Calculates the prediction for a given ssrMLP
}
\usage{
ssrmlp_predict(X, W)
}
\arguments{
  \item{X}{matrix of coordinates.}
  \item{W}{the weight matrices from ssrmlp_train method.}
}
\value{
  \item{Yp}{array with predictions.}
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
Yp <- ssrmlp_predict(X, W)
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
