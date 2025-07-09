% File sisireg/man/ssrnd2_predict.Rd
\name{ssrnd2_predict}
\alias{ssrnd2_predict}
\title{Prediction function for the multi-dimensional SSR model - reworked}
\description{
Calculates the prediction for a given multi-dimensional SSR model
}
\usage{
ssrnd2_predict(mdl, xx)
}
\arguments{
  \item{mdl}{data list with previously calculated model.}
  \item{xx}{matrix with coordinates for prediction.}
}
\value{
  \item{z}{array with predictions.}
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
yy <- ssrnd2_predict(mdl, X)
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
