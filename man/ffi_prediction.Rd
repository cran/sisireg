% File sisireg/man/ffi_prediction.Rd
\name{ffi_prediction}
\alias{ffi_prediction}
\title{
Factor-wise Influence Indicator (Prediction-FFI) for a given ssrmlp model
regarding a given input vector}
\description{
The Prediction-FFI depicts the overall influence of the given input factors on the resulting prediction from a trained ssrmlp model. For computation the components of the input vectors a taken separately as input for the model.
}
\usage{
ffi_prediction(W, x)
}
\arguments{
  \item{W}{a trained ssrmlp model}
  \item{x}{a matrix of input vectors}
}
\value{
  \item{ffi}{array of influence indicators}
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
ffi_prediction(W, X)
}

}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
