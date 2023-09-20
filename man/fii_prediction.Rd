% File sisireg/man/fii_prediction.Rd
\name{fii_prediction}
\alias{fii_prediction}
\title{
Factor-wise Influence Indicator (Prediction-fii) for a given ssrmlp model
regarding a given input vector}
\description{
The Prediction-fii depicts the overall influence of the given input factors on the resulting prediction from a trained ssrmlp model. For computation the components of the input vectors a taken separately as input for the model.
}
\usage{
fii_prediction(W, x)
}
\arguments{
  \item{W}{a trained ssrmlp model}
  \item{x}{a matrix of input vectors}
}
\value{
  \item{fii}{array of influence indicators}
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
fii_prediction(W, X)
}

}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
