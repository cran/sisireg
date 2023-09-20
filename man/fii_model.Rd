% File sisireg/man/fii_model.Rd
\name{fii_model}
\alias{fii_model}
\title{
Factor-wise Influence Indicator (Model-fii) for a given ssrmlp model
}
\description{
The Model-fii depicts the overall influence of the input factors on the resulting trained ssrmlp model. For computation a unit matrix is used to accumulate the weights for each factor separately.
}
\usage{
fii_model(W)
}
\arguments{
  \item{W}{a trained ssrmlp model}
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
fii_model(W)
}

}
\author{
Dr. Lars Metzner
}

\keyword{models}
\keyword{ts}
\keyword{nonparametric}
